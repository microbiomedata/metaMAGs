#!/usr/bin/env python3
import argparse
import collections
import glob
import json
import os
import sqlite3
import subprocess
import sys
import time
from datetime import datetime


def line_count(fname):
    count = len(open(fname).readlines(  ))
    return count

def count_fasta(fname):
    """
    Count total contigs, and contigs less than 3000 bp.
    Returns tuple: (total, less_than_3000)
    """
    total = 0
    short = 0
    length = 0
    with open(fname, "r") as f:
        for line in f:
            if line.startswith(">"):
                if total > 0 and length < 3000:
                    short += 1
                length = 0
                total += 1
            else:
                length += len(line.strip())
        # Final contig at EOF
        if total > 0 and length < 3000:
            short += 1
    return total, short

def write_mags_tsv(mag_list, outtsvfile, prok_bin_method, euk_bin_method):
    """
    Write MAGs_stats.tsv from a list of MAG metadata dicts (mag_list).
    Each element of mag_list is the OrderedDict created in mag_meta().
    """
    with open(outtsvfile, "w") as f:
        # Header comments
        f.write("#Prok bin_methods: {}\n".format(prok_bin_method))
        f.write("#Euk bin_methods: {}\n".format(euk_bin_method))
        # Column header
        f.write(
            "\t".join(
                [
                    "Bin ID",
                    "Bin Quality",
                    "GTDB-TK lineage",
                    "Bin Completeness",
                    "Bin Contamintation",
                    "Total Number of Bases",
                    "Number of genes",
                    "Num of 5s rRNA",
                    "Num of 16s rRNA",
                    "Num of 23s rRNA",
                    "Num of tRNA",
                    "euk_ncbi_lineage_tax_ids",
                    "euk_completeness",
                    "euk_contamination",
                    "scaffold IDs of members",
                ]
            )
            + "\n"
        )
        
        # Data rows: use values from the dict rather than the original SQL row
        for mag in mag_list:
            gtdb_lineage = ";".join(
                [
                    str(mag.get("gtdbtk_domain", "")),
                    str(mag.get("gtdbtk_phylum", "")),
                    str(mag.get("gtdbtk_class", "")),
                    str(mag.get("gtdbtk_order", "")),
                    str(mag.get("gtdbtk_family", "")),
                    str(mag.get("gtdbtk_genus", "")),
                ]
            )
            # Default empty values for new fields
            ncbi_lineage_tax_ids = ""
            euk_completeness = ""
            euk_contamination = ""
            if mag.get("bin_quality", "") == "LQ" and "eukaryotic_evaluation" in mag:
                euk_eval = mag["eukaryotic_evaluation"]
                ncbi_lineage_tax_ids = str(euk_eval.get("ncbi_lineage_tax_ids", ""))
                euk_completeness = (
                    str(euk_eval.get("completeness", "")) if "completeness" in euk_eval else ""
                )
                euk_contamination = (
                    str(euk_eval.get("contamination", "")) if "contamination" in euk_eval else ""
                )
            members = ",".join(mag.get("members_id", []))

            f.write(
                "\t".join(
                    [
                        str(mag.get("bin_name", "")),
                        str(mag.get("bin_quality", "")),
                        gtdb_lineage,
                        str(mag.get("completeness", "")),
                        str(mag.get("contamination", "")),
                        str(mag.get("total_bases", "")),
                        str(mag.get("gene_count", "")),
                        str(mag.get("num_5s", "")),
                        str(mag.get("num_16s", "")),
                        str(mag.get("num_23s", "")),
                        str(mag.get("num_t_rna", "")),
                        ncbi_lineage_tax_ids,
                        euk_completeness,
                        euk_contamination,
                        members,
                    ]
                )
                + "\n"
            )

def mag_meta(dbname):
    conn = sqlite3.connect(dbname)
    c = conn.cursor()

    qsql = "select bin_name, asm_scaffold_id, scaffold_id from bin_scaffolds;"
    rvals = c.execute(qsql)
    d = collections.defaultdict()
    d_list = collections.defaultdict(list)
    total_bin_contig = 0
    for row in rvals:
        d_list[row[0]].append(row[2])
        d[row[0]] = d[row[0]] + 1 if row[0] in d else 1
        total_bin_contig += 1

    # Get bin_method from bin table with LIMIT 1
    c.execute(
        "SELECT bin_methods FROM bin WHERE bin_quality IN ('MQ', 'HQ') LIMIT 1;"
    )
    prok_bin_method_row = c.fetchone()
    prok_bin_method = prok_bin_method_row[0] if prok_bin_method_row else "N/A"

    c.execute("SELECT bin_methods FROM euk_bin limit 1;")
    euk_bin_method_row = c.fetchone()
    euk_bin_method = euk_bin_method_row[0] if euk_bin_method_row else "N/A"

    qsql = """
        select * from bin
        order by case
            when bin_quality = 'HQ' then 1
            when bin_quality = 'MQ' then 2
            when bin_quality = 'LQ' then 3
        end,
        completeness DESC;
    """
    rvals = c.execute(qsql)
    out_list = []

    for row in rvals:
        tmp_d = collections.OrderedDict()
        tmp_d["bin_name"] = row[0]
        tmp_d["number_of_contig"] = d[row[0]]
        tmp_d["completeness"] = row[8]
        tmp_d["contamination"] = row[9]
        tmp_d["total_bases"] = row[10]
        tmp_d["gene_count"] = row[11]
        tmp_d["bin_quality"] = row[12]
        tmp_d["num_16s"] = row[13]
        tmp_d["num_5s"] = row[14]
        tmp_d["num_23s"] = row[15]
        tmp_d["num_t_rna"] = row[16]
        tmp_d["gtdbtk_domain"] = row[20]
        tmp_d["gtdbtk_phylum"] = row[21]
        tmp_d["gtdbtk_class"] = row[22]
        tmp_d["gtdbtk_order"] = row[23]
        tmp_d["gtdbtk_family"] = row[24]
        tmp_d["gtdbtk_genus"] = row[25]
        tmp_d["gtdbtk_species"] = row[26]
        tmp_d["members_id"] = d_list[row[0]]
        out_list.append(tmp_d)

    conn.close()
    # return the MAG dicts plus the methods/summary info needed for TSV
    return out_list, total_bin_contig, prok_bin_method, euk_bin_method


def parse_eukcc(dbname, bin_list):
    conn = sqlite3.connect(dbname)
    c = conn.cursor()
    #qsql="select bin_name, count(scaffold_id) from bin_scaffolds group by bin_name;"
    qsql="select bin_name, completeness, contamination, ncbi_lng, bin_methods from euk_bin;"
    rvals = c.execute(qsql)
    eukcc=collections.defaultdict(dict)
    for row in rvals:
        eukcc[row[0]]['completeness'] = float(row[1])
        eukcc[row[0]]['contamination'] = float(row[2])
        eukcc[row[0]]['ncbi_lineage_tax_ids'] = row[3]
        eukcc[row[0]]['bin_methods'] = row[4]
    update_bin_list=[]
    
    for bin in bin_list:
        if bin['bin_name'] in eukcc:
            if 'eukaryotic_evaluation' not in bin:
                bin['eukaryotic_evaluation']= dict()
            bin['eukaryotic_evaluation']['completeness'] = float(eukcc[bin['bin_name']]['completeness'])
            bin['eukaryotic_evaluation']['contamination'] = float(eukcc[bin['bin_name']]['contamination'])
            bin['eukaryotic_evaluation']['ncbi_lineage_tax_ids'] = eukcc[bin['bin_name']]['ncbi_lineage_tax_ids']
            bin['eukaryotic_evaluation']['bin_methods'] = eukcc[bin['bin_name']]['bin_methods']
        update_bin_list.append(bin)
    conn.close()
    return update_bin_list

def main():
    parser = argparse.ArgumentParser(description="MAGs meta/stat collector")
    parser.add_argument("--sdb", required=True, help="Input .sdb file")
    parser.add_argument("--fasta", required=True, help="Input FASTA file")
    parser.add_argument("--outdir", default=os.getcwd(), help="Output directory")
    args = parser.parse_args()

    outdir = args.outdir
    outfile = os.path.join(outdir, "MAGs_stats.json")
    outtsvfile = os.path.join(outdir, "MAGs_stats.tsv")

    metadata = collections.OrderedDict()

    # Count FASTA contigs
    total_contigs, short_contigs = count_fasta(args.fasta)
    metadata['input_contig_num'] = total_contigs
    metadata['too_short_contig_num'] = short_contigs

    # Extract MAG stats from .sdb (no TSV here)
    mag_list, total_bin_contig_num, prok_bin_method, euk_bin_method = mag_meta(
        args.sdb
    )
    mag_list = parse_eukcc(args.sdb, mag_list)

    metadata['binned_contig_num'] = total_bin_contig_num
    metadata['unbinned_contig_num'] = total_contigs - total_bin_contig_num - short_contigs

    metadata['mags_list'] = mag_list
    with open(outfile, 'w') as outjson:
        json.dump(metadata, outjson, indent=4)
    # Write TSV from the same JSON-style dicts
    write_mags_tsv(mag_list, outtsvfile, prok_bin_method, euk_bin_method)


if __name__ == "__main__":
    main()
