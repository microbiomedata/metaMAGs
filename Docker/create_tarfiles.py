#!/usr/bin/env python

import sys
from zipfile import ZipFile
import os
import tarfile
import os.path
import json
import shutil
from gff2txt import parse_cog_tigr_cathfunfam_smart_supfam_input_gff_files


mapping = {
        "_cog.gff": ".cog.txt",
        "_proteins.faa": ".faa",
        "_product_names.tsv": ".gene_product.txt",
        "_ko.tsv": ".ko.txt",
        "_gene_phylogeny.tsv": ".phylodist.txt",
        "_ec.tsv": ".ec.txt",
        "_pfam.gff": ".pfam.txt",
        "_tigrfam.gff": ".tigr.txt",
        "_cath_funfam.gff": ".cathfunfam.txt",
        "_smart.gff": ".smart.txt",
        "_supfam.gff": ".supfam.txt",
        "_functional_annotation.gff": ".gff"
        }

def filter_gff(inf, outf, contig_ids):
    of = open(outf, "w")
    cid = ""
    with open(inf) as f:
        for line in f:
            cid=line.rstrip().split()[0]
            if cid in contig_ids:
                of.write(line)
    of.close()


def filter_faa(inf, outf, contig_ids):
    of = open(outf, "w")
    cid = ""
    with open(inf) as f:
        for line in f:
            if line.startswith(">"):
                fid=line[1:].rstrip().split()[0]
                cid = "_".join(fid.split("_")[0:-2])
            if cid in contig_ids:
                of.write(line)
    of.close()

def filter_inp(inf, outf, contig_ids):
    of = open(outf, "w")
    with open(inf) as f:
        for line in f:
            fid=line.rstrip().split()[0]
            cid = "_".join(fid.split("_")[0:-2])
            if cid in contig_ids:
                of.write(line)
    of.close()

def filter_txt(inf, outf, contig_ids):
    tout = "tmp.gff"
    filter_inp(inf, tout, contig_ids)
    parse_cog_tigr_cathfunfam_smart_supfam_input_gff_files(tout, outf, "bogus", "bogus")
    os.unlink(tout)

def find_ext(inf):
    found = False
    for pattern, ext in mapping.items():
        if pattern in inf:
            return ext
    return None

def create_new_zip(prefix, binid, binf, inputs, contig_ids):
    outd = "%s_%s" % (prefix, binid)
    if not os.path.exists(outd):
        os.mkdir(outd)
    # Create the contigs file
    outfn = "%s_%s%s" % (prefix, binid, ".fna")
    shutil.copy(binf, os.path.join(outd, outfn))
    outs = [outfn]
    for inf in inputs:
        # Find the input type
        ext = find_ext(inf)
        if not ext:
            sys.stderr.write("Unknown type: %s\n" % (inf))
            continue
        outfn = "%s_%s%s" % (prefix, binid, ext)
        outf = "%s/%s" % (outd, outfn)
        outs.append(outfn)
        if ext.endswith(".faa"):
            filter_faa(inf, outf, contig_ids)
        elif ext.endswith(".gff"):
            filter_gff(inf, outf, contig_ids)
        elif inf.endswith(".gff") and ext.endswith(".txt"):
            filter_txt(inf, outf, contig_ids)
        else:
            filter_inp(inf, outf, contig_ids)

    tgzf = "%s_%s.tar.gz" % (prefix, binid)
    with tarfile.open(tgzf, "w:gz") as tar:
        for ofn in outs:
            tar.add("%s/%s" % (outd, ofn), arcname=ofn)


if __name__ == "__main__":
    data = None
    inputs = []
    bin_files = dict()
    # Extract files from inputs
    for f in sys.argv[1:]:
        fn = os.path.basename(f)
        if fn.endswith(".json"):
            data = json.load(open(f))
        elif fn.startswith("bins") and fn.endswith(".fa"):
            binid = fn.replace('.fa', '')
            bin_files[binid] = f
        else:
            inputs.append(f)
    prefix = inputs[0].split('/')[-1].replace('_proteins.faa','')
    for d in data['mags_list']:
        if d['bin_quality'] in ['MQ', 'HQ']:
            print("Processing {bin_name}".format(**d))
            bin_id = d['bin_name']
            contigs = d['members_id']
            bin_file = bin_files[bin_id]
            create_new_zip(prefix, bin_id, bin_file, inputs, contigs)
