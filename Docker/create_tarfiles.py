#!/usr/bin/env python3
import sys
import os
import json
import shutil
import tarfile
import glob
import subprocess
import shlex
import pandas as pd
from pathlib import Path
import fitz
from gff2txt import parse_cog_tigr_cathfunfam_smart_supfam_input_gff_files
from multiprocessing import Pool
from time import time
from collections import OrderedDict



__version__ = "0.7.0"


# File extension mapping
EXTENSION_MAPPING = {
    "cog.gff": ".cog.txt",
    "proteins.faa": ".faa",
    "product_names.tsv": ".gene_product.txt",
    "products.tsv": ".gene_product.txt",
    "ko.tsv": ".ko.txt",
    "gene_phylogeny.tsv": ".phylodist.txt",
    "ec.tsv": ".ec.txt",
    "crispr.tsv": ".crisprs.txt",
    "pfam.gff": ".pfam.txt",
    "tigrfam.gff": ".tigr.txt",
    "cath_funfam.gff": ".cathfunfam.txt",
    "smart.gff": ".smart.txt",
    "supfam.gff": ".supfam.txt",
    "functional_annotation.gff": ".gff"
}


def get_contig_faa(line, contig_id):
    if line.startswith(">"):
        file_id = line[1:].rstrip().split()[0]
        contig_id = "_".join(file_id.split("_")[0:-2])
    return contig_id


def get_contig_gff(line, contig_id):
    return line.rstrip().split()[0]


def get_contig_tsv(line, contig_id):
    file_id = line.rstrip().split()[0]
    return "_".join(file_id.split("_")[0:-2])

def filter_one_pass(input_file, prefix, mags_data, ext, filter_func,
                    post=None):
    """
    This function reads an input and splits it out into files for each
    MAG.

    Inputs:
    input_file: input filename
    prefix: prefix to use for each output file
    mags_data: list of dictionary objects with data for each MAG
    ext: extension to use
    filter_func: Function that will be called on each line to extract
                 the contig_id
    post: optional function to call at the end given a list of output_files
    """
    contig_fd_map = dict()
    fd_list = []
    output_files = []
    # This will open an output file for each bin and then
    # create a map for the contigs to these open files
    for mag in mags_data:
        bin_id = mag['bin_name']
        contig_ids = mag['members_id']
        output_dir = mag['output_dir']
        output_filename = f"{prefix}_{bin_id}{ext}"
        output_file = f"{output_dir}/{output_filename}"
        output_files.append(output_file)
        fd = open(output_file, 'w')
        fd_list.append(fd)
        for contig_id in contig_ids:
            if contig_id not in contig_fd_map:
                contig_fd_map[contig_id] = []
            contig_fd_map[contig_id].append(fd)
    contig_id = None
    # Now we go through the input file and demultiplex
    # the output.
    with open(input_file) as in_file:
        for line in in_file:
            # Extract the contig_id using the provided
            # function.
            contig_id = filter_func(line, contig_id)
            if contig_id in contig_fd_map:
                for fd in contig_fd_map[contig_id]:
                    fd.write(line)
    # Cleanup
    for fd in fd_list:
        fd.close()

    # Call optional post function
    if post:
        post(output_files)


def find_extension(input_file):
    for pattern, extension in EXTENSION_MAPPING.items():
        if pattern in input_file:
            return extension
    return None


def parse_gffs(output_files):
    """
    post function to run on gff->txt outputs
    """
    for temp_out in output_files:
        output_file = temp_out.replace(".tmp", "")
        parse_cog_tigr_cathfunfam_smart_supfam_input_gff_files(temp_out,
                                                               output_file,
                                                               "bogus",
                                                               "bogus")
        os.unlink(temp_out)


def write_kos(output_files):
    """
    Post function to extract the list of KOs
    """
    for in_file in output_files:
        bin_id = in_file.split(".")[-3]
        output_file = f"bins.{bin_id}.ko"
        write_ko_list(in_file, output_file)


def write_ko_list(input_file, output_file):
    with open(output_file, "w") as out_file:
        with open(input_file) as in_file:
            for line in in_file:
                ko_id = line.rstrip().split()[2]
                out_file.write(ko_id.replace("KO:", "") + "\n")


def rewrite_files(prefix, inputs, mags):
    for input_file in inputs:
        # Find the input type
        extension = find_extension(input_file)
        if not extension:
            sys.stderr.write(f"Unknown type: {input_file}\n")
            continue

        start = time()
        post = None
        if input_file.endswith(".faa"):
            filter_func = get_contig_faa
        elif input_file.endswith(".gff") and extension.endswith(".txt"):
            filter_func = get_contig_tsv
            extension += ".tmp"
            post = parse_gffs
        elif extension.endswith(".gff"):
            filter_func = get_contig_gff
        elif input_file.endswith("ko.tsv"):
            filter_func = get_contig_tsv
            post = write_kos
        elif input_file.endswith("crispr.tsv"):
            filter_func = get_contig_gff
        else:
            filter_func = get_contig_tsv
        filter_one_pass(input_file, prefix, mags, extension, filter_func,
                        post=post)
        
        print(f" - {input_file.split('/')[-1]}: {time()-start:.3f}s")   


def ko_analysis(prefix):
    ko_list = glob.glob("*.ko")
    if ko_list:
        cmd = ["ko_mapper.py", "-i"] + sorted(ko_list) + ["-p", prefix]
        proc = subprocess.Popen(shlex.split(" ".join(cmd)), shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        outs, errs = proc.communicate()
        if proc.returncode == 0:
            if os.path.exists(f"{prefix}_barplot.pdf"):
                pdf_to_png(f"{prefix}_barplot.pdf")
            if os.path.exists(f"{prefix}_heatmap.pdf"):
                pdf_to_png(f"{prefix}_heatmap.pdf")
            return f"{prefix}_module_completeness.tab"
        else:
            print(errs.decode().rstrip(), file=sys.stderr)
            sys.exit()


def pdf_to_png(pdf):
    prefix = Path(pdf).stem
    doc = fitz.open(pdf)  # open document
    mat = fitz.Matrix(2, 2)   # zoom factor 2 in each dimension
    for page in doc:  # iterate through the pages
        pix = page.get_pixmap(matrix=mat)  # render page to an image
        pix.save(f"{prefix}.png")  # store image as a PNG
    return True


def krona_plot(ko_result, prefix):
    ko_list = glob.glob("*.ko")
    if ko_list:
        df = pd.read_csv(ko_result, sep="\t")
        krona_list = []
        for ko in ko_list:
            krona_list.append(f"{ko}.krona")
            df[[ko, 'pathway group', 'name']].to_csv(f"{ko}.krona", sep="\t",
                                                     index=False,
                                                     header=False)
        cmd = ["ktImportText"] + krona_list + ["-o", f"{prefix}_ko_krona.html"]
        proc = subprocess.Popen(shlex.split(" ".join(cmd)), shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        outs, errs = proc.communicate()
        if proc.returncode == 0:
            return f"{prefix}_ko_krona.html"
        else:
            print(errs.decode().rstrip(), file=sys.stderr)
            return

def gene_count(bin_dirs):
    mags_list=[]
    for bin_dir in bin_dirs:
        bin_data=bin_dirs[bin_dir]
        for output_file_name in glob.glob(f"{bin_dir}/*", recursive=True):
            if output_file_name.endswith(".fna"):
                total_bases=0
                with open (output_file_name, "r") as f :
                    for line in f:
                        if(line[0] == ">") :
                            continue
                        total_bases += len(line)-1
                bin_data['total_bases'] = total_bases
            if output_file_name.endswith(".gff"):
                with open(output_file_name, 'r') as fp:
                    lines = len(fp.readlines())
                bin_data['gene_count'] = lines
        mags_list.append(bin_data)
    return mags_list

def create_tar_file(bin_dir):
    tar_file_name = f"{bin_dir}.tar.gz"
    with tarfile.open(tar_file_name, "w:gz") as tar:
        for output_file_name in glob.glob(f"{bin_dir}/*", recursive=True):
            tar.add(f"{output_file_name}", arcname=output_file_name)


def create_tarfiles(bin_dirs, threads):
    """
    This parallelize the creation of the tar files
    """
    def error_cb(e):
        raise e

    p = Pool(threads)
    for bin_dir in bin_dirs:
        p.apply_async(create_tar_file, args=(bin_dir, ),
                      error_callback=error_cb)
    p.close()
    p.join()


if __name__ == "__main__":
    data = None
    input_files = []
    bin_files_dict = {}
    bin_dirs = OrderedDict()
    threads = int(os.environ.get("THREADS", "32"))
    prefix = sys.argv[1]
    for file in sys.argv[2:]:
        file_name = os.path.basename(file)
        if file_name.endswith(".json"):
            data = json.load(open(file))
        elif file_name.startswith("bins") and file_name.endswith(".fa"):
            bin_id = file_name.replace('.fa', '')
            bin_files_dict[bin_id] = file
        else:
            input_files.append(file)
    for bin_data in data['mags_list']:
        if bin_data['bin_quality'] in ['MQ', 'HQ', 'LQ']:
            bin_id = bin_data['bin_name']
            contig_ids = bin_data['members_id']
            bin_file = bin_files_dict[bin_id]
            output_dir = f"{prefix}_{bin_id}_{bin_data['bin_quality']}"
            if not os.path.exists(output_dir):
                os.mkdir(output_dir)
            bin_data['output_dir'] = output_dir
            bin_dirs[output_dir] = bin_data
            output_filename = f"{prefix}_{bin_id}.fna"
            shutil.copy(bin_file, os.path.join(output_dir, output_filename))
    print(f"Processing {len(bin_dirs)} mags")
    rewrite_files(prefix, input_files, data['mags_list'])
    print("Generating KOs")
    ko_result = ko_analysis(prefix)
    print("Generating Krona Plot")
    krona_plot(ko_result, prefix)
    print("Count Total base and Gene for bins")
    mags_list = gene_count(bin_dirs)
    print(f"Update {prefix}_stats.json")
    data['mags_list'] = mags_list
    with open(f"{prefix}_stats.json", "w") as of:
        json.dump(data,of,indent = 4)

    print("Generating zip")
    create_tarfiles(bin_dirs, threads)
