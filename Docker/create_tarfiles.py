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
from zipfile import ZipFile
from gff2txt import parse_cog_tigr_cathfunfam_smart_supfam_input_gff_files

# File extension mapping
EXTENSION_MAPPING = {
    "cog.gff": ".cog.txt",
    "proteins.faa": ".faa",
    "product_names.tsv": ".gene_product.txt",
    "products.tsv": ".gene_product.txt",
    "ko.tsv": ".ko.txt",
    "gene_phylogeny.tsv": ".phylodist.txt",
    "ec.tsv": ".ec.txt",
    "pfam.gff": ".pfam.txt",
    "tigrfam.gff": ".tigr.txt",
    "cath_funfam.gff": ".cathfunfam.txt",
    "smart.gff": ".smart.txt",
    "supfam.gff": ".supfam.txt",
    "functional_annotation.gff": ".gff"
}

# Filter functions
def filter_gff(input_file, output_file, contig_ids):
    with open(output_file, "w") as out_file:
        with open(input_file) as in_file:
            for line in in_file:
                contig_id = line.rstrip().split()[0]
                if contig_id in contig_ids:
                    out_file.write(line)

def filter_faa(input_file, output_file, contig_ids):
    with open(output_file, "w") as out_file:
        with open(input_file) as in_file:
            for line in in_file:
                if line.startswith(">"):
                    file_id = line[1:].rstrip().split()[0]
                    contig_id = "_".join(file_id.split("_")[0:-2])
                if contig_id in contig_ids:
                    out_file.write(line)

def filter_inp(input_file, output_file, contig_ids):
    with open(output_file, "w") as out_file:
        with open(input_file) as in_file:
            for line in in_file:
                file_id = line.rstrip().split()[0]
                contig_id = "_".join(file_id.split("_")[0:-2])
                if contig_id in contig_ids:
                    out_file.write(line)

def filter_txt(input_file, output_file, contig_ids):
    temp_out = "tmp.gff"
    filter_inp(input_file, temp_out, contig_ids)
    parse_cog_tigr_cathfunfam_smart_supfam_input_gff_files(temp_out, output_file, "bogus", "bogus")
    os.unlink(temp_out)

def find_extension(input_file):
    for pattern, extension in EXTENSION_MAPPING.items():
        if pattern in input_file:
            return extension
    return None

def get_ko_list(input_file,output_file,contig_ids):
    with open(output_file, "w") as out_file:
        with open(input_file) as in_file:
            for line in in_file:
                ko_id = line.rstrip().split()[2]
                file_id = line.rstrip().split()[0]
                contig_id = "_".join(file_id.split("_")[0:-2])
                if contig_id in contig_ids:
                    out_file.write(ko_id.replace("KO:","") + "\n")

def get_bin_annotations(prefix, bin_id, bin_file, inputs, contig_ids, output_dir):
    
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    # Create the contigs file
    output_filename = f"{prefix}_{bin_id}.fna"
    shutil.copy(bin_file, os.path.join(output_dir, output_filename))
    output_files = [output_filename]

    for input_file in inputs:
        # Find the input type
        extension = find_extension(input_file)
        if not extension:
            sys.stderr.write(f"Unknown type: {input_file}\n")
            continue
        
        output_filename = f"{prefix}_{bin_id}{extension}"
        output_file = f"{output_dir}/{output_filename}"
        output_files.append(output_filename)
        if extension.endswith(".faa"):
            filter_faa(input_file, output_file, contig_ids)
        elif extension.endswith(".gff"):
            filter_gff(input_file, output_file, contig_ids)
        elif input_file.endswith(".gff") and extension.endswith(".txt"):
            filter_txt(input_file, output_file, contig_ids)
        elif input_file.endswith("ko.tsv"):
            get_ko_list(input_file, f"{bin_id}.ko", contig_ids)
            filter_inp(input_file, output_file, contig_ids)
        else:
            filter_inp(input_file, output_file, contig_ids)

def ko_analysis(prefix):
    ko_list = glob.glob("*.ko")
    if ko_list:
        cmd = ["ko_mapper.py" , "-i"] +  sorted(ko_list) + [ "-p" , prefix] 
        proc = subprocess.Popen(shlex.split(" ".join(cmd)), shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        outs, errs = proc.communicate()
        if proc.returncode == 0:
            return f"{prefix}_module_completeness.tab"
        else:
            print(errs.decode().rstrip())
            return
    
def krona_plot(ko_result,prefix):
    ko_list = glob.glob("*.ko")
    if ko_list:
        df = pd.read_csv(ko_result,sep="\t")
        krona_list=[]
        for ko in ko_list:
            krona_list.append(f"{ko}.krona")
            df[[ko,'pathway group','name']].to_csv(f"{ko}.krona",sep="\t",index=False,header=False)
        cmd=["ktImportText"] + krona_list + [ "-o" , f"{prefix}_ko_krona.html"]
        proc = subprocess.Popen(shlex.split(" ".join(cmd)), shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        outs, errs = proc.communicate()
        if proc.returncode == 0:
            return f"{prefix}_ko_krona.html"
        else:
            print(errs.decode().rstrip())
            return
    
def create_new_zip(bin_dirs):
    for dir in bin_dirs:
        tar_file_name = f"{dir}.tar.gz"
        with tarfile.open(tar_file_name, "w:gz") as tar:
            for output_file_name in glob.glob(f"{dir}/*",recursive=True):
                tar.add(f"{output_file_name}", arcname=output_file_name)
		
if __name__ == "__main__":
    data = None
    input_files = []
    bin_files_dict = {}
    bin_dirs = []
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
            print(f"Processing {bin_data['bin_name']}")
            bin_id = bin_data['bin_name']
            contig_ids = bin_data['members_id']
            bin_file = bin_files_dict[bin_id]
            output_dir = f"{prefix}_{bin_id}_{bin_data['bin_quality']}"
            bin_dirs.append(output_dir)
            get_bin_annotations(prefix, bin_id, bin_file, input_files, contig_ids, output_dir)
    
    ko_result = ko_analysis(prefix)
    krona_plot(ko_result,prefix)
    create_new_zip(bin_dirs)
