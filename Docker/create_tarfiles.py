import sys
import os
import json
import shutil
import tarfile
from zipfile import ZipFile
from gff2txt import parse_cog_tigr_cathfunfam_smart_supfam_input_gff_files

# File extension mapping
EXTENSION_MAPPING = {
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

def create_new_zip(prefix, bin_id, bin_file, inputs, contig_ids):
    output_dir = f"{prefix}_{bin_id}"
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
        
        output_file = f"{output_dir}/{output_filename}"
        output_files.append(output_filename)
        if extension.endswith(".faa"):
            filter_faa(input_file, output_file, contig_ids)
        elif extension.endswith(".gff"):
            filter_gff(input_file, output_file, contig_ids)
        elif input_file.endswith(".gff") and extension.endswith(".txt"):
            filter_txt(input_file, output_file, contig_ids)
        else:
            filter_inp(input_file, output_file, contig_ids)

    tar_file_name = f"{prefix}_{bin_id}.tar.gz"
    with tarfile.open(tar_file_name, "w:gz") as tar:
        for output_file_name in output_files:
            tar.add(f"{output_dir}/{output_file_name}", arcname=output_file_name)


if __name__ == "__main__":
    data = None
    input_files = []
    bin_files_dict = {}
    for file in sys.argv[1:]:
        file_name = os.path.basename(file)
        if file_name.endswith(".json"):
            data = json.load(open(file))
        elif file_name.startswith("bins") and file_name.endswith(".fa"):
            bin_id = file_name.replace('.fa', '')
            bin_files_dict[bin_id] = file
        else:
            input_files.append(file)
    prefix = input_files[0].split('/')[-1].replace('_proteins.faa','')
    for bin_data in data['mags_list']:
        if bin_data['bin_quality'] in ['MQ', 'HQ']:
            print(f"Processing {bin_data['bin_name']}")
            bin_id = bin_data['bin_name']
            contig_ids = bin_data['members_id']
            bin_file = bin_files_dict[bin_id]
            create_new_zip(prefix, bin_id, bin_file, input_files, contig_ids)