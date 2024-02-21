#  Metagenome assembled genomes generation workflow

## Summary

The workflow is based on [IMG MAGs pipeline](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6323987/)<sup>1</sup> for metagenome assembled genomes generation. It takes assembled contigs, bam file from reads mapping to contigs and [contigs annotations](https://github.com/microbiomedata/mg_annotation) result to to associate groups of contigs as deriving from a seemingly coherent microbial species (binning) and evaluted by checkM and gtdb-tk. 

## Required Database

* [CheckM](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4484387/)<sup>2</sup> database is 275MB contains the databases used for the Metagenome Binned contig quality assessment. (requires 40GB+ of memory, included in the image) 
    * https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz

* [GTDB-Tk](https://doi.org/10.1093/bioinformatics/btz848)<sup>3</sup> requires ~78G of external data that need to be downloaded and unarchived. (requires ~150GB of memory)
    * https://data.gtdb.ecogenomic.org/releases/release214/214.0/auxillary_files/gtdbtk_r214_data.tar.gz

* Prepare the Database

```bash
    wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
    tar -xvzf checkm_data_2015_01_16.tar.gz
    mkdir -p refdata/CheckM_DB && tar -xvzf checkm_data_2015_01_16.tar.gz -C refdata/CheckM_DB
    rm checkm_data_2015_01_16.tar.gz
```

```bash
    wget https://data.gtdb.ecogenomic.org/releases/release214/214.0/auxillary_files/gtdbtk_r214_data.tar.gz
    mkdir -p refdata/GTDBTK_DB && tar -xvzf gtdbtk_r214_data.tar.gz 
    mv release214 refdata/GTDBTK_DB
    rm gtdbtk_r214_data.tar.gz
```

## Running Workflow in Cromwell

Description of the files:
 - `.wdl` file: the WDL file for workflow definition
 - `.json` file: the example input for the workflow
 - `.conf` file: the conf file for running Cromwell.
 - `.sh` file: the shell script for running the example workflow (sbatch)

## The Docker image 

[microbiomedata/nmdc_mbin](https://hub.docker.com/r/microbiomedata/nmdc_mbin)

## Input files

A json files with following entries:

1. Project Name
2. Metagenome Assembled Contig fasta file
3. Sam/Bam file from reads mapping back to contigs.
4. Contigs functional annotation result in gff format
5. Contigs functional annotated protein FASTA file
6. Tab delimited file for [COG](http://reusabledata.org/cogs) annotation.
7. Tab delimited file for [EC](https://reusabledata.org/kegg-ftp) annotation.
8. Tab delimited file for [KO](https://reusabledata.org/kegg-ftp) annotation.
9. Tab delimited file for [PFAM](http://reusabledata.org/pfam) annotation.
10. Tab delimited file for [TIGRFAM](http://reusabledata.org/tigrfams) annotation.
11. Tab delimited file for [CATH FUNFAM](http://reusabledata.org/cath) annotation.
12. Tab delimited file for [SMART](https://reusabledata.org/smart) annotation.
13. Tab delimited file for [SUPER FAMILY](https://reusabledata.org/supfam) annotation.
14. Tab delimited file for Gene Product name assignment.
15. Tab delimited file for Gene Phylogeny assignment.
16. Tab delimited file for Contig/Scaffold lineage.
17. GTDBTK Database
18. CheckM Database
19. (optional) nmdc_mags.threads: The number of threads used by metabat/samtools/checkm/gtdbtk. default: 64
20. (optional) nmdc_mags.pthreads: The number of threads used by pplacer (Use lower number to reduce the memory usage) default: 1

```
{
    "nmdc_mags.proj_name": "nmdc_wfmgan-xx-xxxxxxxx",
    "nmdc_mags.contig_file": "/path/to/Assembly/nmdc_wfmgas-xx-xxxxxxx_contigs.fna",
    "nmdc_mags.sam_file": "/path/to/Assembly/nmdc_wfmgas-xx-xxxxxxx_pairedMapped_sorted.bam",
    "nmdc_mags.gff_file": "/path/to/Annotation/nmdc_wfmgas-xx-xxxxxxx_functional_annotation.gff",
    "nmdc_mags.proteins_file": "/path/to/Annotation/nmdc_wfmgas-xx-xxxxxxx_proteins.faa",
    "nmdc_mags.cog_file": "/path/to/Annotation/nmdc_wfmgas-xx-xxxxxxx_cog.gff",
    "nmdc_mags.ec_file": "/path/to/Annotation/nmdc_wfmgas-xx-xxxxxxx_ec.tsv",
    "nmdc_mags.ko_file": "/path/to/Annotation/nmdc_wfmgas-xx-xxxxxxx_ko.tsv",
    "nmdc_mags.pfam_file": "/path/to/Annotation/nmdc_wfmgas-xx-xxxxxxx_pfam.gff",
    "nmdc_mags.tigrfam_file": "/path/to/Annotation/nmdc_wfmgas-xx-xxxxxxxtigrfam.gff",
    "nmdc_mags.cath_funfam_file": "/path/to/Annotation/nmdc_wfmgas-xx-xxxxxxx_cath_funfam.gff",
    "nmdc_mags.smart_file": "/path/to/Annotation/nmdc_wfmgas-xx-xxxxxxx_smart.gff",
    "nmdc_mags.supfam_file": "/path/to/Annotation/nmdc_wfmgas-xx-xxxxxxx.gff",
    "nmdc_mags.product_names_file": "/path/to/Annotation/nmdc_wfmgas-xx-xxxxxxx_product_names.tsv",
    "nmdc_mags.gene_phylogeny_file": "/path/to/Annotation/nmdc_wfmgas-xx-xxxxxxx_gene_phylogeny.tsv",
    "nmdc_mags.lineage_file": "/path/to/Annotation/nmdc_wfmgas-xx-xxxxxxx_scaffold_lineage.tsv",
    "nmdc_mags.gtdbtk_db": "refdata/GTDBTK_DB",
    "nmdc_mags.checkm_db": "refdata/CheckM_DB"
}
```

## Output files

The output will have a bunch of output directories, files, including statistical numbers, status log and a shell script to reproduce the steps etc. 

The final [MiMAG](https://www.nature.com/articles/nbt.3893#Tab1) output includes following files.

```
|-- project_name_mags_stats.json
|-- project_name_hqmq_bin.zip
|-- project_name_lq_bin.zip
|-- project_name_bin.info
|-- project_name_bins.lowDepth.fa
|-- project_name_bins.tooShort.fa
|-- project_name_bins.unbinned.fa
|-- project_name_checkm_qa.out
|-- project_name_gtdbtk.ar122.summary.tsv
|-- project_name_gtdbtk.bac122.summary.tsv
|-- project_name_heatmap.pdf
|-- project_name_barplot.pdf
|-- project_name_kronaplot.html
```

### Citation
1. Chen IA, Chu K, Palaniappan K, et al. IMG/M v.5.0: an integrated data management and comparative analysis system for microbial genomes and microbiomes. Nucleic Acids Res. 2019;47(D1):D666‐D677. [doi:10.1093/nar/gky901](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6323987/)
2. Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome Res. 2015;25(7):1043‐1055. [doi:10.1101/gr.186072.114](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4484387/)
3. Pierre-Alain Chaumeil, Aaron J Mussig, Philip Hugenholtz, Donovan H Parks, GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database, Bioinformatics, Volume 36, Issue 6, 15 March 2020, Pages 1925–1927, [https://doi.org/10.1093/bioinformatics/btz848](https://doi.org/10.1093/bioinformatics/btz848)
