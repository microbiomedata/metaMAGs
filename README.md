#  Metagenome assembled genomes generation workflow

## Summary

The workflow is based on [IMG MAGs pipeline](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6323987/)<sup>1</sup> for metagenome assembled genomes generation. It takes assembled contigs, bam file from reads mapping to contigs and [contigs annotations](https://github.com/microbiomedata/mg_annotation) result to associate groups of contigs as deriving from a seemingly coherent microbial species (binning) and evaluted by checkM, gtdb-tk and eukcc. 

## Required Database

* [CheckM2](https://pubmed.ncbi.nlm.nih.gov/37500759/)<sup>2</sup> database is ~3GB contains the databases used for the Metagenome Binned contig quality assessment. (requires 90GB+ of memory) 
    * https://zenodo.org/records/14897628

* [GTDB-Tk](https://doi.org/10.1093/bioinformatics/btz848)<sup>3</sup> requires ~78G of external data that need to be downloaded and unarchived. (requires ~150GB of memory)
    * https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz
    * https://portal.nersc.gov/cfs/m3408/databases/mash_sketch_db_r220.msh

* [EuKCC2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02155-4)<sup>4</sup> requires ~12G of external data that need to be downloaded and unarchived.
    * http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.2.tar.gz


* Prepare the Database

```bash
    wget https://zenodo.org/records/14897628/files/checkm2_database.tar.gz
    tar -xvzf checkm2_database.tar.gz
    mkdir -p refdata/CheckM_DB && tar -xvzf checkm2_database.tar.gz -C refdata/CheckM_DB
    rm checkm2_database.tar.gz
```

```bash
    wget https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz
    mkdir -p refdata/GTDBTK_DB && tar -xvzf gtdbtk_r220_data.tar.gz
    mv release220 refdata/GTDBTK_DB
    rm gtdbtk_r220_data.tar.gz
    wget https://portal.nersc.gov/cfs/m3408/databases/mash_sketch_db_r220.msh
    mv mash_sketch_db_r220.msh refdata/GTDBTK_DB
```

```bash
    wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.2.tar.gz
    tar -xvzf eukcc2_db_ver_1.2.tar.gz
    mkdir -p refdata/
    mv eukcc2_db_ver_1.2 refdata/EUKCC2_DB
    rm eukcc2_db_ver_1.2.tar.gz
```

## Running Workflow in Cromwell

Description of the files:
 - `.wdl` file: the WDL file for workflow definition
 - `.json` file: the example input for the workflow
 - `.conf` file: the conf file for running Cromwell.
 - `.sh` file: the shell script for running the example workflow (sbatch)

## The Docker images

[njvarghese/semibin:v2.2.0-pyt](https://hub.docker.com/r/njvarghese/semibin)

[checkm2:1.0.2](https://quay.io/repository/biocontainers/checkm2?tab=tags&tag=1.0.2--pyh7cba7a3_0)

[ecogenomic/gtdbtk:2.4.0](https://hub.docker.com/r/ecogenomic/gtdbtk)

[doejgi/eukcc:2.1.2](https://hub.docker.com/r/doejgi/eukcc)

[microbiomedata/nmdc_mbin_vis](https://hub.docker.com/r/microbiomedata/nmdc_mbin_vis)

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
11. Tab delimited file for [CRISPR](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-209) annotation.
12. Tab delimited file for Gene Product name assignment.
13. Tab delimited file for Gene Phylogeny assignment.
14. Tab delimited file for Contig/Scaffold lineage.
15. nmdc_mags.map_file: MAP file containing mapping of contig headers to annotation IDs 
16. seqtype: "short_reads" or "long_reads"
17. GTDBTK Database
18. GTDBTK Mash Database 
19. CheckM2 Database
20. EukCC Database
21. nmdc_mags.use_gpu: Use GPU or CPU
22. (optional) nmdc_mags.threads: The number of threads used by metabat/samtools/checkm/gtdbtk. default: 64


```
{
    "nmdc_mags.proj_name": "nmdc_wfmag-xx-xxxxxxxx",
    "nmdc_mags.contig_file": "/path/to/Assembly/nmdc_wfmgas-xx-xxxxxxx_contigs.fna",
    "nmdc_mags.sam_file": "/path/to/Assembly/nmdc_wfmgas-xx-xxxxxxx_pairedMapped_sorted.bam",
    "nmdc_mags.gff_file": "/path/to/Annotation/nmdc_wfmgan-xx-xxxxxxx_functional_annotation.gff",
    "nmdc_mags.proteins_file": "/path/to/Annotation/nmdc_wfmgan-xx-xxxxxxx_proteins.faa",
    "nmdc_mags.cog_file": "/path/to/Annotation/nmdc_wfmgan-xx-xxxxxxx_cog.gff",
    "nmdc_mags.ec_file": "/path/to/Annotation/nmdc_wfmgan-xx-xxxxxxx_ec.tsv",
    "nmdc_mags.ko_file": "/path/to/Annotation/nmdc_wfmgan-xx-xxxxxxx_ko.tsv",
    "nmdc_mags.pfam_file": "/path/to/Annotation/nmdc_wfmgan-xx-xxxxxxx_pfam.gff",
    "nmdc_mags.tigrfam_file": "/path/to/Annotation/nmdc_wfmgan-xx-xxxxxxxtigrfam.gff",
    "nmdc_mags.crispr_file": "/path/to/Annotation/nmdc_wfmgan-xx-xxxxxxx_crt.crisprs,
    "nmdc_mags.product_names_file": "/path/to/Annotation/nmdc_wfmgan-xx-xxxxxxx_product_names.tsv",
    "nmdc_mags.gene_phylogeny_file": "/path/to/Annotation/nmdc_wfmgan-xx-xxxxxxx_gene_phylogeny.tsv",
    "nmdc_mags.lineage_file": "/path/to/Annotation/nmdc_wfmgan-xx-xxxxxxx_scaffold_lineage.tsv",
    "nmdc_mags.map_file": "/path/to/Annotation/nmdc_wfmgan-xx-xxxxxxx_contig_names_mapping.tsv",
    "nmdc_mags.seqtype": "short_reads",
    "nmdc_mags.gtdbtk_db": "refdata/GTDBTK_DB/release220",
    "nmdc_mags.gtdbtk_mash_db": "refdata/GTDBTK_DB/mash_sketch_db_r220.msh",
    "nmdc_mags.checkm2_db": "refdata/CheckM2_database/uniref100.KO.1.dmnd",
    "nmdc_mags.eukcc2_db": "refdata/EUKCC2_DB/eukcc2_db_ver_1.2",
    "nmdc_mags.use_gpu": true
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
|-- project_name_checkm_qa.out
|-- project_name_gtdbtk.json
|-- project_name_heatmap.pdf  (The Heatmap presents the pdf file containing the KO analysis results for metagenome bins)
|-- project_name_barplot.pdf  (The Bar chart presents the pdf file containing the KO analysis results for metagenome bins
|-- project_name_kronaplot.html  (The Krona plot presents the HTML file containing the KO analysis results for metagenome bins)
```

### Citation
1. Chen IA, Chu K, Palaniappan K, et al. IMG/M v.5.0: an integrated data management and comparative analysis system for microbial genomes and microbiomes. Nucleic Acids Res. 2019;47(D1):D666‐D677. [doi:10.1093/nar/gky901](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6323987/)
2. Chklovski A, Parks DH, Woodcroft BJ, Tyson GW. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. Nat Methods. 2023 Aug;20(8):1203-1212. [doi: 10.1038/s41592-023-01940-w.](https://pubmed.ncbi.nlm.nih.gov/37500759/)
3. Pierre-Alain Chaumeil, Aaron J Mussig, Philip Hugenholtz, Donovan H Parks, GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database, Bioinformatics, Volume 36, Issue 6, 15 March 2020, Pages 1925–1927, [doi.org/10.1093/bioinformatics/btz848](https://doi.org/10.1093/bioinformatics/btz848)
4. Saary, Paul, Alex L. Mitchell, and Robert D. Finn. "Estimating the quality of eukaryotic genomes recovered from metagenomic analysis with EukCC." Genome biology 21.1 (2020): 1-21. [doi.org/10.1186/s13059-020-02155-4](https://doi.org/10.1186/s13059-020-02155-4)
5. Pan S, Zhu C, Zhao XM, Coelho LP. A deep siamese neural network improves metagenome-assembled genomes in microbiome datasets across different environments. Nat Commun. 2022 Apr 28;13(1):2326. [doi: 10.1038/s41467-022-29843-y.](https://pubmed.ncbi.nlm.nih.gov/35484115/)
