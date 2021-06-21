#  Metagenome assembled genomes generation workflow

## Summary

The workflow is based on [IMG MAGs pipeline](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6323987/)<sup>1</sup> for metagenome assembled genomes generation. It takes assembled contigs, reads mapping result bam file and contigs annotations result to to associate groups of contigs as deriving from a seemingly coherent microbial species (binning) and evaluted by checkM and gtdb-tk. 

## Required Database

* [CheckM](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4484387/)<sup>2</sup> database is 275MB contains the databases used for the Metagenome Binned contig quality assessment. (requires 40GB+ of memory, included in the image) 
    * https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz

* [GTDB-Tk](https://doi.org/10.1093/bioinformatics/btz848)<sup>3</sup> requires ~27G of external data that need to be downloaded and unarchived. (requires ~100GB of memory)
    * https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/gtdbtk_r89_data.tar.gz

* Prepare the GTDB-Tk Database

```bash
    
    wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/gtdbtk_r89_data.tar.gz
    tar -xvzf gtdbtk_r89_data.tar.gz
    mv release89 refdata/GTDBTK_DB
    
    rm gtdbtk_r89_data.tar.gz
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

1. Number of CPUs,
2. The number of threads used by pplacer (Use lower number to reduce the memory usage)
3. Output directory
4. Project name
5. Metagenome Assembled Contig fasta file
6. Sam/Bam file from reads mapping back to contigs.
7. Contigs functional annotation result in gff format
8. Optioal: Tab-delimited text file which containing mapping of headers between SAM and FNA (ID in SAM/FNA<tab>ID in GFF). A two column tab-delimited file. When the annotation and assembly are performed using different identifiers for contigs. The map file is to link the gff file content and mapping result bam file content to the assembled contigs ID.
9. The database directory path which includes `checkM_DB` and `GTDBTK_DB` subdirectories. 
10. (optional) scratch_dir: use --scratch_dir for gtdbtk disk swap to reduce memory usage but longer runtime

```
{
  "nmdc_mags.cpu":32,
  "nmdc_mags.pplacer_cpu":1,
  "nmdc_mags.outdir":"/global/cfs/cdirs/m3408/aim2/metagenome/MAGs/output",
  "nmdc_mags.proj_name":"3300037552",
  "nmdc_mags.contig_file":"/global/cfs/cdirs/m3408/aim2/metagenome/MAGs/mbin-nmdc-test-dataset/3300037552.a.fna",
  "nmdc_mags.sam_file":"/global/cfs/cdirs/m3408/aim2/metagenome/MAGs/mbin-nmdc-test-dataset/3300037552.bam.sorted.bam",
  "nmdc_mags.gff_file":"/global/cfs/cdirs/m3408/aim2/metagenome/MAGs/mbin-nmdc-test-dataset/3300037552.a.gff",
  "nmdc_mags.map_file":"/global/cfs/cdirs/m3408/aim2/metagenome/MAGs/mbin-nmdc-test-dataset/3300037552.a.map.txt",
  "nmdc_mags.gtdbtk_database":"/path/to/GTDBTK_DB"
}
```

## Output files

The output will have a bunch of output directories, files, including statistical numbers, status log and a shell script to reproduce the steps etc. 

The final [MiMAG](https://www.nature.com/articles/nbt.3893#Tab1) output is in `hqmq-metabat-bins` directory and its corresponding lineage result in `gtdbtk_output` directory.

```
|-- MAGs_stats.json
|-- 3300037552.bam.sorted
|-- 3300037552.depth
|-- 3300037552.depth.mapped
|-- bins.lowDepth.fa
|-- bins.tooShort.fa
|-- bins.unbinned.fa
|-- checkm-out
|   |-- bins/
|   |-- checkm.log
|   |-- lineage.ms
|   `-- storage
|-- checkm_qa.out
|-- gtdbtk_output
|   |-- align/
|   |-- classify/
|   |-- identify/
|   |-- gtdbtk.ar122.classify.tree -> classify/gtdbtk.ar122.classify.tree
|   |-- gtdbtk.ar122.markers_summary.tsv -> identify/gtdbtk.ar122.markers_summary.tsv
|   |-- gtdbtk.ar122.summary.tsv -> classify/gtdbtk.ar122.summary.tsv
|   |-- gtdbtk.bac120.classify.tree -> classify/gtdbtk.bac120.classify.tree
|   |-- gtdbtk.bac120.markers_summary.tsv -> identify/gtdbtk.bac120.markers_summary.tsv
|   |-- gtdbtk.bac120.summary.tsv -> classify/gtdbtk.bac120.summary.tsv
|   `-- ..etc 
|-- hqmq-metabat-bins
|   |-- bins.11.fa
|   |-- bins.13.fa
|   `-- ... etc 
|-- mbin-2020-05-24.sqlite
|-- mbin-nmdc.20200524.log
|-- metabat-bins
|   |-- bins.1.fa
|   |-- bins.10.fa
|   `-- ... etc 
```

### Citation
1. Chen IA, Chu K, Palaniappan K, et al. IMG/M v.5.0: an integrated data management and comparative analysis system for microbial genomes and microbiomes. Nucleic Acids Res. 2019;47(D1):D666‐D677. [doi:10.1093/nar/gky901](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6323987/)
2. Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome Res. 2015;25(7):1043‐1055. [doi:10.1101/gr.186072.114](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4484387/)
3. Pierre-Alain Chaumeil, Aaron J Mussig, Philip Hugenholtz, Donovan H Parks, GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database, Bioinformatics, Volume 36, Issue 6, 15 March 2020, Pages 1925–1927, [https://doi.org/10.1093/bioinformatics/btz848](https://doi.org/10.1093/bioinformatics/btz848)
