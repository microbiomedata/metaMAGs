# Validation Workflow

## Summary 

The validation workflow is meant to compare test data in json format available on the NMDC github against data generated by user test input. The purpose is to ensure user test data matches to or falls within an acceptable range of the NMDC test data in an effort to push towards reproducible data standards. The **test_mags.wdl** file contains most of the neccesary inputs needed for the MAGS worfklow and will grab user generated json file, but will require user input for the path of the database in their local environment and local user paths for test data to be updated in the input.json. 

For the metaMAGS workflow test validation, the user will need to get and unpack the test data:
```
mkdir -p test_data
wget https://portal.nersc.gov/cfs/m3408/test_data/metaMAGs_small_test_data.tgz
tar xvzf metaMAGs_small_test_data.tgz -C test_data
rm metaMAGs_small_test_data.tgz
```
sam, contig, gff, and map files will be derived from this directory for testing

## Input files for rqc_test.wdl
1. Binning contianer(provided),  
2. Comparejson container(provided), 
3. Outdir (specified by user),
4. Sam File (user specified),
5. Contig File (user specified),
6. GFF File (user specified),
7. Poject prefix (provied, can be changed by user),
8. url for NMDC test json (provided for small test data),
9. CPU 
10. pplacer cpu allocation
```
workflow test_mags {
    String container="microbiomedata/nmdc_mbin:0.1.2"
    String validate_container="microbiomedata/comparejson:0.1"
    String gtdbtk_database="/vol_b/nmdc_workflows/data/refdata/GTDBTK_GB"
    String outdir="/vol_b/nmdc_workflows/test_nmdc/metaMAGS/test_dir"
    File sam_file="/vol_b/nmdc_workflows/data/public/pairedMapped_sorted.bam"
    File contig_file="/vol_b/nmdc_workflows/data/public/Ga0482263_contigs.fna"
    File gff_file="/vol_b/nmdc_workflows/data/public/Ga0482263_functional_annotation.gff"
    File? map_file="/vol_b/nmdc_workflows/data/public/Ga0482263_contig_names.map.txt"
    String proj_name= "test"
    File ref_json="/vol_b/nmdc_workflows/test_nmdc/metaMAGs/small_test/small_output/MAGs_stats.json"
    Int cpu=32
    Int pplacer_cpu=1
```
## Docker contianers can be found here:
Bin: [microbiomedata/nmdc_mbin:0.1.2](https://hub.docker.com/r/microbiomedata/nmdc_mbin)
Comparjson: [microbiomedata/comparejson:0.1](https://hub.docker.com/r/microbiomedata/comparejson)

## Running Testing Validation Workflow

The command for running test validation is similar to that found in the submit.sh file, with the exception of switching out mbin_nmdc.wdl for test_mags.wdl.

 - `rqc_test.wdl` file: the WDL file for workflow definition
 - `input.json` file: the test input for the workflow
 - `cromwell.conf` file: the conf file for running Cromwell.
 -  `cromwell.jar` file: the jar file for running Cromwell.
 -  `metadata_out.json` file: file collects run data, will be created after run of command
 
Example:
```
java -Dconfig.file=cromwell.conf -jar cromwell.jar run -m metadata_out.json -i input.json test_mags.wdl
```

## Validation Metric
Validation metric is determined through a printed command line statement that will read:
```
"test.validate.result": ["No differences detected: test validated"]
```
or 
```
"test.validate.result": ["Test Failed"]
```

If test fails, please check inputs or contact local system administrators to ensure there are no system issues causing discrepency in results. 