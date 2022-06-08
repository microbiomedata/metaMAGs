import "mbin_nmdc.wdl" as bins

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

  call bins.mbin_nmdc as mags {
    input:  name=proj_name,
            fasta=contig_file,
            sam=sam_file,
            gff=gff_file,
            map=map_file,
            cpu=cpu,
            pplacer_cpu=pplacer_cpu,
            database=gtdbtk_database,
            container=container
  }
  call bins.make_output as out {
   input: outdir= outdir, mbin_nmdc_output=mags.runScript
    }
   
  call validate {
    input: container=validate_container,
           refjson=ref_json,
           user_json=out.json_stats
  }
}
task validate {
   String container
   File refjson
   File user_json
   command {
       compare_json.py -i ${refjson} -f ${user_json}
   }
   output {
       Array[String] result = read_lines(stdout())
   }
   runtime {
     memory: "1 GiB"
     cpu:  1
     maxRetries: 1
     docker: container
   }
}
