workflow checkM {
    File contigs_file_zip
    String container = "microbiomedata/nmdc_mbin:0.1.4"
    Int cpu = 16
    parameter_meta {
        contigs_file_zip: "contigs fasta files in a zip file"
        checkm: "metabat bin checkm result"
    }
    call run_checkM {
        input: contigs_file_zip = contigs_file_zip, container=container, cpu = cpu
    }
    output{
        File? checkm = run_checkM.checkm
    }
}

task run_checkM {
    File contigs_file_zip
    String container
    Int cpu 
    runtime {
        docker: container
        memory: "60 GiB"
        cpu:  cpu
    }

    command {
       mkdir -p contigs_dir
       unzip ${contigs_file_zip} -d contigs_dir
       # set TMPDIR to avoid AF_UNIX path too long error
       export TMPDIR=/tmp
       checkm lineage_wf --pplacer_threads 4 -x fa -t ${cpu} contigs_dir checkm_out
       checkm qa -f checkm_qa.out checkm_out/lineage.ms checkm_out
    }
    output {
        File? checkm = "checkm_qa.out"
    }
}
