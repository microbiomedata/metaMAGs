workflow nmdc_mags {
    String? outdir
    String  proj_name
    File contig_file
    File sam_file
    File gff_file
    String container = "microbiomedata/nmdc_mbin:0.1.4"
    File? map_file
    File? domain_file
    String? scratch_dir
    Int cpu=32
    Int pplacer_cpu=1
    String? proj = "MAGs" 
    String? informed_by = "${proj}"  # "nmdc:xxxxxxxxxxxxxxx"
    String resource = "NERSC - Cori"
    String url_root = "https://data.microbiomedata.org/data/"
    String git_url = "https://github.com/microbiomedata/metaMAGs/releases/tag/1.0.4"
    
    String gtdbtk_database="/refdata/GTDBTK_DB"

    call mbin_nmdc{
         input:  name=proj_name, 
                 fasta=contig_file, 
                 sam=sam_file, 
                 gff=gff_file, 
                 map=map_file, 
                 domain=domain_file,
                 cpu=cpu,
                 pplacer_cpu=pplacer_cpu,
                 scratch_dir=scratch_dir,
                 database=gtdbtk_database,
                 container=container
    }
    call generate_objects {
         input: container="microbiomedata/workflowmeta:1.0.5.1",
                proj=proj,
                start = mbin_nmdc.start,
                informed_by = "${informed_by}",
                resource = "${resource}",
                url_base = "${url_root}",
                git_url = "${git_url}",
                fasta = "${contig_file}",
                bam = "${sam_file}",
                functional_gff ="${gff_file}",
                lowdepth = mbin_nmdc.low,
                unbinned = mbin_nmdc.unbinned,
                short = mbin_nmdc.short,
                checkm = mbin_nmdc.checkm,
                json_stats = mbin_nmdc.json_stats,
                bac_summary = mbin_nmdc.bacsum,
                ar_summary = mbin_nmdc.arcsum,
                metabat_bin_fasta_files = mbin_nmdc.bin_fasta_files,
                hqmq_bin_fasta_files = mbin_nmdc.hqmq_bin_fasta_files
    }
    if (defined(outdir)){
        call make_output {
            input: container="microbiomedata/workflowmeta:1.0.5.1",
                   activity_json=generate_objects.activity_json,
                   object_json=generate_objects.data_object_json,
                   short=mbin_nmdc.short,
                   low=mbin_nmdc.low,
                   unbinned=mbin_nmdc.unbinned,
                   json_stats=mbin_nmdc.json_stats,
                   hqmq_bin_fasta_zip=generate_objects.hqmq_bin_fasta_zip,
                   bin_fasta_zip=generate_objects.metabat_bin_fasta_zip,
                   checkm=mbin_nmdc.checkm,
                   gtdbtk_bac_summary=mbin_nmdc.bacsum,
                   gtdbtk_ar_summary=mbin_nmdc.arcsum,
                   outdir=outdir
                   
        }
    }
    output {
        Array[File] hqmq_bin_fasta_files = mbin_nmdc.hqmq_bin_fasta_files
        Array[File] bin_fasta_files = mbin_nmdc.bin_fasta_files
        File? final_hqmq_bins = make_output.hqmq_bin_fa_zip
        File? metabat_bins = make_output.metabat_bin_fa_zip
        File  activityjson=generate_objects.activity_json
        File  objectjson=generate_objects.data_object_json
        File? final_checkm = make_output.checkm_output
        File? final_gtdbtk_bac_summary = make_output.bac_summary
        File? final_gtdbtk_ar_summary = make_output.ar_summary
        File? final_tooShort_fa = make_output.tooShort_fa
        File? final_lowDepth_fa = make_output.lowDepth_fa
        File? final_unbinned_fa = make_output.unbinned_fa
        File? final_stats = make_output.stats
        File short = mbin_nmdc.short
        File low = mbin_nmdc.low
        File unbinned = mbin_nmdc.unbinned
        File? checkm = mbin_nmdc.checkm

    }
    parameter_meta {
        cpu: "number of CPUs"
        pplacer_cpu: "number of threads used by pplacer"
        outdir: "the final output directory path"
        scratch_dir: "use --scratch_dir for gtdbtk disk swap to reduce memory usage but longer runtime"
        proj_name: "project name"
        contig_file: "input assembled contig fasta file"
        sam_file: "Sam/Bam file from reads mapping back to contigs. [sam.gz or bam]"
        gff_file: "contigs functional annotation result in gff format"
        map_file: "text file which containing mapping of headers between SAM and FNA (ID in SAM/FNA ID in GFF)"
        database: "GTDBTK_DB database directory path"
        final_hqmq_bins: "high quality and medium quality bin fasta output"
        metabat_bins: "initial metabat bining result fasta output"
        final_checkm: "metabat bin checkm result"
        final_gtdbtk_bac_summary: "gtdbtk bacterial assignment result summary table"
        final_gtdbtk_ar_summary: "gtdbtk archaea assignment result summary table"
        final_tooShort_fa: "tooShort (< 3kb) filtered contigs fasta file by metaBat2"
        final_lowDepth_fa: "lowDepth (mean cov <1 )  filtered contigs fasta file by metabat2"
        final_unbinned_fa: "unbinned fasta file from metabat2"
        final_stats: "statistics summary in json format"
        activityjson: "nmdc activity json file"
        objectjson: "nmdc data object json file"
    }
    meta {
        author: "Chienchi Lo, B10, LANL"
        email: "chienchi@lanl.gov"
        version: "1.0.4"
    }

}

task mbin_nmdc {
    String name
    File fasta
    File sam
    File gff
    File? map
    File? domain
    String database
    Int cpu
    Int pplacer_cpu
    String? scratch_dir
    String container
    String filename_outlog="stdout.log"
    String filename_errlog="stderr.log"
    String filename_stat="checkm_qa.out"
    String dollar="$"
    runtime {
        docker: container
        memory: "120 GiB"
        cpu:  cpu
        database: database
    }   

     command {
        export TIME="time result\ncmd:%C\nreal %es\nuser %Us \nsys  %Ss \nmemory:%MKB \ncpu %P"
        set -eo pipefail
        # Capture the start time
        date --iso-8601=seconds > start.txt
        # set TMPDIR to avoid AF_UNIX path too long error 
        export TMPDIR=/tmp
        #export TMPDIR=$HOME 
        export GTDBTK_DATA_PATH=/databases
        mbin_nmdc.py ${"--map " + map} ${"--domain " + domain} ${"--scratch_dir " + scratch_dir} --pplacer_cpu ${pplacer_cpu} --cpu ${cpu} ${name} ${fasta} ${sam} ${gff}
        mbin_stats.py $PWD
     }
     output {
        File runScript = "script"
        File? stat = filename_stat
        String start = read_string("start.txt")
        File short = "bins.tooShort.fa"
        File low = "bins.lowDepth.fa"
        File unbinned = "bins.unbinned.fa"
        File? checkm = "checkm_qa.out"
        File json_stats= "MAGs_stats.json"
        File? bacsum = "gtdbtk_output/gtdbtk.bac120.summary.tsv"
        File? arcsum = "gtdbtk_output/gtdbtk.ar122.summary.tsv"
        Array[File] hqmq_bin_fasta_files = glob("hqmq-metabat-bins/*fa")
        Array[File] bin_fasta_files = glob("metabat-bins/*fa")
     }
}

task generate_objects{
    String container
    String proj
    String start
    String informed_by
    String resource
    String url_base
    String git_url
    File fasta
    File bam
    File? checkm
    File functional_gff
    File short
    File lowdepth
    File unbinned
    File json_stats
    
    File? bac_summary
    File? ar_summary
    Array[File] metabat_bin_fasta_files
    Array[File] hqmq_bin_fasta_files
    String dollar="$"
   
    command<<<
        set -e
        end=`date --iso-8601=seconds`
        ### set IFS to avoid space separate string and save into the outputs array elements
        IFS=""
        outputs=()
        zip -j hqmq-metabat-bins.zip ${sep=" " hqmq_bin_fasta_files} || true
        zip -j metabat-bins.zip ${sep=" " metabat_bin_fasta_files} || true
        [ -e hqmq-metabat-bins.zip ] && outputs+=( "hqmq-metabat-bins.zip" ) && outputs+=( "high quality and medium quality bin fasta output" )
        [ -e maetabat-bins.zip ] && outputs+=( "metabat-bins.zip" ) && outputs+=( "initial metabat bining result fasta output" )

        /scripts/generate_objects.py --type "MAGs" --id ${informed_by} \
             --name "MAGs Analysis Activity for ${proj}" --part ${proj} \
             --start ${start} --end $end \
             --extra ${json_stats} \
             --resource '${resource}' --url ${url_base}${proj}/MAGs/ --giturl ${git_url} \
             --inputs ${fasta} ${bam} ${functional_gff} \
             --outputs ${dollar}{outputs[@]} \
             ${short} "tooShort (< 3kb) filtered contigs fasta file by metaBat2" \
             ${lowdepth} "lowDepth (mean cov <1 )  filtered contigs fasta file by metabat2" \
             ${unbinned} "unbinned fasta file from metabat2" \
             ${" " + checkm + " \"metabat2 bin checkm quality assessment result\""} \
             ${" " + bac_summary + " \"gtdbtk bacterial assignment result summary table\""} \
             ${" " + ar_summary + " \"gtdbtk archaea assignment result summary table\""}  
    >>>
    runtime {
        docker: container
        memory: "10 GiB"
        cpu:  1
    }
    output{
        File activity_json = "activity.json"
        File data_object_json = "data_objects.json"
        File? metabat_bin_fasta_zip = "hqmq-metabat-bins.zip"
        File? hqmq_bin_fasta_zip = "metabat-bins.zip"
    }
}


task make_output{
    String outdir
    File short
    File low
    File unbinned
    File? hqmq_bin_fasta_zip
    File? bin_fasta_zip
    File? checkm
    File json_stats
    File? gtdbtk_bac_summary
    File? gtdbtk_ar_summary
    String container
    File activity_json
    File object_json
 
    command{
        mkdir -p ${outdir}
        cp ${short} ${low} ${unbinned} ${json_stats} ${checkm} \
                   ${gtdbtk_bac_summary} ${gtdbtk_ar_summary} \
                   ${activity_json} ${object_json} \
                   ${outdir}
        # These may not exist
        ${"cp " + hqmq_bin_fasta_zip + " " + outdir} 
        ${"cp " + bin_fasta_zip + " " + outdir}
        chmod 755 -R ${outdir}
    }
    output {
        File? hqmq_bin_fa_zip = "${outdir}/hqmq-metabat-bins.zip"
        File? metabat_bin_fa_zip = "${outdir}/metabat-bins.zip"
        File? checkm_output = "${outdir}/checkm_qa.out"
        File? bac_summary = "${outdir}/gtdbtk.bac120.summary.tsv"
        File? ar_summary = "${outdir}/gtdbtk.ar122.summary.tsv"
        File unbinned_fa = "${outdir}/bins.unbinned.fa"
        File tooShort_fa = "${outdir}/bins.tooShort.fa"
        File lowDepth_fa = "${outdir}/bins.lowDepth.fa"
        File stats = "${outdir}/MAGs_stats.json"
        File outactivity = "${outdir}/activity.json"
        File outobject = "${outdir}/data_objects.json"
    }
    runtime {
        docker: container
        memory: "1 GiB"
        cpu:  1
    }
}

