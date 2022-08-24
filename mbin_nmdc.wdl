workflow nmdc_mags {
    String? outdir
    String  proj_name
    File contig_file
    File sam_file
    File gff_file
    File proteins_file
    File cog_file
    File ec_file
    File ko_file
    File pfam_file
    File tigrfam_file
    File cath_funfam_file
    File smart_file
    File supfam_file
    File product_names_file
    File gene_phylogeny_file

    String container = "microbiomedata/nmdc_mbin:0.1.6"
    File? map_file
    File? domain_file
    String? scratch_dir
    Int cpu=32
    Int pplacer_cpu=1

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

    call package {
         input:  bins=mbin_nmdc.hqmq_bin_fasta_files,
                 json_stats=mbin_nmdc.stats_json,
                 gff_file=gff_file,
                 proteins_file=proteins_file,
                 cog_file=cog_file,
                 ec_file=ec_file,
                 ko_file=ko_file,
                 pfam_file=pfam_file,
                 tigrfam_file=tigrfam_file,
                 cath_funfam_file=cath_funfam_file,
                 smart_file=smart_file,
                 supfam_file=supfam_file,
                 product_names_file=product_names_file
    }

    if (defined(outdir)){
        call make_output {
            input: short=mbin_nmdc.short,
                   low=mbin_nmdc.low,
                   unbinned=mbin_nmdc.unbinned,
                   hqmq_bin_fasta_files=mbin_nmdc.hqmq_bin_fasta_files,
                   bin_fasta_files=mbin_nmdc.bin_fasta_files,
                   checkm=mbin_nmdc.checkm,
                   gtdbtk_bac_summary=mbin_nmdc.bacsum,
                   gtdbtk_ar_summary=mbin_nmdc.arcsum,
                   json_stats=mbin_nmdc.stats_json,
                   tsv_stats=mbin_nmdc.stats_tsv,
                   outdir=outdir
        }
    }

    output {
        File? final_hqmq_bins_zip = make_output.hqmq_bin_fasta_zip
        File? metabat_bins_zip = make_output.metabat_bin_fasta_zip
        File? final_checkm = make_output.checkm_output
        File? final_gtdbtk_bac_summary = make_output.bac_summary
        File? final_gtdbtk_ar_summary = make_output.ar_summary
        File? final_tooShort_fa = make_output.tooShort_fa
        File? final_lowDepth_fa = make_output.lowDepth_fa
        File? final_unbinned_fa = make_output.unbinned_fa
        File short = mbin_nmdc.short
        File low = mbin_nmdc.low
        File unbinned = mbin_nmdc.unbinned
        File? checkm = mbin_nmdc.checkm
        File stats_json = mbin_nmdc.stats_json
        File stats_tsv = mbin_nmdc.stats_tsv
        Array[File] hqmq_bin_fasta_files = mbin_nmdc.hqmq_bin_fasta_files
        Array[File] bin_fasta_files = mbin_nmdc.bin_fasta_files
        Array[File] hqmq_bin_tarfiles = package.hqmq_bin_tarfiles
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
        stats_json: "statistics summary in json format"
        stats_tsv: "statistics summary in tsv format"
    }
    meta {
        author: "Chienchi Lo, B10, LANL"
        email: "chienchi@lanl.gov"
        version: "1.0.1"
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
     runtime {
         docker: container
         memory: "120 GiB"
         cpu:  cpu
         database: database
     }

     command {
         export TIME="time result\ncmd:%C\nreal %es\nuser %Us \nsys  %Ss \nmemory:%MKB \ncpu %P"
         set -eo pipefail
         # set TMPDIR to avoid AF_UNIX path too long error
         export TMPDIR=/tmp
         export GTDBTK_DATA_PATH=${database}
         mbin_nmdc.py ${"--map " + map} ${"--domain " + domain} ${"--scratch_dir " + scratch_dir} --pplacer_cpu ${pplacer_cpu} --cpu ${cpu} ${name} ${fasta} ${sam} ${gff}
         mbin_stats.py $PWD
         touch MAGs_stats.tsv
     }
     output {
         File runScript = "script"
         File? stat = filename_stat
         File short = "bins.tooShort.fa"
         File low = "bins.lowDepth.fa"
         File unbinned = "bins.unbinned.fa"
         File? checkm = "checkm_qa.out"
         File? bacsum = "gtdbtk.bac120.summary.tsv"
         File? arcsum = "gtdbtk.ar122.summary.tsv"
         File stats_json = "MAGs_stats.json"
         File stats_tsv = "MAGs_stats.tsv"
         Array[File] hqmq_bin_fasta_files = glob("hqmq-metabat-bins/*fa")
         Array[File] bin_fasta_files = glob("metabat-bins/*fa")
     }
}

task package{
     Array[File] bins
     File json_stats
     File gff_file
     File proteins_file
     File cog_file
     File ec_file
     File ko_file
     File pfam_file
     File tigrfam_file
     File cath_funfam_file
     File smart_file
     File supfam_file
     File product_names_file
     String container="microbiomedata/nmdc_mbin_test:0.0.0"

     command {
             python3 /opt/conda/bin/create_tarfiles.py \
                     ${json_stats} ${gff_file} ${proteins_file} ${cog_file} \
                     ${ec_file} ${ko_file} ${pfam_file} ${tigrfam_file} \
                     ${cath_funfam_file} ${smart_file} ${supfam_file} \
                     ${product_names_file} \
                     ${sep=" " bins}
     }
     output {
         Array[File] hqmq_bin_tarfiles = glob("*tar.gz")
     }
     runtime {
         docker: container
         memory: "1 GiB"
         cpu:  1
     }
}

task make_output{
     File short
     File low
     File unbinned
     Array[File] hqmq_bin_fasta_files
     Array[File] bin_fasta_files
     File? checkm
     File? gtdbtk_bac_summary
     File? gtdbtk_ar_summary
     File json_stats
     File tsv_stats
     String? outdir
     String container="scanon/nmdc-meta:v0.0.2"

     command{
         zip hqmq-metabat-bins.zip ${sep=" " hqmq_bin_fasta_files}
         zip metabat-bins.zip ${sep=" " bin_fasta_files}
         if [ ! -z ${outdir} ] ; then
             mkdir -p ${outdir}
             cp ${short} ${low} ${unbinned} ${json_stats} ${tsv_stats}\
                ${gtdbtk_bac_summary} ${gtdbtk_ar_summary} \
                ${outdir}
             # These may not exist
             [ -e hqmq-metabat-bins.zip ] && cp hqmq-metabat-bins.zip ${outdir}
             [ -e metabat-bins.zip ] && cp metabat-bins.zip ${outdir}
             chmod 764 -R ${outdir}
         fi
    }
    output {
        String checkm_output = "${outdir}/checkm_qa.out"
        String bac_summary = "${outdir}/gtdbtk.bac120.summary.tsv"
        String ar_summary = "${outdir}/gtdbtk.ar122.summary.tsv"
        String unbinned_fa = "${outdir}/bins.unbinned.fa"
        String tooShort_fa = "${outdir}/bins.tooShort.fa"
        String lowDepth_fa = "${outdir}/bins.lowDepth.fa"
        String hqmq_bin_fasta_zip = "${outdir}/hqmq-metabat-bins.zip"
        String metabat_bin_fasta_zip = "${outdir}/metabat-bins.zip"
        File? hqmq_bin_zip = "hqmq-metabat-bins.zip"
        File? metabat_bin_zip = "metabat-bins.zip"
        File stats = "${outdir}/MAGs_stats.json"
        File tsvstats = "${outdir}/MAGs_stats.tsv"
    }
    runtime {
        docker: container
        memory: "1 GiB"
        cpu:  1
    }
}

