workflow nmdc_mags {
    String informed_by
    String resource
    String url_root
    String git_url
    String  proj_name
    String contig_file
    String sam_file
    String gff_file
    String proteins_file
    String cog_file
    String ec_file
    String ko_file
    String pfam_file
    String tigrfam_file
    String cath_funfam_file
    String smart_file
    String supfam_file
    String product_names_file
    String gene_phylogeny_file
    String container = "microbiomedata/nmdc_mbin:0.1.6"
    File? map_file
    File? domain_file
    String? scratch_dir
    Int cpu=32
    Int pplacer_cpu=1
    String gtdbtk_database="/refdata/GTDBTK_DB"

    call stage {
        input:
            container=container,
            contig_file=contig_file,
            sam_file=sam_file,
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
            product_names_file=product_names_file,
            gene_phylogeny_file=gene_phylogeny_file
    }

    call mbin_nmdc{
         input:  name=proj_name,
                 fasta=stage.contig,
                 sam=stage.sam,
                 gff=stage.gff,
                 #map=stage.map_file,
                 #domain=stage.domain_file,
                 cpu=cpu,
                 pplacer_cpu=pplacer_cpu,
                 scratch_dir=scratch_dir,
                 database=gtdbtk_database,
                 container=container
    }

    call package {
         input:  bins=mbin_nmdc.hqmq_bin_fasta_files,
                 json_stats=mbin_nmdc.stats_json,
                 gff_file=stage.gff,
                 proteins_file=stage.proteins,
                 cog_file=stage.cog,
                 ec_file=stage.ec,
                 ko_file=stage.ko,
                 pfam_file=stage.pfam,
                 tigrfam_file=stage.tigrfam,
                 cath_funfam_file=stage.cath_funfam,
                 smart_file=stage.smart,
                 supfam_file=stage.supfam,
                 product_names_file=stage.product_names
    }

    call finish_mags {
        input:
        container="microbiomedata/workflowmeta:1.1.1",
        contigs=stage.contig,
        anno_gff=stage.gff,
        sorted_bam=stage.sam,
        proj=proj_name,
        start=stage.start,
        resource=resource,
        url_root=url_root,
        git_url=git_url,
        informed_by=informed_by,
        checkm = mbin_nmdc.checkm,
        bacsum= mbin_nmdc.bacsum,
        arcsum = mbin_nmdc.arcsum,
        short = mbin_nmdc.short,
        low = mbin_nmdc.low,
        unbinned = mbin_nmdc.unbinned,
        checkm = mbin_nmdc.checkm,
        stats_json = mbin_nmdc.stats_json,
        stats_tsv = mbin_nmdc.stats_tsv,
        hqmq_bin_fasta_files = mbin_nmdc.hqmq_bin_fasta_files,
        bin_fasta_files = mbin_nmdc.bin_fasta_files,
        hqmq_bin_tarfiles = package.hqmq_bin_tarfiles,

    }

    output {
        File? final_hqmq_bins_zip = finish_mags.final_hqmq_bins_zip
        File? final_gtdbtk_bac_summary = finish_mags.final_gtdbtk_bac_summary
        File? final_gtdbtk_ar_summary = finish_mags.final_gtdbtk_ar_summary
        File short = finish_mags.final_short
        File low = finish_mags.final_lowDepth_fa
        File final_unbinned_fa  = finish_mags.final_unbinned_fa
        File final_checkm = finish_mags.final_checkm
        File final_stats_json = finish_mags.final_stats_json
        File stats_tsv = mbin_nmdc.stats_tsv
        File mags_objects = finish_mags.objects
        # Array[File] hqmq_bin_fasta_files = mbin_nmdc.hqmq_bin_fasta_files
        # Array[File] bin_fasta_files = mbin_nmdc.bin_fasta_files
        # Array[File] hqmq_bin_tarfiles = package.hqmq_bin_tarfiles
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

task stage {
    String container
    String contig_file
    String sam_file
    String gff_file
    String proteins_file
    String cog_file
    String ec_file
    String ko_file
    String pfam_file
    String tigrfam_file
    String cath_funfam_file
    String smart_file
    String supfam_file
    String product_names_file
    String gene_phylogeny_file
    String contigs_out="contigs.fasta"
    String bam_out="pairedMapped_sorted.bam"
    String gff_out="functional_annotation.gff"
    String proteins_out="proteins.faa"
    String cog_out="cog.gff"
    String ec_out="ec.tsv"
    String ko_out="ko.tsv"
    String pfam_out="pfam.gff"
    String tigrfam_out="tigrfam.gff"
    String cath_funfam_out="cath_funfam.gff"
    String smart_out="smart.gff"
    String supfam_out="supfam.gff"
    String products_out="products.tsv"
    String gene_phylogeny_out="gene_phylogeny.tsv"

   command<<<

       set -e

        function stage() {
            in=$1
            out=$2
            if [ $( echo $in |egrep -c "https*:") -gt 0 ] ; then
                wget $in -O $out
            else
                ln $in $out || cp $in $out
            fi
        }

        stage ${contig_file} ${contigs_out}
        stage ${sam_file} ${bam_out}
        stage ${gff_file} ${gff_out}
        stage ${proteins_file} ${proteins_out}
        stage ${cog_file} ${cog_out}
        stage ${ec_file} ${ec_out}
        stage ${ko_file} ${ko_out}
        stage ${pfam_file} ${pfam_out}
        stage ${tigrfam_file} ${tigrfam_out}
        stage ${cath_funfam_file} ${cath_funfam_out}
        stage ${smart_file} ${smart_out}
        stage ${supfam_file} ${supfam_out}
        stage ${product_names_file} ${products_out}
        stage ${gene_phylogeny_file} ${gene_phylogeny_out}

       date --iso-8601=seconds > start.txt

    >>>

   output{
        File contig = "contigs.fasta"
        File sam = "pairedMapped_sorted.bam"
        File gff = "functional_annotation.gff"
        File proteins = "proteins.faa"
        File cog = "cog.gff"
        File ec = "ec.tsv"
        File ko = "ko.tsv"
        File pfam = "pfam.gff"
        File tigrfam = "tigrfam.gff"
        File cath_funfam = "cath_funfam.gff"
        File smart = "smart.gff"
        File supfam = "supfam.gff"
        File product_names = "products.tsv"
        File gene_phylogeny = "gene_phylogeny.tsv"
        String start = read_string("start.txt")
   }
   runtime {
     memory: "1 GiB"
     cpu:  2
     maxRetries: 1
     docker: container
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

     command <<<
         export TIME="time result\ncmd:%C\nreal %es\nuser %Us \nsys  %Ss \nmemory:%MKB \ncpu %P"
         set -eo pipefail
         # set TMPDIR to avoid AF_UNIX path too long error
         export TMPDIR=/tmp
         export GTDBTK_DATA_PATH=${database}
         mbin_nmdc.py ${"--map " + map} ${"--domain " + domain} ${"--scratch_dir " + scratch_dir} --pplacer_cpu ${pplacer_cpu} --cpu ${cpu} ${name} ${fasta} ${sam} ${gff}
         mbin_stats.py $PWD
         touch MAGs_stats.tsv


        if [ -f  gtdbtk.bac120.summary.tsv ]; then
            echo "bacterial summary exists."
        else
            echo "No Results" > gtdbtk.bac120.summary.tsv
        fi

        if [ -f  gtdbtk.ar122.summary.tsv ]; then
            echo "archeal summary exists."
        else
            echo "No Results" > gtdbtk.ar122.summary.tsv
        fi

     >>>
     output {
         File runScript = "script"
         File? stat = filename_stat
         File short = "bins.tooShort.fa"
         File low = "bins.lowDepth.fa"
         File unbinned = "bins.unbinned.fa"
         File? checkm = "checkm_qa.out"
         File bacsum = "gtdbtk.bac120.summary.tsv"
         File arcsum = "gtdbtk.ar122.summary.tsv"
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


task finish_mags {
    String container
    File contigs
    File anno_gff
    File sorted_bam
    String proj
    String prefix=sub(proj, ":", "_")
    String start
    String informed_by
    String resource
    String url_root
    String git_url
    File bacsum
    File arcsum
    File? short
    File? low
    File? unbinned
    File? checkm
    Array[File] hqmq_bin_fasta_files
    Array[File] bin_fasta_files
    Array[File] hqmq_bin_tarfiles
    File stats_json
    File stats_tsv
    Int n_hqmq=length(hqmq_bin_tarfiles)
    Int n_bin=length(bin_fasta_files)

    command {
        set -e
        end=`date --iso-8601=seconds`

        ln ${low} ${prefix}_bins.lowDepth.fa
        ln ${short} ${prefix}_bins.tooShort.fa
        ln ${unbinned} ${prefix}_bins.unbinned.fa
        ln ${checkm} ${prefix}_checkm_qa.out
        ln ${stats_json} ${prefix}_mags_stats.json
        ln ${bacsum} ${prefix}_gtdbtk.bac122.summary.tsv
        ln ${arcsum} ${prefix}_gtdbtk.ar122.summary.tsv

        # cp all tarfiles, zip them under prefix, if empty touch no_mags.txt
        mkdir -p hqmq
        if [ ${n_hqmq} -gt 0 ] ; then
            (cd hqmq && cp ${sep=" " hqmq_bin_tarfiles} .)
            (zip ../${prefix}_hqmq_bin.zip *tar.gz)
        else
            (cd hqmq && touch no_hqmq_mags.txt)
            (cd hqmq && zip ../${prefix}_hqmq_bin.zip *.txt)
        fi

        # Fix up attribute name
        cat ${stats_json} | \
           sed 's/lowDepth_/low_depth_/' > stats.json

        /scripts/generate_object_json.py \
                 --type "nmdc:MagsAnalysisActivity" \
                 --set mags_activity_set \
                 --part ${proj} \
                 -p "name=MAGS Activity for ${proj}" \
                    was_informed_by=${informed_by} \
                    started_at_time=${start} \
                    ended_at_time=$end \
                    execution_resource=${resource} \
                    git_url=${git_url} \
                    version="v1.0.4-beta" \
                 --url ${url_root}${proj}/mags/ \
                 --extra ./stats_json \
                 --inputs ${contigs} \
                        ${anno_gff} \
                        ${sorted_bam} \
                 --outputs \
                ${prefix}_checkm_qa.out "CheckM statistics report" "CheckM Statistics" "CheckM for ${proj}" \
                ${prefix}_hqmq_bin.zip "Metagenome bin tarfiles archive" "Metagenome Bins" "Metagenome Bins for ${proj}" \
                ${prefix}_gtdbtk.bac122.summary.tsv "GTDBTK bacterial summary" "GTDBTK Bacterial Summary" "Bacterial Summary for ${proj}" \
                ${prefix}_gtdbtk.ar122.summary.tsv "GTDBTK archaeal summary" "GTDBTK Archaeal Summary" "Archael Summary for ${proj}"

    }

    output {
        File objects = "objects.json"
        File final_checkm = "${prefix}_checkm_qa.out"
        File final_hqmq_bins_zip = "${prefix}_hqmq_bin.zip"
        File final_stats_json = "${prefix}_mags_stats.json"
        File final_gtdbtk_bac_summary = "${prefix}_gtdbtk.bac122.summary.tsv"
        File final_gtdbtk_ar_summary = "${prefix}_gtdbtk.ar122.summary.tsv"
        File final_lowDepth_fa = "${prefix}_bins.lowDepth.fa"
        File final_unbinned_fa = "${prefix}_bins.unbinned.fa"
        File final_short = "${prefix}_bins.tooShort.fa"
    }

    runtime {
        memory: "10 GiB"
        cpu:  4
        maxRetries: 1
        docker: container
    }
}
