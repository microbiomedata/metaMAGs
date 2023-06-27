workflow mbin{
    
    String informed_by
    String resource
    String url_root
    String git_url
    String proj_name
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
    File? map_file
    File? domain_file
    String? scratch_dir
    Int cpu=32
    Int threads=64
    Int pthreads=1
    String gtdbtk_db="/refdata/GTDBTK_DB/gtdbtk_release207_v2"
    String checkm_db="/refdata/CheckM_DB/checkm_data_2015_01_16"
    String container = "microbiomedata/nmdc_mbin@sha256:7587e3b777fa6a37c7355e56300a689776235e86494be0ab454e0cb789594765"

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

    call mbin_nmdc {
        input:  
                name=proj_name,
                fna = stage.contig,
                aln = stage.sam,
                gff = stage.gff,
                #map = stage.map_file,
                threads =  threads,
                pthreads = pthreads,
                gtdbtk_env = gtdbtk_db,
                checkm_env = checkm_db,
                mbin_container = container
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
                 product_names_file=stage.product_names,
                 container=container
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
        hqmq_bin_fasta_files = mbin_nmdc.hqmq_bin_fasta_files,
        bin_fasta_files = mbin_nmdc.bin_fasta_files,
        hqmq_bin_tarfiles = package.hqmq_bin_tarfiles

    }

    output {
        File final_hqmq_bins_zip = finish_mags.final_hqmq_bins_zip
        File final_gtdbtk_bac_summary = finish_mags.final_gtdbtk_bac_summary
        File final_gtdbtk_ar_summary = finish_mags.final_gtdbtk_ar_summary
        File short = finish_mags.final_short
        File low = finish_mags.final_lowDepth_fa
        File final_unbinned_fa  = finish_mags.final_unbinned_fa
        File final_checkm = finish_mags.final_checkm
        File mags_objects = finish_mags.objects
        File mags_version = finish_mags.final_version
    }


}

task mbin_nmdc {
    
    File fna
    File aln
    File gff
    String name
    Int? threads
    Int? pthreads
    String gtdbtk_env
    String checkm_env
    String mbin_container
    

    command<<<
        export GTDBTK_DATA_PATH=${gtdbtk_env}
        export CHECKM_DATA_PATH=${checkm_env}
        mbin.py ${"--threads " + threads} ${"--pthreads " + pthreads} --fna ${fna} --gff ${gff} --aln ${aln}
        mbin_stats.py $PWD
        mbin_versions.py > mbin_nmdc_versions.log
        touch MAGs_stats.tsv
    
        if [ -f  gtdbtk-output/gtdbtk.bac120.summary.tsv ]; then
            echo "bacterial summary exists."
        else
            echo "No Bacterial Results for ${name}" > gtdbtk_output/gtdbtk.bac120.summary.tsv
        fi

        if [ -f  gtdbtk-output/gtdbtk.ar122.summary.tsv ]; then
            echo "archaeal summary exists."
        else
            echo "No Archaeal Results for ${name}" > gtdbtk_output/gtdbtk.ar122.summary.tsv
        fi
    >>>

    runtime{
        docker : mbin_container
        memory : "60 G"
	    time : "2:00:00"
        cpu : threads
    }

    output{
        File short = "bins.tooShort.fa"
        File low = "bins.lowDepth.fa"
        File unbinned = "bins.unbinned.fa"
        File checkm = "checkm-qa.out"
        File stats_json = "MAGs_stats.json"
        File stats_tsv = "MAGs_stats.tsv"
        File mbin_version = "mbin_nmdc_version.log"
        File bacsum = "gtdbtk-output/gtdbtk.bac120.summary.tsv"
        File arcsum = "gtdbtk-output/gtdbtk.ar122.summary.tsv"
        Array[File] hqmq_bin_fasta_files = glob("hqmq-metabat-bins/*fa")
        Array[File] bin_fasta_files = glob("metabat-bins/*fa")
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
     String container 

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
    File mbin_version
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
        ln ${mbin_version} ${prefix}_bin.info
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
           sed 's/: null/: "null"/g' | \
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
                 --extra ./stats.json \
                 --inputs ${contigs} \
                        ${anno_gff} \
                        ${sorted_bam} \
                 --outputs \
                ${prefix}_checkm_qa.out "CheckM statistics report" "CheckM Statistics" "CheckM for ${proj}" \
                ${prefix}_hqmq_bin.zip "Metagenome bin tarfiles archive" "Metagenome Bins" "Metagenome Bins for ${proj}" \
                ${prefix}_gtdbtk.bac122.summary.tsv "GTDBTK bacterial summary" "GTDBTK Bacterial Summary" "Bacterial Summary for ${proj}" \
                ${prefix}_gtdbtk.ar122.summary.tsv "GTDBTK archaeal summary" "GTDBTK Archaeal Summary" "Archaeal Summary for ${proj}" \
                ${prefix}_bin.info "File containing version information on the binning workflow" "Metagenome Bins Info File" "Metagenome Bins Info File for ${proj}"

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
        File final_version = "${prefix}_bin.info"
    }

    runtime {
        memory: "10 GiB"
        cpu:  4
        maxRetries: 1
        docker: container
    }
}
