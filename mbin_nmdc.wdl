version 1.0
workflow nmdc_mags {
    input {
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
        String lineage_file
        File? map_file
        File? domain_file
        String? scratch_dir
        Int cpu=32
        Int threads=64
        Int pthreads=1
        String gtdbtk_db="/refdata/GTDBTK_DB/gtdbtk_release207_v2"
        String checkm_db="/refdata/checkM_DB/checkm_data_2015_01_16"
        String package_container = "microbiomedata/nmdc_mbin_vis:0.1"
        String container = "microbiomedata/nmdc_mbin@sha256:c8df293e80698627ce66df7cd07f6b10e9112184e3bf1379e615d10123f7bc64"
    }
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
            gene_phylogeny_file=gene_phylogeny_file,
            lineage_file=lineage_file
    }

    call mbin_nmdc {
        input:  
                name=proj_name,
                fna = stage.contig,
                aln = stage.sam,
                gff = stage.gff,
                lineage=stage.lineage_tsv,
                threads =  threads,
                pthreads = pthreads,
                gtdbtk_env = gtdbtk_db,
                checkm_env = checkm_db,
                mbin_container = container
    }
    call package {
         input:  proj = proj_name,
                 bins=flatten([mbin_nmdc.hqmq_bin_fasta_files,mbin_nmdc.lq_bin_fasta_files]),
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
                 container=package_container
    }

    call finish_mags {
        input:
        container="microbiomedata/workflowmeta:1.1.1",
        contigs=stage.contig,
        anno_gff=stage.gff,
        sorted_bam=stage.sam,
        proj=proj_name,
        start=stage.start,
        checkm = mbin_nmdc.checkm,
        bacsum= mbin_nmdc.bacsum,
        arcsum = mbin_nmdc.arcsum,
        short = mbin_nmdc.short,
        low = mbin_nmdc.low,
        unbinned = mbin_nmdc.unbinned,
        checkm = mbin_nmdc.checkm,
        mbin_sdb = mbin_nmdc.mbin_sdb,
        mbin_version = mbin_nmdc.mbin_version,
        stats_json = mbin_nmdc.stats_json,
        stats_tsv = mbin_nmdc.stats_tsv,
        hqmq_bin_fasta_files = mbin_nmdc.hqmq_bin_fasta_files,
        bin_fasta_files = mbin_nmdc.lq_bin_fasta_files,
        hqmq_bin_tarfiles = package.hqmq_bin_tarfiles,
        lq_bin_tarfiles = package.lq_bin_tarfiles,
        barplot = package.barplot,
        heatmap = package.heatmap,
        kronaplot = package.kronaplot,
        ko_matrix = package.ko_matrix
    }

    output {
        File final_hqmq_bins_zip = finish_mags.final_hqmq_bins_zip
        File final_lq_bins_zip = finish_mags.final_lq_bins_zip
        File final_gtdbtk_bac_summary = finish_mags.final_gtdbtk_bac_summary
        File final_gtdbtk_ar_summary = finish_mags.final_gtdbtk_ar_summary
        File short = finish_mags.final_short
        File low = finish_mags.final_lowDepth_fa
        File final_unbinned_fa  = finish_mags.final_unbinned_fa
        File final_checkm = finish_mags.final_checkm
        File mags_version = finish_mags.final_version
        File final_stats_json = finish_mags.final_stats_json
        File barplot = finish_mags.final_barplot
        File heatmap = finish_mags.final_heatmap
        File kronaplot = finish_mags.final_kronaplot
    }


}

task mbin_nmdc {
    input{
        File fna
        File aln
        File gff
        File lineage
        String name
        Int? threads
        Int? pthreads
        String gtdbtk_env
        String checkm_env
        String mbin_container
    }

    command<<<
        set -euo pipefail
        export GTDBTK_DATA_PATH=~{gtdbtk_env}
        export CHECKM_DATA_PATH=~{checkm_env}
        mbin.py ~{"--threads " + threads} ~{"--pthreads " + pthreads} --fna ~{fna} --gff ~{gff} --aln ~{aln} --lintsv ~{lineage}
        mbin_stats.py $PWD
        mbin_versions.py > mbin_nmdc_versions.log
        touch MAGs_stats.tsv
    
        if [ -f  gtdbtk-output/gtdbtk.bac120.summary.tsv ]; then
            echo "bacterial summary exists."
        else
            mkdir -p gtdbtk-output
            echo "No Bacterial Results for ~{name}" > gtdbtk-output/gtdbtk.bac120.summary.tsv
        fi

        if [ -f  gtdbtk-output/gtdbtk.ar122.summary.tsv ]; then
            echo "archaeal summary exists."
        else
            mkdir -p gtdbtk-output
            echo "No Archaeal Results for ~{name}" > gtdbtk-output/gtdbtk.ar122.summary.tsv
        fi

        if [ -f checkm-qa.out ]; then
            echo "checkm summary exists."
        else
            mkdir -p gtdbtk-output
            echo "No Checkm Results for ~{name}" > checkm-qa.out
        fi

        if [ -f mbin.sdb ]; then
            echo "mbin.sdb exists."
        else
            mkdir -p gtdbtk-output
            echo "Mbin Sdb Could not be created for ~{name}" > mbin.sdb
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
        File mbin_sdb = "mbin.sdb"
        File mbin_version = "mbin_nmdc_versions.log"
        File bacsum = "gtdbtk-output/gtdbtk.bac120.summary.tsv"
        File arcsum = "gtdbtk-output/gtdbtk.ar122.summary.tsv"
        Array[File] hqmq_bin_fasta_files = glob("hqmq-metabat-bins/*fa")
        Array[File] lq_bin_fasta_files = glob("filtered-metabat-bins/*fa")
    }    
}


task stage {
    input{
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
        String lineage_file
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
        String lineage_out="lineage.tsv"
    }
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

        stage ~{contig_file} ~{contigs_out}
        stage ~{sam_file} ~{bam_out}
        stage ~{gff_file} ~{gff_out}
        stage ~{proteins_file} ~{proteins_out}
        stage ~{cog_file} ~{cog_out}
        stage ~{ec_file} ~{ec_out}
        stage ~{ko_file} ~{ko_out}
        stage ~{pfam_file} ~{pfam_out}
        stage ~{tigrfam_file} ~{tigrfam_out}
        stage ~{cath_funfam_file} ~{cath_funfam_out}
        stage ~{smart_file} ~{smart_out}
        stage ~{supfam_file} ~{supfam_out}
        stage ~{product_names_file} ~{products_out}
        stage ~{gene_phylogeny_file} ~{gene_phylogeny_out}
        stage ~{lineage_file} ~{lineage_out}

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
        File lineage_tsv = "lineage.tsv"
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
    input{
        String proj
        String prefix=sub(proj, ":", "_")
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
    }
     command<<<
         set -e
         create_tarfiles.py ~{prefix} \
                     ~{json_stats} ~{gff_file} ~{proteins_file} ~{cog_file} \
                     ~{ec_file} ~{ko_file} ~{pfam_file} ~{tigrfam_file} \
                     ~{cath_funfam_file} ~{smart_file} ~{supfam_file} \
                     ~{product_names_file} \
                     ~{sep=" " bins}

        if [ -f ~{prefix}_heatmap.pdf ]; then
            echo "KO analysis plot exists."
        else
            echo "No KO analysis result for ~{proj}" > ${prefix}_heatmap.pdf
            echo "No KO analysis result for ~{proj}" > ${prefix}_barplot.pdf
            echo "No KO analysis result for ~{proj}" > ${prefix}_ko_krona.html
            echo "No KO analysis result for ~{proj}" > ${prefix}_module_completeness.tab
        fi
     >>>
     output {
         Array[File] hqmq_bin_tarfiles = flatten([glob("*_HQ.tar.gz"), glob("*_MQ.tar.gz")])
         Array[File] lq_bin_tarfiles = glob("*_LQ.tar.gz") 
         File barplot = prefix + "_barplot.pdf"
         File heatmap = prefix + "_heatmap.pdf"
         File kronaplot = prefix + "_ko_krona.html"
         File ko_matrix = prefix + "_module_completeness.tab"
     }
     runtime {
         docker: container
         memory: "1 GiB"
         cpu:  1
     }
}

task finish_mags {
    input{
        String container
        File contigs
        File anno_gff
        File sorted_bam
        File mbin_sdb
        File mbin_version
        String proj
        String prefix=sub(proj, ":", "_")
        String start
        File bacsum
        File arcsum
        File? short
        File? low
        File? unbinned
        File? checkm
        Array[File] hqmq_bin_fasta_files
        Array[File] bin_fasta_files
        Array[File] hqmq_bin_tarfiles
        Array[File] lq_bin_tarfiles
        File stats_json
        File stats_tsv
        Int n_hqmq=length(hqmq_bin_tarfiles)
        Int n_bin=length(bin_fasta_files)
        File barplot
        File heatmap
        File kronaplot
        File ko_matrix
    }
    command<<<
        set -e
        end=`date --iso-8601=seconds`

        ln ~{low} ~{prefix}_bins.lowDepth.fa
        ln ~{short} ~{prefix}_bins.tooShort.fa
        ln ~{unbinned} ~{prefix}_bins.unbinned.fa
        ln ~{checkm} ~{prefix}_checkm_qa.out
        ln ~{mbin_version} ~{prefix}_bin.info
        ln ~{bacsum} ~{prefix}_gtdbtk.bac122.summary.tsv
        ln ~{arcsum} ~{prefix}_gtdbtk.ar122.summary.tsv
        ln ~{barplot} ~{prefix}_barplot.pdf
        ln ~{heatmap} ~{prefix}_heatmap.pdf
        ln ~{kronaplot} ~{prefix}_kronaplot.html
        ln ~{ko_matrix} ~{prefix}_ko_matrix.txt

        # cp all tarfiles, zip them under prefix, if empty touch no_mags.txt
        mkdir -p hqmq
        if [ ~{n_hqmq} -gt 0 ] ; then
            (cd hqmq && cp ~{sep=" " hqmq_bin_tarfiles} .)
            (cd hqmq && cp ~{mbin_sdb} .)
            (cd hqmq && zip -j ../~{prefix}_hqmq_bin.zip *tar.gz mbin.sdb ../*pdf ../*kronaplot.html ../*ko_matrix.txt)
        else
            (cd hqmq && touch no_hqmq_mags.txt)
            (cd hqmq && cp ~{mbin_sdb} .)
            (cd hqmq && zip ../~{prefix}_hqmq_bin.zip *.txt mbin.sdb)
        fi

        mkdir -p lq
        if [ ~{n_bin} -gt 0 ] ; then
            (cd lq && cp ~{sep=" " lq_bin_tarfiles} .)
            (cd lq && cp ~{mbin_sdb} .)
            (cd lq && zip -j ../~{prefix}_lq_bin.zip *tar.gz mbin.sdb ../*pdf ../*kronaplot.html ../*ko_matrix.txt)
        else
            (cd lq && touch no_lq_mags.txt)
            (cd lq && cp ~{mbin_sdb} .)
            (cd lq && zip ../~{prefix}_lq_bin.zip *.txt mbin.sdb)
        fi

        # Fix up attribute name
        cat ~{stats_json} | \
           sed 's/: null/: "null"/g' | \
           sed 's/lowDepth_/low_depth_/' > ~{prefix}_mags_stats.json

    >>>

    output {
        File final_checkm = "~{prefix}_checkm_qa.out"
        File final_hqmq_bins_zip = "~{prefix}_hqmq_bin.zip"
        File final_lq_bins_zip = "~{prefix}_lq_bin.zip"
        File final_stats_json = "~{prefix}_mags_stats.json"
        File final_gtdbtk_bac_summary = "~{prefix}_gtdbtk.bac122.summary.tsv"
        File final_gtdbtk_ar_summary = "~{prefix}_gtdbtk.ar122.summary.tsv"
        File final_lowDepth_fa = "~{prefix}_bins.lowDepth.fa"
        File final_unbinned_fa = "~{prefix}_bins.unbinned.fa"
        File final_short = "~{prefix}_bins.tooShort.fa"
        File final_version = "~{prefix}_bin.info"
        File final_kronaplot = "~{prefix}_kronaplot.html"
        File final_heatmap = "~{prefix}_heatmap.pdf"
        File final_barplot = "~{prefix}_barplot.pdf"
    }

    runtime {
        memory: "10 GiB"
        cpu:  4
        maxRetries: 1
        docker: container
    }
}
