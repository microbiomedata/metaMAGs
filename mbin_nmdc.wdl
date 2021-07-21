workflow nmdc_mags {
    String? outdir
    String  proj_name
    File contig_file
    File sam_file
    File gff_file
    String container = "microbiomedata/nmdc_mbin:0.1.2"
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
  
    call make_output {
       	input: outdir= outdir, mbin_nmdc_output=mbin_nmdc.runScript, container=container
    }

    output {
	Array[File] final_hqmq_bins = make_output.hqmq_bin_fasta_files
	Array[File] metabat_bins = make_output.metabat_bin_fasta_files
	File? final_checkm = make_output.checkm_output
	File? final_gtdbtk_bac_summary = make_output.gtdbtk_bac_summary
	File? final_gtdbtk_ar_summary = make_output.gtdbtk_ar_summary
	File? final_tooShort_fa = make_output.tooShort_fa
	File? final_lowDepth_fa = make_output.lowDepth_fa
	File? final_unbinned_fa = make_output.unbinned_fa
	File? final_stats = make_output.json_stats
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
	# set TMPDIR to avoid AF_UNIX path too long error 
	export TMPDIR=/tmp
	export GTDBTK_DATA_PATH=/databases
	mbin_nmdc.py ${"--map " + map} ${"--domain " + domain} ${"--scratch_dir " + scratch_dir} --pplacer_cpu ${pplacer_cpu} --cpu ${cpu} ${name} ${fasta} ${sam} ${gff}
	mbin_stats.py $PWD
     }
     output {
	File runScript = "script"
	#File stderr = filename_errlog
	File? stat = filename_stat
     }
}

task make_output{
 	String outdir
	String mbin_nmdc_output
	String container
 
 	command{
		mbin_nmdc_path=`dirname ${mbin_nmdc_output}`
		mkdir -p ${outdir}
		mv -f $mbin_nmdc_path/* ${outdir}/
 		chmod 764 -R ${outdir}
 	}
	output {
		Array[File] hqmq_bin_fasta_files = glob("${outdir}/hqmq-metabat-bins/*fa")
		Array[File] metabat_bin_fasta_files = glob("${outdir}/metabat-bins/*fa")
		File? checkm_output = "${outdir}/checkm_qa.out"
		File? gtdbtk_bac_summary = "${outdir}/gtdbtk_output/gtdbtk.bac120.summary.tsv"
		File? gtdbtk_ar_summary = "${outdir}/gtdbtk_output/gtdbtk.ar122.summary.tsv"
		File? unbinned_fa = "${outdir}/bins.unbinned.fa"
		File? tooShort_fa = "${outdir}/bins.tooShort.fa"
		File? lowDepth_fa = "${outdir}/bins.lowDepth.fa"
		File? json_stats = "${outdir}/MAGs_stats.json"
	}
	runtime {
	    docker: container
            memory: "1 GiB"
            cpu:  1
        }
}

