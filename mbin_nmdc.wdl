workflow nmdc_mags 
    String? outdir
    String  proj_name
    File contig_file
    File sam_file
    File gff_file
    String container = "microbiomedata/nmdc_mbin:0.1.1"
    File? map_file
    File? domain_file
    Int cpu=32
    Int pplacer_cpu=1
    
    ## Need have db name with checkM_DB and GTDBTK_DB in the database directory
    String database="/refdata"

    call mbin_nmdc{
         input:  name=proj_name, 
                 fasta=contig_file, 
                 sam=sam_file, 
                 gff=gff_file, 
                 map=map_file, 
                 domain=domain_file, 
		 cpu=cpu,
                 pplacer_cpu=pplacer_cpu,
                 database=database,
	         container=container
    }
  
    call make_output {
       	input: outdir= outdir, mbin_nmdc_output=mbin_nmdc.runScript
    }

    output {
	Array[File?] final_hqmq_bins = make_output.hqmq_bin_fasta_files
	Array[File?] metabat_bins = make_output.metabat_bin_fasta_files
	File? final_checkm = make_output.checkm_output
	File? final_gtdbtk_bac_summary = make_output.gtdbtk_bac_summary
	File? final_gtdbtk_ar_summary = make_output.gtdbtk_ar_summary
	File? final_tooShort_fa = make_output.tooShort_fa
	File? final_lowDepth_fa = make_output.lowDepth_fa
	File? final_unbinned_fa = make_output.unbinned_fa
    }
    parameter_meta {
	cpu: "number of CPUs"
        pplacer_cpu: "number of threads used by pplacer"
	outdir: "the final output directory path"
	proj_name: "project name"
	contig_file: "input assembled contig fasta file"
	sam_file: "Sam/Bam file from reads mapping back to contigs. [sam.gz or bam]"
	gff_file: "contigs functional annotation result in gff format"
	map_file: "text file which containing mapping of headers between SAM and FNA (ID in SAM/FNA ID in GFF)"
	database: "database directory path which includes checkM_DB and GTDBTK_DB subdirectories"
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
        String container
	String filename_outlog="stdout.log"
	String filename_errlog="stderr.log"
	String filename_stat="checkm_qa.out"
	String dollar="$"
	runtime {
                docker: container
		mem: "120 GiB"
		cpu:  cpu
		database: database
 	}   

     command {
	export TIME="time result\ncmd:%C\nreal %es\nuser %Us \nsys  %Ss \nmemory:%MKB \ncpu %P"
	set -eo pipefail
	# set TMPDIR to avoid AF_UNIX path too long error 
	export TMPDIR=/tmp
	mbin_nmdc.py ${"--map " + map} ${"--domain " + domain} --pplacer_cpu ${pplacer_cpu} --cpu ${cpu} ${name} ${fasta} ${sam} ${gff}
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
 
 	command{
		mbin_nmdc_path=`dirname ${mbin_nmdc_output}`
		mkdir -p ${outdir}
		mv -f $mbin_nmdc_path/* ${outdir}/
 		chmod 764 -R ${outdir}
 	}
	output {
		Array [String] hqmq_bin_fasta_files = glob("${outdir}/hqmq-metabat-bins/*fa")
		Array [String] metabat_bin_fasta_files = glob("${outdir}/metabat-bins/*fa")
		String checkm_output = "${outdir}/checkm_qa.out"
		String gtdbtk_bac_summary = "${outdir}/gtdbtk_output/gtdbtk.bac120.summary.tsv"
		String gtdbtk_ar_summary = "${outdir}/gtdbtk_output/gtdbtk.ar122.summary.tsv"
		String unbinned_fa = "${outdir}/gtdbtk_output/bins.unbinned.fa"
		String tooShort_fa = "${outdir}/gtdbtk_output/bins.tooShort.fa
		String lowDepth_fa = "${outdir}/gtdbtk_output/bins.lowDepth.fa"
	}
	runtime {
            mem: "1 GiB"
            cpu:  1
        }
}

