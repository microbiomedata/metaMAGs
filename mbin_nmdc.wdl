workflow nmdc_mags {
    String? outdir
    String  proj_name
    File contig_file
    File sam_file
    File gff_file
    File? map_file
    File? domain_file
    Int cpu=32
    
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
                 database=database
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
    }
    parameter_meta {
	cpu: "number of CPUs"
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
    }
    meta {
        author: "Chienchi Lo, B10, LANL"
        email: "chienchi@lanl.gov"
        version: "1.0.0"
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
	String filename_outlog="stdout.log"
	String filename_errlog="stderr.log"
	String filename_stat="checkm_qa.out"
	String dollar="$"
	runtime {
		memory: "120 GiB"
		cpu:  cpu
		database: database
 	}   

     command {
	export TIME="time result\ncmd:%C\nreal %es\nuser %Us \nsys  %Ss \nmemory:%MKB \ncpu %P"
	set -eo pipefail
	# set TMPDIR to avoid AF_UNIX path too long error 
	export TMPDIR=/tmp
	mbin_nmdc.py ${"--map " + map} ${"--domain " + domain} --cpu ${cpu} ${name} ${fasta} ${sam} ${gff}
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
	}
	runtime {
            memory: "1 GiB"
            cpu:  1
        }
}

