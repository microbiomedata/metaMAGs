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
       	input: outdir= outdir, mbin_nmdc_output=mbin_nmdc.stat
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
     }
     output {
	#File stdout = filename_outlog
	#File stderr = filename_errlog
	File stat = filename_stat
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
}

