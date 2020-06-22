workflow nmdc_mags {
    String? outdir
    String  proj_name
    File contig_file
    File sam_file
    File gff_file
    File? map_file
    File? domain_file
    Int cpu=32
    
    String mbin_container="microbiomedata/nmdc_mbin:0.1.0"
    ## Need have db name with checkM_DB and GTDBTK_DB in the database directory
    String database="/global/cfs/projectdirs/m3408/aim2/database"

    call mbin_nmdc{
         input:  name=proj_name, 
                 fasta=contig_file, 
                 sam=sam_file, 
                 gff=gff_file, 
                 map=map_file, 
                 domain=domain_file, 
		 cpu=cpu,
                 container=mbin_container, 
                 database=database
    }
  
    call make_output {
       	input: outdir= outdir, mbin_nmdc_output=mbin_nmdc.stat
    }

}

task mbin_nmdc {
	String name
	File fasta
	File sam
	File gff
	File? map
	File? domain
	String container
        String database
	Int cpu
	String filename_outlog="stdout.log"
	String filename_errlog="stderr.log"
	String filename_stat="checkm_qa.out"
	String dollar="$"
     runtime{ mem: "115GB"
             cpu: cpu
            time: "10:00:00"
        jobname: "nmdc_mags"
     }
    # runtime {
     #       backend: "Local"
      #      docker: container
       #     memory: "120 GiB"
	#    cpu:  16
         #   database: database
     #}

     command {
	export TIME="time result\ncmd:%C\nreal %es\nuser %Us \nsys  %Ss \nmemory:%MKB \ncpu %P"
	set -eo pipefail
	# set TMPDIR to avoid AF_UNIX path too long error
	export TMPDIR=/tmp
	shifter --image=${container} -V ${database}:/databases mbin_nmdc.py ${"--map " + map} ${"--domain " + domain} --cpu ${cpu} ${name} ${fasta} ${sam} ${gff}
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
	runtime{ mem: "2GB"
                cpu: 1
		time: "1:00:00"
		jobname: "mbin_nmdc_output"
	}
}

