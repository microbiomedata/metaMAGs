#! /usr/bin/env python3
#Author: njvarghese|2024

import os, sys, re, argparse, subprocess
from ete3 import NCBITaxa
import logfile

def run(bins_dir, eukcc_db, num_threads, log):
 
 odir = 'eukcc_output'
 ofile = odir + '/' + 'eukcc.csv'
 final_ofile = odir + '/' + 'eukcc.csv.final'
 eukcc_log = odir + '/' + 'eukcc.log'
 ncbi_db = eukcc_db + '/etetoolkit/taxa.sqlite'
 
 #check paths and dbs
 log.info("Checking input files and databases.")
 if not os.path.isdir(eukcc_db): 
  log.info("ERROR. EukCC dabatabase path is incorrect. Please check.")
  exit()

 if not os.path.isfile(ncbi_db):
  log.info("ERROR. NCBI taxa.sqlite missing from dir etetoolkit in EukCC DB.")
  exit()

 if not os.path.isdir(bins_dir):
  log.info("ERROR. Input bins directory path is incorrect. Please check.")
  exit()

 #check if bins dir is empty
 input_bins = os.listdir(bins_dir)
 if len(input_bins) < 1 : 
  log.info("ERROR. Bins directory is empty. Please check.")
  exit()
 else:
  log.info("Found " + str(len(input_bins)) + " bins to process.")
 
 #mkdir
 if not os.path.isdir(odir): os.mkdir(odir)

 #run eukcc
 log.info("Running EukCC.")
 with open(eukcc_log, 'w') as fh:
  try:
   # Run the subprocess and capture return code
   subprocess.check_call(
            ['eukcc', 'folder', '--out', odir, '--threads', str(num_threads), '--db', eukcc_db, bins_dir],
            stdout=fh, stderr=subprocess.STDOUT  # Redirect stderr to stdout to capture all logs
   )
   log.info("EukCC completed successfully.")
  except subprocess.CalledProcessError as e:
   log.error(f"EukCC failed with return code {e.returncode}. Check the log file for details.")
   sys.exit('EukCC failed. Please check\n')
  except Exception as e:
   log.error(f"An unexpected error occurred while running EukCC: {str(e)}")
   sys.exit('EukCC failed. Please check\n')

 #run taxonomy conversion
 ouput_with_lineage = taxlookup(ofile, ncbi_db)
 wh = open(final_ofile, 'w')
 wh.write(ouput_with_lineage)
 wh.close()

 log.info("EukCC output: " + ofile)
 log.info("EBIP output: " + final_ofile)

def taxlookup(eukcc_csv, ncbi_db):
 ncbi = NCBITaxa(dbfile = ncbi_db)
 ovals = "bin" + "\t" + "completeness" + "\t" + "contamination" + "\t" + "ncbi_lineage_taxIDs" + "\t" + "ncbi_lineage" + '\n'

 with open(eukcc_csv) as fh:
  for line in fh:
   binid,comp,cont,lineage = line.rstrip().split('\t')

   if re.search('ncbi',lineage) : continue
   if re.match('NA',lineage): continue

   splits = lineage.split('-')
   #print(splits)

   lineage_dict = ncbi.get_taxid_translator(splits)
   lineage_list = []
   for each_split in splits:
    lineage_list.append(lineage_dict[int(each_split)])

   ovals += binid + '\t' + comp + '\t' + cont + '\t' + lineage + '\t' + ','.join(lineage_list) + '\n'

 return(ovals)

if __name__ == '__main__':
 parser = argparse.ArgumentParser(
  prog= "mbin_ebip.py" ,
  description='Eukaryotic bin identification pipeline. Author: Neha J. Varghese.')

 parser.add_argument("--threads", type=int, default=1,  help = "Num of threads/cores (default:1)")
 requiredNamed = parser.add_argument_group('Required arguments')
 requiredNamed.add_argument("--bins", type=str, help= "Directory with bins in fasta format", required=True)
 requiredNamed.add_argument("--db", type=str, help= "EukCC database path", required=True)
 
 args = parser.parse_args()
 log = logfile.startlog('ebip-docker')
 run(args.bins, args.db, args.threads, log)
