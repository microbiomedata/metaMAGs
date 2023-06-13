#! /usr/bin/env python3
import os, sys, re, subprocess, sqlite3
from Bio import SeqIO

def checkm(clean_bins_dir, checkm_dir, checkm_qa_out, numcores, pplacer_cores):

 cl_rcode = subprocess.check_call(["checkm", "lineage_wf", "--pplacer_threads", str(pplacer_cores) ,"-x", "fa", "-t", str(numcores), clean_bins_dir, checkm_dir])
 if cl_rcode !=0 : sys.exit('Checkm lineage_wf failed. Please check\n')

 cq_rcode = subprocess.check_call(["checkm","qa", "-f", checkm_qa_out, checkm_dir + "/lineage.ms", checkm_dir])
 if cq_rcode !=0 : sys.exit('Checkm qa failed. Please check\n')

def gtdbtk(bins_dir, godir, sdb, numcores, pplacer_threads, log):
 
 rcode = subprocess.check_call(["gtdbtk", "classify_wf", "--pplacer_cpus", str(pplacer_threads), "--cpus", str(numcores), "--genome_dir", bins_dir, "--out_dir", godir, "--extension", "fa"])
 if rcode !=0 : sys.exit('GTDBK-TK failed.please check\n')

 conn = sqlite3.connect(sdb)
 c = conn.cursor()
 gtdbouts = [godir + '/gtdbtk.bac120.summary.tsv', godir + '/gtdbtk.ar122.summary.tsv']
 for gtdbout in gtdbouts:
  if os.path.isfile(gtdbout):
   with open(gtdbout) as rh:
    for line in rh:
     if line.startswith('user'): continue
     vals = line.rstrip().split('\t')
     binid = vals[0]

     lineage_arr = []
     long_lineage = vals[1].split(';') 
     for elineage in long_lineage:
      level,name = elineage.split('__')
      if name: lineage_arr.append(name)
     
     cmd = None
     if len(lineage_arr)>0 : cmd = "gtdbtk_domain = \"" + lineage_arr[0] + "\""
     if len(lineage_arr)>1 : cmd+= ",gtdbtk_phylum = \"" + lineage_arr[1] + "\""
     if len(lineage_arr)>2: cmd+= ",gtdbtk_class = \"" + lineage_arr[2] + "\""
     if len(lineage_arr)>3: cmd+= ",gtdbtk_order = \"" + lineage_arr[3] + "\""
     if len(lineage_arr)>4: cmd+= ",gtdbtk_family = \"" + lineage_arr[4] + "\""
     if len(lineage_arr)>5: cmd+= ",gtdbtk_genus = \"" + lineage_arr[5] + "\""
     if len(lineage_arr)>6: cmd+= ",gtdbtk_species = \"" + lineage_arr[6] + "\""

     if cmd:
      usql = 'update bin set ' + cmd + ' where bin_name = \'' + str(binid) + '\''
      log.debug(usql)
      c.execute(usql)

 conn.commit()
 conn.close() 

def phylodist(lineage_sdb, bins_dir, sdb, log):
 #based on imachen's phylodist script
 bininfo = dict()
 conn = sqlite3.connect(lineage_sdb)
 c = conn.cursor()

 for efile in os.listdir(bins_dir):
  lineage_h = {}
  epath = bins_dir + '/' + efile
  log.debug("phylodist on file:" + epath)
  scaffs = []
  for record in SeqIO.parse(epath, "fasta"):
   scaffs.append(record.id)

  for escaff in scaffs:
   query = "select lineage from contig_lin where scaffold_oid = '" + escaff + "'" 
   c.execute(query)
   for result in c:
    if len(result) > 0:
     scaf_lineage = result[0]
     lineage_arr = scaf_lineage.split(";")
     lineage_str = ""
     for s3 in lineage_arr:
      if len(lineage_str) > 0: lineage_str += ";" + s3
      else: lineage_str = s3
      if lineage_str in lineage_h: lineage_h[lineage_str] += 1
      else: lineage_h[lineage_str] = 1
 
  bin_lineage = ""
  ## threshold is at least 51% of contigs
  lineage_count = len(scaffs)
  threshold = len(scaffs) * 0.51
  ## however, if there are less than 10 contigs, threshold is 3
  if lineage_count < 10: threshold = 3
  ## however, if there is only 1 contig, threshold is 1
  if lineage_count == 1 : threshold = 1  

  for s in lineage_h:
   s_count = lineage_h[s]
   if s_count >= threshold:
    if s_count < lineage_count:
     bin_lineage = s
     lineage_count = s_count
    elif s_count == lineage_count and len(re.findall(';',s)) > len(re.findall(';',bin_lineage)):
     bin_lineage = s

  bin_lineage = bin_lineage.replace("'", "''")
  bininfo[efile] = bin_lineage

 c.close()
 conn.close()

 sconn = sqlite3.connect(sdb)
 sc = sconn.cursor()
 for ebin in bininfo:
  binname = re.sub(r'\.fa','',ebin)
  lvals= bininfo[ebin].split(';')
  uvals = ''
  if len(lvals)>0 : uvals += 'bin_domain = \'' + str(lvals[0]) + '\''
  if len(lvals)>1 : uvals += ',bin_phylum = \'' + str(lvals[1]) + '\''
  if len(lvals)>2 : uvals += ',bin_class = \'' + str(lvals[2]) + '\''
  if len(lvals)>3 : uvals += ',bin_order = \'' + str(lvals[3]) + '\''
  if len(lvals)>4 : uvals += ',bin_family = \'' + str(lvals[4]) + '\''
  if len(lvals)>5 : uvals += ',bin_genus = \'' + str(lvals[5]) + '\''
  if len(lvals)>6 : uvals += ',bin_species = \'' + str(lvals[6]) + '\''

  isql = 'update bin set ' + uvals + ' where bin_name = \'' + binname + '\''
  sc.execute(isql)

 sconn.commit()
 sconn.close()
