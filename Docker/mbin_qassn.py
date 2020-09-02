#! /usr/bin/env python
import os,sys, re
import sqlite3, subprocess
import collections
from collections import defaultdict
#njvarghese 2020

def run(checkm_qa_out, bins_dir, hqmq_bins_dir, sdb_name):

 conn = sqlite3.connect(sdb_name)
 c = conn.cursor()
 bmetadata = collections.defaultdict()

 with open(checkm_qa_out) as rh:
  for line in rh:
   if line.startswith('-'): continue
   if re.search(r'#',line): continue
   evals = re.split('\s{2,}',line.rstrip())
   binid = evals[1]
   comp = float(evals[12])
   cont = float(evals[13])

   if binid not in bmetadata: bmetadata[binid] = dict()
   bmetadata[binid]['completeness'] = comp
   bmetadata[binid]['contamination'] = cont

   if comp > 90 and cont < 5 :
    qsql = 'select num_16s, num_5s, num_23s, num_tRNA from bin where bin_name = "' + binid + '"'
    sys.stderr.write(qsql + '\n')
    c.execute(qsql)
    for result in c:
     num_16s, num_5s, num_23s, num_tRNA = result
     if num_16s == 0 or num_5s == 0 or num_23s == 0 or num_tRNA < 18 : bmetadata[binid]['quality'] = 'MQ'
     else: bmetadata[binid]['quality'] = 'HQ'
   elif comp >=50 and cont < 10 : bmetadata[binid]['quality'] = 'MQ'
   else : bmetadata[binid]['quality'] = 'LQ'

 i = conn.cursor() #insert cursor
 if not os.path.isdir(hqmq_bins_dir): subprocess.check_call(['mkdir',hqmq_bins_dir])

 for binid in bmetadata:
  #update sqlite table
  usql = 'update bin set completeness = ' + str(bmetadata[binid]['completeness']) + ', contamination = ' + str(bmetadata[binid]['contamination']) + ', bin_quality = "' + str(bmetadata[binid]['quality']) + '" where bin_name = "' + str(binid) + '"'
  #print usql
  i.execute(usql)

  #move to separate directory
  if bmetadata[binid]['quality'] == 'HQ' or bmetadata[binid]['quality'] == 'MQ':
   fpath = bins_dir+'/'+binid+'.fa'
   #sys.stderr.write(fpath + '\n')
   subprocess.check_call(['mv', fpath, hqmq_bins_dir])

 conn.commit()
 conn.close()
