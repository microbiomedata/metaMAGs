#! /usr/bin/env python3
#njvarghese
import os, sys, sqlite3, subprocess
from Bio import SeqIO

def run(sdb, bins_dir, out_dir, log):
 if not os.path.isdir(out_dir):
  log.info("Specfied out dir does not exist. Exit module.")
  return()

 #clean metabat bins
 conn = sqlite3.connect(sdb)
 c = conn.cursor()

 #per bin modify or copy-----------------------------
 for efile in os.listdir(bins_dir):
  domains = dict()
  epath = bins_dir + '/' + efile
  log.debug(epath)
 
  scaffs = dict()
  for record in SeqIO.parse(epath, "fasta"):
   scaffs[record.id] = record.seq

  query = "select scaf_domain, scaffold_id from bin_scaffolds where scaffold_id in (\"" + "\",\"".join(scaffs.keys()) + "\")"
  c.execute(query)
  for result in c:
   domain, escaff = result
   if not domain : domain = 'unknown'
   if domain in domains : domains[domain].append(escaff)
   else: domains[domain] = [escaff]

  #check if all domains are same
  reqlen = 1
  if 'unknown' in domains : reqlen = 2

  if len(domains.keys()) > reqlen:

   #get domain with longest set of values
   bestdomain = None
   for edomain in domains:
    if edomain != 'unknown' :
     if not bestdomain : bestdomain = edomain
     elif len(domains[edomain]) > len(domains[bestdomain]) : bestdomain = edomain
   log.info('modifying ' + efile + " to retain 'unknown' and majority domain:" + bestdomain)

   #remove scaffolds that do not match the best domain
   bad_soids = list(scaffs)
   ofn = out_dir + '/' + efile
   wh = open(ofn,'w')
   for esoid in domains[bestdomain]:
    output = '>' + str(esoid) + '\n' + str(scaffs[esoid]) + '\n'
    wh.write(output)
    bad_soids.remove(esoid)
   if 'unknown' in domains:
    for esoid in domains['unknown']:
     output = '>' + str(esoid) + '\n' + str(scaffs[esoid]) + '\n'
     wh.write(output)
     bad_soids.remove(esoid)

   wh.close()
   log.debug(ofn)
 
   #delete bad soids from sqlite db
   r = conn.cursor()
   rsql = 'delete from bin_scaffolds where scaffold_id in (\"' + "\",\"".join(bad_soids) + "\")" 
   sys.stderr.write(rsql + '\n')
   r.execute(rsql)
   conn.commit()

  else:
   npath = out_dir + '/' + efile
   subprocess.check_output(['cp',epath,npath])
   log.debug(npath)

 conn.close()

