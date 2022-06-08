#!/usr/bin/env python
import os,sys,sqlite3
import mbin_sdb

def query_slineage(sdb,domainfn):
 #check if file exists
 if os.path.isfile(domainfn): pass
 else:
  sys.stderr.write('Domain file does not exist. Skipping lineage query.\n')
  return()
 
 #parse file to extract domain informatin and add to sqlite database
 with open(domainfn) as rh:
  for line in rh:
   soid,domain = line.rstrip().split('\t')
   sql = "update bin_scaffolds set domain =" + str(domain)  + " where scaffold_id = \"" + str(esoid) + "\""
   mbin_sdb.update(sql)

def parse_annotation(osdb, gff):
 soids = set()
 isqls = []
 ngenes = dict()
 num_16s = dict()
 num_5s = dict()
 num_23s = dict()
 num_tRNA = dict()

 #connect to current sdb
 oconn = sqlite3.connect(osdb)
 oc = oconn.cursor()
 osql = "select bin_name,scaffold_id from bin_scaffolds"
 oc.execute(osql)
 for result in oc:
  bname, soid = result
  soids.add(soid)
  ngenes[soid] = 0
  num_16s[soid] = 0
  num_5s[soid] = 0
  num_23s[soid] = 0
  num_tRNA[soid] = 0
 oconn.close() 

 #Ga0310925_0000001       GeneMark.hmm-2 v1.05    CDS     7       1149    53.99   +       0       ID=Ga0310925_0000001_7_1149;
 with open(gff) as rh:
  for line in rh:
   (soid,source,ftype,start,stop,score,strand,phase,attr) = line.rstrip().split('\t')
   if soid in soids: 
    ngenes[soid] += 1
    if ftype == 'rRNA' : 
     #D=Ga0310925_0000412_1638_3189;e-value=0;model=RF00177;accession=SSU_rRNA_bacteria;model_start=1;model_end=1533;product=16S ribosomal RNAa
     for atr_vals in attr.split(";"):
       key,val = atr_vals.split("=")
       if key == 'product': 
        if val.startswith('16S'):num_16s[soid] += 1
        if val.startswith('5S'):num_5s[soid] += 1
        if val.startswith('23S'):num_23s[soid] += 1
    elif ftype == 'tRNA' : num_tRNA[soid] += 1 

 for escaff in soids:
  sql = "update bin_scaffolds set gene_count = " + str(ngenes[escaff]) + ",scaf_num_16s = " + str(num_16s[escaff]) + ",scaf_num_5s = " + str(num_5s[escaff]) + ",scaf_num_23s = " + str(num_23s[escaff]) + ",scaf_num_tRNA = " + str(num_tRNA[escaff]) + " where scaffold_id = \"" + str(escaff) + "\""
  isqls.append(sql)

 #write metadata to sqlite database
 mbin_sdb.write_sdb(osdb, isqls)

def query_bin_metadata(osdb, taxon_oid):
 #use existing scaffold information to fill bin metadata
 bins = dict()
 isqls = []

 #connect to current sdb
 oconn = sqlite3.connect(osdb)
 oc = oconn.cursor()
 osql = "select bin_name,scaffold_id from bin_scaffolds"
 oc.execute(osql)
 for result in oc:
  bname, soid = result
  if bname not in bins: bins[bname] = set()
  bins[bname].add(soid)

 for ebin in bins:
  total_bases = 0
  gene_count = 0
  bin_num_16s = 0
  bin_num_5s = 0
  bin_num_23s = 0
  bin_num_trna = 0

  slist = "(\"" + "\",\"".join(bins[ebin]) + "\")"
  msql = 'select scaffold_length,gene_count,scaf_num_16s,scaf_num_5s,scaf_num_23s,scaf_num_trna from bin_scaffolds where scaffold_id in ' + slist
  mvals = oc.execute(msql)
  for mval in mvals:
   slen, sgenecount, s16s, s5s, s23s, strna = mval
   total_bases += slen
   gene_count += sgenecount
   bin_num_16s += s16s
   bin_num_5s += s5s
   bin_num_23s += s23s
   bin_num_trna += strna

  sql = "update bin set num_16s = " + str(bin_num_16s) + ",num_5s = " + str(bin_num_5s) + ",num_23s = " + str(bin_num_23s) + ",num_tRNA = " + str(bin_num_trna) + ",total_bases = " + str(total_bases) + ",gene_count=" + str(gene_count) + " where bin_name = \"" + str(ebin) + "\""
  #sql = "update bin set num_16s = " + str(bin_num_16s) + ",num_5s = " + str(bin_num_5s) + ",num_23s = " + str(bin_num_23s) + ",num_tRNA = " + str(bin_num_trna) + ",gene_count=" + str(gene_count) + " where bin_name = \"" + str(ebin) + "\""
  #print sql
  isqls.append(sql)

 oconn.close()
 #write metadata to sqlite database
 mbin_sdb.write_sdb(osdb, isqls)
