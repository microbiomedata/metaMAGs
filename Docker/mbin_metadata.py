#!/usr/bin/env python3
#njvarghese
import sys,sqlite3
import mbin_sdb


#if lineage is provided as lineage SDB (post loading IMG process)
def query_slineage(osdb, lineage_sdb, log):
 soids = set()
 isqls = []

 log.info('query sdb to get bin to scaffold connections.')
 #connect to current sdb
 oconn = sqlite3.connect(osdb)
 oc = oconn.cursor()
 osql = "select bin_name,scaffold_id from bin_scaffolds"
 oc.execute(osql)
 for result in oc: 
  bname, soid = result
  soids.add(soid)
 oconn.close()

 #connect to IMG lineage sdb
 log.info('query IMG lineage sdb to get per scaffold lineage.')
 lconn = sqlite3.connect(lineage_sdb)
 lc = lconn.cursor()
 for escaff in soids:
  query = "select lineage from contig_lin where scaffold_oid = \"" + escaff + "\""
  lc.execute(query)
  for result in lc:
   lineage_arr = result[0].split(";")
   cmd = None
   if len(lineage_arr)>0 : cmd = "scaf_domain = \"" + lineage_arr[0] + "\""
   if len(lineage_arr)>1 : cmd+= ",scaf_phylum = \"" + lineage_arr[1] + "\""
   if len(lineage_arr)>2: cmd+= ",scaf_class = \"" + lineage_arr[2] + "\""
   if len(lineage_arr)>3: cmd+= ",scaf_order = \"" + lineage_arr[3] + "\""
   if len(lineage_arr)>4: cmd+= ",scaf_family = \"" + lineage_arr[4] + "\""
   if len(lineage_arr)>5: cmd+= ",scaf_genus = \"" + lineage_arr[5] + "\""
   if len(lineage_arr)>6: cmd+= ",scaf_species = \"" + lineage_arr[6] + "\""
  
   if cmd: 
    sql = "update bin_scaffolds set " + cmd + " where scaffold_id = \"" + str(escaff) + "\""  
    isqls.append(sql)

 lconn.close()

 log.info('write per scaffold lineage to sdb.')
 #write metadata to sqlite database
 mbin_sdb.write_sdb(osdb, isqls)

def parse_annotation(osdb, gff, log):
 soids = set()
 sql_params = []
 num_16s = dict()
 num_5s = dict()
 num_23s = dict()
 num_tRNA = dict()

 #connect to current sdb
 log.info('query sdb to get required scaffold list.')
 oconn = sqlite3.connect(osdb)
 oc = oconn.cursor()
 osql = "select bin_name,scaffold_id from bin_scaffolds"
 oc.execute(osql)
 for result in oc:
  bname, soid = result
  soids.add(soid)
  num_16s[soid] = 0
  num_5s[soid] = 0
  num_23s[soid] = 0
  num_tRNA[soid] = 0
 oconn.close() 

 log.info('read gff to get num of rRNA and tRNA.')
 #Ga0310925_0000001       GeneMark.hmm-2 v1.05    CDS     7       1149    53.99   +       0       ID=Ga0310925_0000001_7_1149;
 with open(gff) as rh:
  for line in rh:
   (soid,source,ftype,start,stop,score,strand,phase,attr) = line.rstrip().split('\t')
   if soid in soids: 
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
  sql_params.append((str(num_16s[escaff]),str(num_5s[escaff]),str(num_23s[escaff]),str(num_tRNA[escaff]),str(escaff)))

 log.info('update sdb with per scaffold counts.')
 sql =  "update bin_scaffolds set scaf_num_16s = ? ,scaf_num_5s = ? ,scaf_num_23s = ? ,scaf_num_tRNA = ? where scaffold_id = ?"
 #write metadata to sqlite database
 mbin_sdb.write_sdb_executemany(osdb, sql,sql_params)

def query_bin_metadata(osdb):
 #use existing scaffold information to fill bin metadata
 bins = dict()
 sql_params = []

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
  bin_num_16s = 0
  bin_num_5s = 0
  bin_num_23s = 0
  bin_num_trna = 0

  slist = "(\"" + "\",\"".join(bins[ebin]) + "\")"
  msql = 'select scaffold_length,scaf_num_16s,scaf_num_5s,scaf_num_23s,scaf_num_trna from bin_scaffolds where scaffold_id in ' + slist
  mvals = oc.execute(msql)
  for mval in mvals:
   slen, s16s, s5s, s23s, strna = mval
   if slen: total_bases += slen
   bin_num_16s += s16s
   bin_num_5s += s5s
   bin_num_23s += s23s
   bin_num_trna += strna

  sql_params.append((str(bin_num_16s),str(bin_num_5s),str(bin_num_23s),str(bin_num_trna) ,str(total_bases),str(ebin)))

 oconn.close()
 #write metadata to sqlite database
 sql = "update bin set num_16s =? ,num_5s =? ,num_23s =? ,num_tRNA =? ,total_bases =? where bin_name =?"

 mbin_sdb.write_sdb_executemany(osdb, sql, sql_params)

#if lineage is provided as TSV : pre-loading IMG process
def query_slineage_tsv(osdb, lineage_tsv, log):
 soids = set()
 sql_params = []

 log.info('query sdb to get bin to scaffold connections.')
 #connect to current sdb
 oconn = sqlite3.connect(osdb)
 oc = oconn.cursor()
 osql = "select bin_name,scaffold_id from bin_scaffolds"
 oc.execute(osql)
 for result in oc: 
  bname, soid = result
  soids.add(str(soid))
 oconn.close()

 #connect to IMG lineage sdb
 log.info('parse IMG lineage tsv to get per scaffold lineage.')
 
 with open(lineage_tsv) as rh:
  for line in rh:
   scaf_domain=""
   scaf_phylum=""
   scaf_class=""
   scaf_order=""
   scaf_family=""
   scaf_genus=""
   scaf_species=""
   linevals = line.rstrip().split('\t')
   escaff = str(linevals[0])
   if escaff not in soids : continue
   lineage_arr = linevals[1].split(";")

   if len(lineage_arr)>0 : scaf_domain=lineage_arr[0] 
   if len(lineage_arr)>1 : scaf_phylum=lineage_arr[1] 
   if len(lineage_arr)>2 : scaf_class=lineage_arr[2]
   if len(lineage_arr)>3 : scaf_order=lineage_arr[3]
   if len(lineage_arr)>4 : scaf_family=lineage_arr[4]
   if len(lineage_arr)>5 : scaf_genus=lineage_arr[5]
   if len(lineage_arr)>6 : scaf_species=lineage_arr[6]
   
   sql_params.append((scaf_domain,scaf_phylum,scaf_class,scaf_order,scaf_family,scaf_genus,scaf_species,str(escaff)))


 log.info('write per scaffold lineage to sdb.')
 #write metadata to sqlite database
 sql="update bin_scaffolds set scaf_domain=?, scaf_phylum=?, scaf_class=?, scaf_order=?, scaf_family=?, scaf_genus=?, scaf_species=? where scaffold_id=?"
 mbin_sdb.write_sdb_executemany(osdb, sql, sql_params)
