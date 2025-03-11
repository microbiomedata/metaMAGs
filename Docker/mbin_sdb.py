#! /usr/bin/env python3
#njvarghese
import os, re, sqlite3
from Bio import SeqIO

def write_sdb(osdb, isqls):
 qcount = 0
 oconn = sqlite3.connect(osdb)
 oi = oconn.cursor()
 for esql in isqls:
  oi.execute(esql)
  qcount += 1
  if qcount == 500 :
   oconn.commit()
   qcount = 0
 oconn.commit()
 oconn.close()

def write_sdb_executemany(osdb, sql, params , chunk=10000):
 oconn = sqlite3.connect(osdb)
 oi = oconn.cursor()
 for i in range(0, len(params), chunk):
  oi.executemany(sql,params[i:i + chunk]) 
  oi.execute('COMMIT')
 oconn.close()


def createsdb(dbname,log):
 if os.path.isfile(dbname) : 
  log.info('Sqlite database exists.Exiting module')
  return()

 conn = sqlite3.connect(dbname)
 c = conn.cursor()

 csql ='CREATE TABLE bin('
 csql+='bin_name varchar(255), bin_domain varchar(255), bin_phylum varchar(255), bin_class varchar(255), bin_order varchar(255), bin_family varchar(255), bin_genus varchar(255), bin_species varchar(255), completeness float, contamination float, total_bases integer, gene_count integer, bin_quality varchar(255), num_16s integer, num_5s integer, num_23s integer, num_tRNA integer, gtdbtk_domain varchar(255), gtdbtk_phylum varchar(255),gtdbtk_class varchar(255), gtdbtk_order varchar(255), gtdbtk_family varchar(255), gtdbtk_genus varchar(255), gtdbtk_species varchar(255))'
 log.debug(csql)
 c.execute(csql)
 conn.commit()

 nsql ='CREATE TABLE bin_scaffolds('
 nsql+='bin_name varchar(255), scaffold_id varchar(255), scaffold_length integer, gc float, gene_count integer, scaf_domain varchar(255), scaf_phylum varchar(255), scaf_class varchar(255), scaf_order varchar(255), scaf_family varchar(255), scaf_genus varchar(255), scaf_species varchar(255), scaf_num_16s integer, scaf_num_5s integer, scaf_num_23s integer, scaf_num_trna integer)'
 log.debug(nsql)
 c.execute(nsql)
 conn.commit()

 conn.close()

def createindex(dbname):
 conn = sqlite3.connect(dbname)
 c = conn.cursor()
 csql = "CREATE UNIQUE INDEX idx_scaffold ON bin_scaffolds (scaffold_id);"
 c.execute(csql)
 csql = "CREATE UNIQUE INDEX idx_bin ON bin (bin_name);"
 c.execute(csql)
 conn.commit()
 conn.close()

def query(dbname, qsql):
 conn = sqlite3.connect(dbname)
 c = conn.cursor()
 rvals = c.execute(qsql)
 conn.close()
 return(rvals)

def update(dbname, usql):
 conn = sqlite3.connect(dbname)
 c = conn.cursor()
 c.execute(usql)
 conn.commit()
 conn.close()

def writebins(dbname, bins_dir):
 binname_list_tuple=[]
 biname_id_list_tuple=[]
 sql = 'insert into bin(bin_name) values (?)'
 sql2 = 'insert into bin_scaffolds(bin_name,scaffold_id) values (?,?)'
 for efile in os.listdir(bins_dir):
  epath = bins_dir + '/' + efile
  if os.path.isfile(epath):
   binname = re.sub(r'\.fa*','',efile)
   binname_list_tuple.append((str(binname),))
   for record in SeqIO.parse(epath, 'fasta'):
    biname_id_list_tuple.append((str(binname),str(record.id)))
    #log.debug(rsql)
   
 write_sdb_executemany(dbname,sql,binname_list_tuple)
 write_sdb_executemany(dbname,sql2,biname_id_list_tuple)