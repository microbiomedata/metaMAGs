#! /usr/bin/env python
import os,sys,re
import subprocess, sqlite3
from Bio import SeqIO
#njvarghese 2020

def write_sdb(osdb, isqls):
 #write metadata to sqlite database
 sys.stderr.write('writing to database\n')
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

def createsdb(dbname):
 if os.path.isfile(dbname) : return('Sqlite database exits.Exiting module')

 conn = sqlite3.connect(dbname)
 c = conn.cursor()

 csql ='CREATE TABLE bin('
 csql+='bin_name varchar(255), bin_domain varchar(255), bin_phylum varchar(255), bin_class varchar(255), bin_order varchar(255), bin_family varchar(255), bin_genus varchar(255), bin_species varchar(255), completeness float, contamination float, total_bases integer, gene_count integer, bin_quality varchar(255), num_16s integer, num_5s integer, num_23s integer, num_tRNA integer, gtdbtk_domain varchar(255), gtdbtk_phylum varchar(255),gtdbtk_class varchar(255), gtdbtk_order varchar(255), gtdbtk_family varchar(255), gtdbtk_genus varchar(255), gtdbtk_species varchar(255))'

 c.execute(csql)
 conn.commit()

 nsql ='CREATE TABLE bin_scaffolds('
 nsql+='bin_name varchar(255), scaffold_id varchar(255), scaffold_length integer, gc float, gene_count integer, scaf_domain varchar(255), scaf_phylum varchar(255), scaf_class varchar(255), scaf_order varchar(255), scaf_family varchar(255), scaf_genus varchar(255), scaf_species varchar(255), scaf_num_16s integer, scaf_num_5s integer, scaf_num_23s integer, scaf_num_trna integer)'

 c.execute(nsql)
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
 conn = sqlite3.connect(dbname)
 c = conn.cursor()

 for efile in os.listdir(bins_dir):
  epath = bins_dir + '/' + efile
  if os.path.isfile(epath):
   binname = re.sub(r'\.fa*','',efile)
   isql = 'insert into bin(bin_name) values(\'' + str(binname) + '\')'
   #print isql
   c.execute(isql)

   for record in SeqIO.parse(epath, 'fasta'):
    rsql = 'insert into bin_scaffolds(bin_name,scaffold_id) values (\'' + str(binname) + '\',\'' + str(record.id) + '\')'
    #print rsql
    c.execute(rsql)

 conn.commit()
 conn.close()      
