#!/usr/bin/env python3
import os, json, time, sys, glob
from datetime import datetime
import collections
import subprocess, sqlite3

def line_count(fname):
    count = len(open(fname).readlines(  ))
    return count

def count_fasta(fname):
  n=0
  with open(fname, "rb") as f:
    for line in f:
      if line.startswith(b'>'):
        n += 1
  return n

def mag_meta(dbname, output):
 conn = sqlite3.connect(dbname)
 c = conn.cursor()
 #qsql="select bin_name, count(scaffold_id) from bin_scaffolds group by bin_name;"
 qsql="select bin_name, scaffold_id from bin_scaffolds;"
 rvals = c.execute(qsql)
 d = collections.defaultdict()
 d_list = collections.defaultdict(list)
 total_bin_contig = 0
 for row in rvals:
  d_list[row[0]].append(row[1])  
  d[row[0]] = d[row[0]] + 1 if row[0] in d else 1
  total_bin_contig += 1

 qsql="select * from bin order by case when bin_quality = 'HQ' then 1 when bin_quality = 'MQ' then 2 when bin_quality = 'LQ' then 3 end,completeness DESC;"
 rvals = c.execute(qsql)
 out_list=[]
 f = open(output,'w')
 f.write('\t'.join(['Bin ID','Bin Quality','GTDB-TK lineage','Bin Completeness','Bin Contamintation','Total Number of Bases','Number of genes','Num of 5s rRNA','Num of 16s rRNA','Num of 23s rRNA','Num of tRNA','scaffold IDs of members'])+ "\n")
 for row in rvals:
  tmp_d = collections.OrderedDict()
  tmp_d["bin_name"] = row[0]
  tmp_d["number_of_contig"] = d[row[0]]
  tmp_d["completeness"] = row[8]
  tmp_d["contamination"] = row[9]
  tmp_d["total_bases"] = row[10]
  tmp_d["gene_count"] = row[11]
  tmp_d["bin_quality"] = row[12]
  tmp_d["num_16s"] = row[13]
  tmp_d["num_5s"] = row[14]
  tmp_d["num_23s"] = row[15]
  tmp_d["num_tRNA"] = row[16]
  tmp_d["gtdbtk_domain"] = row[17]
  tmp_d["gtdbtk_phylum"] = row[18]
  tmp_d["gtdbtk_class"] = row[19]
  tmp_d["gtdbtk_order"] = row[20]
  tmp_d["gtdbtk_family"] = row[21]
  tmp_d["gtdbtk_genus"] = row[22]
  tmp_d["gtdbtk_species"] = row[23]
  tmp_d["members_id"] = d_list[row[0]] 
  out_list.append(tmp_d)
  f.write('\t'.join([row[0],str(row[12]),';'.join(str(item) for item in row[17:24]),str(row[8]),str(row[9]), str(row[10]), str(row[11]), str(row[14]),str(row[13]),str(row[15]),str(row[16]),','.join(tmp_d["members_id"])] ) + "\n")
  
 conn.close()
 return(out_list, total_bin_contig)


if __name__ == "__main__":
       
	outdir= os.getcwd() if len(sys.argv) < 2 else sys.argv[1]
	outfile= outdir + "/MAGs_stats.json"
	outtsvfile= outdir + "/MAGs_stats.tsv"

	#input contigs #
	#tooShort (<3000) contigs #
	#lowDepth contigs #
	#unbinned contigs #
	#binned contigs #
	#bins list:

	metadata=collections.OrderedDict()
	mag_list=[]
	total_bin_contig_num = 0

	#script_file = outdir + "/script"
	#script_time = datetime.fromtimestamp(os.path.getmtime(script_file)).strftime("%Y-%m-%d")
	#end_file = outdir + "/complete.mbin"
	#end_time = datetime.fromtimestamp(os.path.getmtime(end_file)).strftime("%Y-%m-%d")
	depth_file= glob.glob(outdir + "/*.depth") 
	input_contig_num = line_count(depth_file[0]) - 1
	too_short_file = outdir + "/bins.tooShort.fa"
	too_short_contig_num = 0 if not os.path.isfile(too_short_file) else count_fasta(too_short_file)
	lowDepth_file = outdir + "/bins.lowDepth.fa"
	lowDepth_contig_num = 0 if not os.path.isfile(lowDepth_file) else count_fasta(lowDepth_file)
	unbinned_file = outdir + "/bins.unbinned.fa"
	unbinned_contig_num = 0 if not os.path.isfile(unbinned_file) else count_fasta(unbinned_file)
	sql_files = glob.glob(outdir + "/*.sdb")
	if len(sql_files) > 0:
		mag_list, total_bin_contig_num = mag_meta(sql_files[0],outtsvfile)
		sqlfile_time=datetime.fromtimestamp(os.path.getmtime(sql_files[0])).strftime("%Y-%m-%d")

	metadata['input_contig_num'] = input_contig_num
	metadata['too_short_contig_num'] = too_short_contig_num
	metadata['lowDepth_contig_num'] = lowDepth_contig_num
	metadata['unbinned_contig_num'] = unbinned_contig_num
	metadata['binned_contig_num'] = total_bin_contig_num
	metadata['mags_list'] = mag_list

	with open(outfile, 'w') as outfile:
		json.dump(metadata, outfile,indent=4)
