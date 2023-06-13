#! /usr/bin/env python3
import os, sys, subprocess
from Bio import SeqIO

def internal(fnafn, mapfn, bam, log, numcores, combined, sortfna, sortaln):
 
 #predefined variables
 depfn = 'metabat.depth'
 mapdepfn = str(depfn) + '.mapped'
 basename = 'bins' 
 sorted_fna = fnafn + '.sorted'
 sorted_bam = os.path.basename(bam) + '.sorted'
 sampath = None

 #Step1 : Generate DEPTH FILE
 if combined == "yes":
  log.info("Combined flag set to yes. Will use --aln/-a value as directory.")

  #check if files end in .gz, then modify sampath accordingly 
  #added Nov 2022 due to multiple format changes by asm group
  for each_file in os.listdir(bam):
    if each_file.endswith(".bam"): #preferred option
      sampath = bam + '/*.bam'
      break
    elif each_file.endswith(".bam.gz"): 
      sampath = bam + '/*.bam.gz'
      break
    elif each_file.endswith(".sam"): sampath = bam + '/*.sam'
    elif each_file.endswith(".sam.gz"): sampath = bam + '/*.sam.gz'
  
  log.info("path with extension to be used:" + sampath)

  #use metabat to generate depth file from multiple sam file
  log.info('generating depth file using metabat.')
  depth_cmd = 'jgi_summarize_bam_contig_depths --outputDepth ' + depfn + ' ' + sampath
  log.debug(depth_cmd)
  ocode = os.system(depth_cmd)
  if ocode>0 : sys.exit('Fatal: Metabat error')

 elif combined == "no":

  # option to sort alignment file since metabat requires sorted aln
  if sortaln == "yes":
   log.info('sorting alignment file using samtools.')
   subprocess.check_output(['samtools', 'sort', '--threads', str(numcores), '-o', sorted_bam, bam])
  else:
   sorted_bam = bam

  #use metabat to generate depth file
  log.info('generating depth file using metabat.')
  subprocess.check_output(['jgi_summarize_bam_contig_depths', '--outputDepth', depfn, sorted_bam])

 #Step 2: map contig names in depth file if map file is provided
 if mapfn and os.path.isfile(mapfn):
  #map contig names between internal and IMG identifiers
  log.info('renaming contig/scaffold headers in depth file')
  ids = dict()
  with open(mapfn) as rh:
   for line in rh:
    sid,cid = line.rstrip().split('\t')
    ids[sid] = cid

  wh = open(mapdepfn,'w')
  with open(depfn) as fh:
   for line in fh:
    if line.startswith('contigName'):
     wh.write(line.rstrip()+'\n')
     continue
    vals = line.strip().split('\t')
    sid = vals.pop(0)
    if sid in ids: wh.write(ids[sid] + "\t" + "\t".join(vals) + '\n')
  wh.close()
  log.info("depth file with img headers generated: " + str(mapdepfn))
 
 else:
  mapdepfn = depfn

 #STEP 3 : generate bins using metabat
 if combined == "yes":

  #run metabat
  m_rcode = subprocess.check_call(["metabat2", "-i", fnafn, "-a", mapdepfn, "-o", basename, "-t", str(numcores), "--unbinned", "-m", "3000","--minS", "80", "--maxEdges","500", "--seed", "1000"])
  if m_rcode !=0 : sys.exit('Metabat failed. Please check\n')

 elif combined == "no":

  #option to sort fasta file
  #added May 2022 due to issues with AWS annotation
  if sortfna == "yes": 
   log.info("sorting fasta file.")
   seqs = dict()
   for record in SeqIO.parse(fnafn, "fasta"):
    seqs[record.id]=record.seq
   wh = open(sorted_fna,"w")
   for header in sorted (seqs.keys()):
    wh.write(">" + header + '\n' + str(seqs[header]) + "\n")
   wh.close()
  else:
   sorted_fna = fnafn
 
  #run metabat
  m_rcode = subprocess.check_call(["metabat2", "-i", sorted_fna, "-a", mapdepfn, "-o", basename, "-t", str(numcores), "--unbinned", "-m", "3000","--minS", "80", "--maxEdges","500", "--seed", "1000"])
  if m_rcode !=0 : sys.exit('Metabat failed. Please check\n')

