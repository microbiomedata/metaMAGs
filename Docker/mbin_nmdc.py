#! /usr/bin/env python
import os,sys,re,shutil
import subprocess, sqlite3, argparse,datetime
from Bio import SeqIO
import mbin_sdb, mbin_cleanup, mbin_qassn
import logfile, mbin_metadata
#njvarghese 2020

def preprocess(log,toid,fna,sam,mapfn,numCPU):
 bamfn = toid + '.bam.sorted'
 depfn = toid + '.depth'
 mapdepfn = depfn + '.mapped'
 sam_filename, sam_file_extension = os.path.splitext(sam)
 tmp_dir = os.path.dirname(os.path.abspath(sam))
 #convert to BAM and sort
 if 'bam' in sam_file_extension.lower():
  log.info('Working on BAM input')
  shutil.copy(sam,bamfn)
 else:
  if not os.path.isfile(bamfn):
   log.info('Converting to BAM and sorting using samtools')
   ps = subprocess.Popen(('/usr/bin/time','samtools', 'view', '-@', str(numCPU), '-u', sam), stdout=subprocess.PIPE)
   ps2 = subprocess.Popen(('/usr/bin/time', 'samtools', 'sort', '-@', str(numCPU), '-o', bamfn, '-T', tmp_dir), stdin=ps.stdout,stdout=subprocess.PIPE)
   ps.stdout.close()
   output=ps2.communicate()[0]

 if not os.path.isfile(depfn):
  #use metabat to generate depth file
  log.info('Running metabat to generate depth file')
  #subprocess.check_output(['/usr/bin/time', 'shifter','--image=docker:metabat/metabat:latest','jgi_summarize_bam_contig_depths', '--outputDepth', depfn, bamfn])
  subprocess.check_output(['/usr/bin/time', 'jgi_summarize_bam_contig_depths', '--outputDepth', depfn, bamfn])

 if mapfn:
  #map contig names
  log.info('Mapping contig names')
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
    sid = str(vals.pop(0))
    if sid in ids : wh.write(ids[sid] + "\t" + "\t".join(vals) + '\n')
  wh.close()
  log.info(mapdepfn)
  return(mapdepfn)
 else:
  return(depfn)


def metabat(log,odir,fnafn, covfn, numCPU):
 #run metabat
 basename = odir + '/bins'
 #m_rcode = subprocess.check_call(["/usr/bin/time","shifter","--image=docker:metabat/metabat:latest", "metabat2", "-i", fnafn, "-a", covfn, "-o", basename, "-t", "32", "--unbinned", "-m", "3000","--minS", "80", "--maxEdges","500"])
 m_rcode = subprocess.check_call(["/usr/bin/time","metabat2", "-i", fnafn, "-a", covfn, "-o", basename, "-t",str(numCPU), "--unbinned", "-m", "3000","--minS", "80", "--maxEdges","500"])
 if m_rcode !=0 : sys.exit('Metabat failed. Please check\n')

def move_bins(log,odir, bins_dir):
 if not os.path.isdir(bins_dir) : subprocess.check_call(["mkdir",bins_dir])
 numbins = 0
 for efn in os.listdir(odir):
  fpath = odir+'/'+efn
  if os.path.isfile(fpath):
   if efn == 'bins.unbinned.fa' : unbinned = fpath
   elif re.match(r'bins\.\d+\.fa',efn):
    subprocess.check_call(['mv',fpath, bins_dir])
    numbins += 1

 #Move unbinned contigs of size >=500000 bp to their own bin
 for record in SeqIO.parse(unbinned, "fasta"):
  if (len(record.seq) > 500000) :
  #write to their own bin
   numbins += 1
   new_fn = bins_dir + '/bins.' + str(numbins) + '.fa'
   log.info('New bin:' + newfn)
   wh = open(bins_dir + '/bins\.' + numbins + '.fa')
   wh.write('>' + record.id + '\n' + record.seq + '\n')

def checkm(log,clean_bins_dir,checkm_dir,checkm_qa_out, numCPU, pplacerCPU):
 #run checkm
 cl_rcode = subprocess.check_call(["/usr/bin/time", "checkm", "lineage_wf", "--pplacer_threads", str(pplacerCPU) , "-x", "fa", "-t", str(numCPU), clean_bins_dir, checkm_dir])
 if cl_rcode !=0 : sys.exit('Checkm lineage_wf failed. Please check\n')

 cq_rcode = subprocess.check_call(["/usr/bin/time", "checkm","qa", "-f", checkm_qa_out, checkm_dir + "/lineage.ms", checkm_dir])
 if cq_rcode !=0 : sys.exit('Checkm qa failed. Please check\n')

def gtdbtk_lineage(log,bins_dir, godir, sdb, numCPU, pplacerCPU, scratchdir):
 if not os.path.isdir(scratchdir) : subprocess.check_call(["mkdir",scratchdir])
 rcode = subprocess.check_call(["/usr/bin/time", "gtdbtk", "classify_wf", "--pplacer_cpus", str(pplacerCPU), "--cpus", str(numCPU) , "--genome_dir", bins_dir, "--out_dir", godir, "--extension", "fa"])
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
     g_domain = g_phylum = g_class = g_order = g_family = g_genus = g_species = ""
     if len(lineage_arr)>0 : cmd = "gtdbtk_domain = \"" + lineage_arr[0] + "\""
     if len(lineage_arr)>1 : cmd+= ",gtdbtk_phylum = \"" + lineage_arr[1] + "\""
     if len(lineage_arr)>2: cmd+= ",gtdbtk_class = \"" + lineage_arr[2] + "\""
     if len(lineage_arr)>3: cmd+= ",gtdbtk_order = \"" + lineage_arr[3] + "\""
     if len(lineage_arr)>4: cmd+= ",gtdbtk_family = \"" + lineage_arr[4] + "\""
     if len(lineage_arr)>5: cmd+= ",gtdbtk_genus = \"" + lineage_arr[5] + "\""
     if len(lineage_arr)>6: cmd+= ",gtdbtk_species = \"" + lineage_arr[6] + "\""

     if cmd:
      usql = 'update bin set ' + cmd + ' where bin_name = \'' + str(binid) + '\''
      #sys.stderr.write(usql + '\n')
      c.execute(usql)
 conn.commit()
 conn.close() 

def run(toid, fna, sam, gff, map_fn, domain_fn, numCPU, pplacerCPU, scratchDir):
 #predefined file names
 odir = os.getcwd()
 bins_dir = odir + '/metabat-bins'
 clean_bins_dir = odir + '/clean-metabat-bins'
 hqmq_bins_dir = odir + '/hqmq-metabat-bins'
 checkm_dir= odir + '/checkm-out'
 checkm_qa_out = odir + '/checkm_qa.out'
 godir = odir + '/' + 'gtdbtk_output'
 today = str(datetime.date.today())
 sdb_name = odir + '/mbin-' + today + '.sqlite'
 no_load_fn = odir + '/noload.mbin'
 mbin_complete_fn = odir + '/complete.mbin'
 log = logfile.startlog('mbin-nmdc')

 #get fna files and coverage files and run metabat
 log.info('..preprocess module')
 dep_fn = preprocess(log, toid, fna, sam, map_fn, numCPU)

 log.info('..run metabat')
 
 #check if any bins were generated by metabat
 if not os.path.isdir(bins_dir) : subprocess.check_call(["mkdir",bins_dir])
 num_bins = len([name for name in os.listdir(bins_dir) if os.path.isfile(os.path.join(bins_dir, name))])
 if num_bins >= 1 : pass
 else:
  metabat(log, odir, fna, dep_fn,numCPU)

  #move bin files to directory
  move_bins(log, odir, bins_dir)

 #check if any bins were generated by metabat
 num_bins = len([name for name in os.listdir(bins_dir) if os.path.isfile(os.path.join(bins_dir, name))])
 if num_bins >= 1 : pass
 else:
  #touch file to avoid pickup
  subprocess.check_call(['touch',no_load_fn])
  subprocess.check_call(['touch',mbin_complete_fn])
  log.info('No bins generated by metabat. Created flag file to avoid pickup:'+ no_load_fn)
  log.info('Completion flag file:' + mbin_complete_fn)
  sys.exit(0)

 #create sqlitedb
 log.info('..creating sdb')
 mbin_sdb.createsdb(sdb_name)

 #write bin informatin from bins_dir to sdb 
 log.info('..write bins to sdb')
 mbin_sdb.writebins(sdb_name,bins_dir)

 #get scaff info
 log.info('..get per scaffold metadata')
 mbin_metadata.parse_annotation(sdb_name,gff)

 #query scaffold lineage
 if domain_fn:
  log.info('..query scaffold lineage and write to sdb')
  mbin_metadata.query_slineage(sdb_name, domain_fn)

  #clean up bins remove any scaffolds that match another phylum
  log.info('..clean bins based on domain information')
  mbin_cleanup.run(sdb_name, bins_dir, clean_bins_dir)
 else:
  clean_bins_dir = bins_dir

 #get rRNA and tRNA count per bin
 log.info('..get per bin metadata')
 mbin_metadata.query_bin_metadata(sdb_name, toid)

 #run checkm
 log.info('..run checkm')
 if not os.path.isfile(checkm_qa_out):
  checkm(log,clean_bins_dir,checkm_dir,checkm_qa_out,numCPU, pplacerCPU)

 #mimag based hq,mq,lq assignment & move hq.mq to sep dir
 log.info('..assign quality and move hq,mq to sep dir')
 mbin_qassn.run(checkm_qa_out, clean_bins_dir, hqmq_bins_dir, sdb_name)

 #check if there are any hq/mq bins
 num_hqmq = len([name for name in os.listdir(hqmq_bins_dir) if os.path.isfile(os.path.join(hqmq_bins_dir, name))])
 if num_hqmq >=1 :
 
  #gtdb-tk lineage
  log.info('..run gtdb-tk on hq,mq bins')
  gtdbtk_lineage(log,hqmq_bins_dir, godir, sdb_name, numCPU, pplacerCPU, scratchDir)
 
 else:
  #touch file to avoid pickup
  subprocess.check_call(['touch',no_load_fn])  
  log.info('No HQ/MQ bins. Created flag file to avoid pickup:'+ no_load_fn)

 #write flag file
 subprocess.check_call(['touch',mbin_complete_fn])
 log.info('Completion flag file:' + mbin_complete_fn)

if __name__ == '__main__':
 parser = argparse.ArgumentParser(description='Metabin binning pipeline for IMG modified for NMDC')
 parser.add_argument("toid", help = "Taxon ID or name of metagenome")
 parser.add_argument("fna", help = "Metagenome fna")
 parser.add_argument("sam", help = "Mapped reads in SAM format or sorted BAM format")
 parser.add_argument("gff", help = "GFF file from IMG annotation pipeline")
 parser.add_argument("--map", help = "MAP file containing mapping of headers between SAM and FNA : ID in FNA<tab>ID in GFF")
 parser.add_argument("--cpu", default=8, type=int, help = "number of CPU [8]")
 parser.add_argument("--pplacer_cpu", default=1, type=int, help = "number of pplacer CPU [1]")
 parser.add_argument("--scratch_dir", type=str, help = "scratch directory of pplacer. use disk instead of memory")
 parser.add_argument("--domain", help = "Per scaffold domain information : soid<tab>domain")
 
 args = parser.parse_args()
 run(args.toid, args.fna, args.sam, args.gff, args.map, args.domain, args.cpu, args.pplacer_cpu, args.scratch_dir)
