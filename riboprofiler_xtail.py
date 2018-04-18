#!/usr/bin/python
import sys
import subprocess
import shlex
from multiprocessing import Pool
import glob
import os
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--species", help="species human or mouse (default human)", type=str)
parser.add_argument("RPF", help="ordered list of RPF samples: rpf1.fastq.gz,rpf2,fastq.gz", type=str)
parser.add_argument("RNA", help="ordered list of RNA samples: rna1.fastq.gz,rna2,fastq.gz", type=str)
parser.add_argument("conditions", help="input for R: cpntrol,control,sample,sample", type=str)
args = parser.parse_args()

RPFs = args.RPF.split(',')
RNAs = args.RNA.split(',')

conditions = args.conditions

gtffile = "/avicenna/genomes/hg19/hg19_genes.gtf"
genfile = "/avicenna/genome/hg19/hg19.fa"
gendir = "/avicenna/genomes/hg19"
contam = "/avicenna/genomes/hg19/contam/RNAcontam"
if (args.species=="mouse"):
  gtffile = "/avicenna/genomes/mm10/mm10_genes.gtf"  
  genfile = "/avicenna/genome/mm10/mm10.fa"
  gendir = "/avicenna/genomes/mm10"
  contam = "/avicenna/genomes/mm10/contam/RNAcontam"

gindex = gtffile
gindex = re.sub("\.gtf$", "", gindex)

genindex = genfile
genindex = re.sub("\.fa$", "", genindex)


indir = os.getcwd()

logfile=indir+"/pipeline.txt"
log= open(logfile, "wt")
print "Step 1: Linker removal and quality trimming"
for fn in RPFs+RNAs:
  ofn = re.sub("\.fastq\.gz$", ".trim.fastq.gz", fn)
  cmdline='cutadapt -q 15 -m 20 -a AGATCGGAAGAGCACACGTCT -o %s %s' % (ofn,fn)
  log.write("%s\n" % (cmdline))
  print cmdline
  subprocess.call(cmdline,shell=True)

print "Step 2: Removing RNA contaminants (rRNA and tRNA)"
for fn in RPFs+RNAs:
  ifn = re.sub("\.fastq\.gz$", ".trim.fastq.gz", fn)
  ofn = re.sub("\.fastq\.gz$", ".trim.uncontam.fastq.gz", fn)
  cmdline='bowtie2 -p 24 --end-to-end --un-gz=%s -x %s -U %s 2>> %s/stats.txt > %s/aln.out;' % (ofn,contam,ifn,indir,indir)
  log.write("%s\n" % (cmdline))
  print cmdline
  subprocess.call(cmdline,shell=True)
  
print "Step 3: Run Tophat"
for fn in RPFs+RNAs:
  ifn = re.sub("\.fastq\.gz$", ".trim.uncontam.fastq.gz", fn)
  thdir = re.sub("\.fastq\.gz$", "", fn)
  if not(os.path.exists(thdir)):
    os.mkdir(thdir)
  
  cmdline='STAR --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --runThreadN 12 --sjdbGTFfile %s  --genomeDir %s --readFilesIn %s --outFileNamePrefix %s_' % (gtffile,gendir,ifn,thdir)
  #cmdline='tophat2 -p 12 -G %s -o %s --library-type fr-secondstrand --transcriptome-only --b2-sensitive --transcriptome-index=%s %s %s' % (gtffile,thdir,gindex,genindex, ifn)
  log.write("%s\n" % (cmdline))
  print cmdline
  subprocess.call(cmdline,shell=True)
  cmdline='mv %s_Aligned.sortedByCoord.out.bam %s.bam' % (thdir,thdir)
  print cmdline
  subprocess.call(cmdline,shell=True)
  cmdline='samtools index %s.bam' % (thdir)
  print cmdline
  subprocess.call(cmdline,shell=True)

print "Step 4.1: Count reads for RPF samples (exluding those at start codon)"
RPF_count_files=[]
for fn in RPFs:
  #bam = re.sub("\.fastq\.gz$", "_tophat/accepted_hits.bam", fn)
  bam = re.sub("\.fastq\.gz$", ".bam", fn)
  ofn = re.sub("\.fastq\.gz$", ".trim.uncontam.cnt", fn)
  RPF_count_files.append(ofn)
  cmdline='python /nvme/bins/riboprofiler/RPF_counts_CDS.py %s %s > %s' % (bam,gtffile,ofn)
  log.write("%s\n" % (cmdline))
  print cmdline
  subprocess.call(cmdline,shell=True)

print "Step 4.2: Count reads for RNA samples"
RNA_count_files=[]
for fn in RNAs:
  #bam = re.sub("\.fastq\.gz$", "_tophat/accepted_hits.bam", fn)
  bam = re.sub("\.fastq\.gz$", ".bam", fn)
  ofn = re.sub("\.fastq\.gz$", ".trim.uncontam.cnt", fn)
  RNA_count_files.append(ofn)
  cmdline='htseq-count --format=bam --stranded=yes %s %s > %s' % (bam,gtffile,ofn)
  log.write("%s\n" % (cmdline))
  print cmdline
  subprocess.call(cmdline,shell=True)

###Run xtail

cmdline = "Rscript /nvme/bins/riboprofiler/xtail_analysis.R %s %s %s %s" % (",".join(RPF_count_files), ",".join(RNA_count_files), conditions, indir)
print cmdline
log.write("%s\n" % (cmdline))
subprocess.call(cmdline,shell=True)

log.close()

