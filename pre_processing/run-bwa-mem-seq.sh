#!/bin/bash

fastq1=$1
fastq2=$2
refseq_path=$3
outdir=$4
prefix=$5
nThreads=$6

if [[ $outdir != '.' ]]
then
  cd $outdir
fi

bwa index $refseq_path

# unzip fastq files
if [[ $fastq1 =~ \.gz$ ]]
then
  cp $fastq1 fastq1.gz
  gunzip fastq1.gz
else
  cp $fastq1 fastq1
fi
  fastq1=fastq1

if [[ $fastq2 =~ \.gz$ ]]
then
  cp $fastq2 fastq2.gz
  gunzip fastq2.gz
else
  cp $fastq2 fastq2
fi
  fastq2=fastq2


#-S            skip mate rescue                                                                               
#       -P            skip pairing; mate rescue performed unless -S also in use
#       -5         for split alignment, take the alignment with the smallest coordinate as primary
#       -M        mark shorter split hits as secondary
# run bwa
bwa mem -t $nThreads -SP5M $refseq_path $fastq1 $fastq2 | samtools view -Shb - > $prefix.bam