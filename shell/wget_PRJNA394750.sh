#!/bin/bash
mkdir data
cd data
for ((i=394, j=4; i<=399; i++,j++));
do 
  mkdir SRR5833${i}; 
  cd SRR5833${i}; 
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR583/00${j}/SRR5833${i}/SRR5833${i}_1.fastq.gz; 
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR583/00${j}/SRR5833${i}/SRR5833${i}_2.fastq.gz; 
  cd ..; 
done
# Originally wrong, as:
# for i in ((i=400, j=0; i<=408; i++,j++)); # don't mix languages!!
for ((i=400, j=0; i<=408; i++,j++));
do 
  mkdir SRR5833${i}; 
  cd SRR5833${i}; 
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR583/00${j}/SRR5833${i}/SRR5833${i}_1.fastq.gz; 
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR583/00${j}/SRR5833${i}/SRR5833${i}_2.fastq.gz; 
  cd ..; 
done
cd .. 