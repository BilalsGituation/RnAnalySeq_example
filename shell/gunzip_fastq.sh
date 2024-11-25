#!/bin/bash
for fn in data/SRR5833{394..408};
do
samp=`basename ${fn}`
echo "Extracting sample ${samp}"
gunzip -vf ${fn}/${samp}_1.fastq.gz ${fn}/${samp}_2.fastq.gz
done