#!/bin/bash
for fn in VC_out/GSM27058{80..94};
do
samp=`basename ${fn}`
echo "Extracting sample ${samp}"
gunzip -vf ${fn}.g.vcf.gz 
done