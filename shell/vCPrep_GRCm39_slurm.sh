#!/bin/bash

#SBATCH --job-name vCPrep_GRCm39 ## The name that will show up in the queue
#SBATCH --time 1-12:00:00	    ## Time for analysis (day-hour:min:sec)
## SBATCH --output    ## Filename of the output; default is slurm-[joblD].out
#SBATCH --ntasks=6              ## Number of tasks (analyses) to run; default = 1
#SBATCH --cpus-per-task  	    ## The num of threads the code will use; default = 1
#SBATCH --mem-per-cpu           ## Memory per allocated CPU

srun /path/to/gatk-version/gatk IndexFeatureFile -I mus_musculus.vcf.gz
srun gunzip -vf mus_musculus.vcf.gz
srun /path/to/gatk-version/gatk SplitNCigarReads -R Mus_musculus.GRCm39.dna_sm.primary_assembly.fa \
    --sequence-dictionary Mus_musculus.GRCm39.dna_sm.primary_assembly.dict \
    -I alignment_out/marked_duplicates.bam \
    -O alignment_out/SplitNCigarReads.bam \
    --tmp-dir $1
srun /path/to/gatk-version/gatk BaseRecalibrator -R Mus_musculus.GRCm39.dna_sm.primary_assembly.fa \
    --known-sites mus_musculus_incl_consequences.vcf \
    -I alignment_out/SplitNCigarReads.bam \
    -O alignment_out/baseRecal.table \
    --tmp-dir $1
srun /path/to/gatk-version/gatk ApplyBQSR -R Mus_musculus.GRCm39.dna_sm.primary_assembly.fa \
    -bqsr alignment_out/baseRecal.table \
    -I alignment_out/SplitNCigarReads.bam \
    -O alignment_out/baseRecal.bam \
    --tmp-dir $1
srun /path/to/gatk-version/gatk AnalyzeCovariates -bqsr alignment_out/baseRecal.table \
    -csv alignment_out/AnalyseCovariates.csv \
    --tmp-dir $1
    
    
