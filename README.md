# RNA-Seq NGS Analysis SoftwaRe: RnAnalySeq (working title)

DESCRIPTION:
Example use of work-in-progress software that analyses samples from mouse RNA-Seq experiments.
The reason for sharing at this early stage is to show my independent self-development in the 
bioinformatics field.

This project is written in R, which was mainly chosen due to its support of exploratory data analysis 
via the IDE RStudio.

In its current iteration, the software is able to preprocess mouse RNA-Seq samples, align them against 
the GRCm39 reference genome, then perform variant calling (and its specific preprocessing), differential 
expression analysis, and gene enrichment analysis. This list of capabilities will be expanded in 
subsequent versions.

Any feedback, especially constructive, is welcome and appreciated!

PLEASE NOTE:
Since the main workflow contains calls to GATK (McKenna et al., 2010) and Rbowtie2 
(Schubert, Lindgreen, and Orlando, 2016), it is highly recommended to run this software in a 
POSIX-compatible system such as Linux and macOS.

PLANNED UPDATES:
Next patch(es) should include:
- Update anything related to ensembl release 112 => 113 and test them
- Update GATK version and validate calls to it
- Remove loading of superfluous packages
- Optimise parallelisation of longer workflow steps
- Operations on VCF files
- Development of enrichment analysis capabilities
- Alternative splicing analysis
- PPI analysis
- Fusion analysis

Planned expansions (minor or major versions):
- Accommodate more genomes
- Accommodate more sequencing strategies and platforms
- Restructure everything into functions and (after expanding the user's options) a version-controlled package

REFERENCES (where R's citation() function will not provide them):
- Ensembl (for reference genome, known variants and annotation file) 
- Credit to Rodolfo Briones (Berlin) for valuable planning and feedback contributions 
  - https://github.com/fitosky1/fitosky1, https://www.linkedin.com/in/rodolfo-briones-phd/
- Youtube@Bioinformagician (BioMart video)
- Love MI, Anders S, Kim V, Huber W. (2015). RNA-Seq workflow: gene-level exploratory analysis and differential expression. 
  F1000Res. Oct 14;4:1070. doi: 10.12688/f1000research.7035.1. PMID: 26674615; PMCID: PMC4670015.
	(Updated: https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html)
- McKenna et al. (2010). Original description of the GATK framework
- Poplin et al. (2017). Detailed description of HaplotypeCaller; best reference for germline joint calling
- Wong CK et al., "The p300 and CBP Transcriptional Coactivators Are Required for β-Cell and α-Cell Proliferation.", Diabetes, 2018 Mar;67(3):412-422
- (adapter sequences for Illumina can be found at) https://support-docs.illumina.com/SHARE/AdapterSequences/Content/SHARE/FrontPages/AdapterSeq.htm
- Ensembl (for reference genome, known variants and annotation file) 



