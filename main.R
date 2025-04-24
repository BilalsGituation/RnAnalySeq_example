# Run in POSIX (this is only tested in Ubuntu)
options(timeout = 1000)
pacman::p_load(#rfiglet, # install below
               Rbowtie2, 
               stringr,
               tidyverse,
               #QuasR part
               QuasR,
               BSgenome,
               Rsamtools,
               rtracklayer,
               GenomicFeatures,
               txdbmaker,
               Gviz,
               BSgenomeForge,
               #BSgenome.Mmusculus.UCSC.mm10.masked,
               Rhisat2,
               readr,
               #Rsubread,
               #DESeq2 part
               biomaRt,
               apeglm,
               DESeq2,
               dplyr,
               genefilter,
               ggplot2,
               ggbeeswarm,
               glmpca,
               IHW,
               magrittr,
               naniar,
               org.Mm.eg.db,
               pheatmap,
               RColorBrewer, 
               PoiClaClu,
               #tximeta,
               vsn,
               #Enrichment analyses part
               clusterProfiler,
               DOSE,
               enrichplot,
               pathview,
               ReactomePA,
               reactome.db,
               #VarCalling
               sys,#     talks to system
               gsalib,
               vcfR,
               #AS
               ASpli,
               #PPI
               rbioapi#,
               #Fusion analysis
               
)
#Bioconductor Package Maintainer <maintainer@bioconductor.org>

citation("rfiglet") # install below
citation("Rbowtie2") 
citation("stringr")
citation("tidyverse")
citation("QuasR")
citation("BSgenome")
citation("Rsamtools")
citation("rtracklayer")
citation("GenomicFeatures")
citation("txdbmaker")
citation("Gviz")
citation("BSgenomeForge")
citation("Rhisat2")
citation("readr")
citation("Rsubread")
citation("biomaRt")
citation("apeglm")
citation("DESeq2")
citation("dplyr")
citation("genefilter")
citation("ggplot2")
citation("ggbeeswarm")
citation("glmpca")
citation("IHW")
citation("magrittr")
citation("naniar")
citation("org.Mm.eg.db")
citation("pheatmap")
citation("RColorBrewer")
citation("PoiClaClu")
citation("vsn")
citation("clusterProfiler")
citation("DOSE")
citation("enrichplot")
citation("pathview")
citation("ReactomePA")
citation("reactome.db")
citation("sys")
citation("vcfR")
citation("ASpli")
citation("rbioapi")


# STARTING MATERIALS (And their wget scripts where possible):
# - gunzipped fastq files

#run in your terminal "./shell/wget_PRJNA394750.sh"
#               then  "./shell/gunzip_fastq.sh", both without quotes
# or in RStudio
exec_wait("./shell/wget_PRJNA394750.sh")
exec_wait("./shell/gunzip_fastq.sh")

# - metadata table from the sra project (may need implementation of unix line endings) 
#     https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP112610&o=acc_s%3Aa

# - flush, delete or change the "alignment_out" directory unless you have saved preprocessed files

# - mm39 ensembl reference genome (soft masked) available here https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/
# Warning: at end of April 2025, it looks like the next releases of auxiliary ensembl files will come out imminently. 
#         So try accessing https://ftp.ensembl.org/pub/release-115 (which is "Forbidden" at time of writing.)
#         even if it does contain the same reference genome as recent previous releases
exec_wait("wget", "https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz")
exec_wait("gunzip", c("-vf", "Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz"))

# - mm39 ensembl annotation file available here https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/
exec_wait("wget", "https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz")
exec_wait("gunzip", c("-vf", "Mus_musculus.GRCm39.113.gtf.gz"))

# - mm39 ensembl variant file available here https://ftp.ensembl.org/pub/release-112/variation/vcf/mus_musculus/
exec_wait("wget", "https://ftp.ensembl.org/pub/release-113/variation/vcf/mus_musculus/mus_musculus.vcf.gz")
# decompress this ^^^ one later! You will be prompted

# This makes it obvious in the console that a long job is complete
# Not cited since github link makes that redundant
remotes::install_github("richfitz/rfiglet", upgrade = FALSE)
library(rfiglet)
figlet("FIGLET RUNS")

# BEFORE MOVING FORWARD
conflicted::conflict_scout()
# CHECK THE INFILES WERE DOWNLOADED AND DECOMPRESSED PROPERLY


cl <- makeCluster(parallel::detectCores())
dir_path <- "data"
dirs <- dir(path = dir_path, full.names = T)
infiles <- dir(path = dirs, pattern="*fastq*", full=TRUE)
infiles # ALWAYS RUN AND READ THIS

adapters_real <- NULL
# Illumina TruSeq Adapters (since project in example used NextSeq 500)
# Variable is named "real" since the adapters are copied from the toolkit,
# rather than computationally inferred
adapters_real[1] <- "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
adapters_real[2] <- "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

td <- tempdir()
trimfiles <- str_replace(basename(infiles),".fastq","_trimmed.fq")
trimfiles <- file.path(td,trimfiles) 
trimfiles
outfiles <- str_replace(basename(infiles),".fastq","_processed.fq")
outfiles <- file.path(td, outfiles) 
outfiles

# CHECK THE INFILES WERE DOWNLOADED AND DECOMPRESSED PROPERLY
infiles
infiles#!!!
i <- 1
for (i in seq(1,length(infiles),2)) 
{remove_adapters(file1=infiles[i], #Rbowtie2
                 paste0("--threads ", parallel::detectCores()), 
                 "--minquality 5", 
                 file2=infiles[i+1],
                 adapter1 = adapters_real[1], 
                 adapter2 = adapters_real[2],
                 output1=trimfiles[i],
                 output2=trimfiles[i+1],
                 basename=file.path(td,"reads.base"),
                 overwrite=TRUE)}

cl

preP <- preprocessReads(filename = trimfiles[seq(1,length(infiles),2)],
                        filenameMate = trimfiles[seq(2,length(infiles),2)],
                        outputFilename = outfiles[seq(1,length(outfiles),2)],
                        outputFilenameMate = outfiles[seq(2,length(outfiles),2)],
                        complexity=0.3, # from fastp defaults
                        truncateEndBases = 3,
                        minLength = 14, 
                        nBases = 2,
                        clObj = cl)
figlet("done")



### MANUALLY GET YOUR METADATA IF YOU HAVEN'T #####################
raw_metadata <- read.csv("SraRunTable.csv")
View(raw_metadata)
Phenotype <- raw_metadata$genotype.variation
metadata_neater <- tibble(raw_metadata$Run,
                          raw_metadata$AvgSpotLen,
                          Phenotype,
                          raw_metadata$genotype.variation,
                          raw_metadata$source_name)
metadata_neater
metadata_cols <- c(Accession ="raw_metadata$Run",
                   AvgSpotLen = "raw_metadata$AvgSpotLen",
                   Genotype = "raw_metadata$genotype.variation", 
#                   Phenotype = "raw_metadata$genotype.variation",
                   Source = "raw_metadata$source_name")
metadata_neater <- rename(metadata_neater,all_of(metadata_cols))
metadata_neater
metadata_neater <- mutate(metadata_neater,
                          across(c(Genotype, Phenotype, 
                                   Source), as.factor)
)
levels(metadata_neater$Phenotype) <- c("WT","WT","CBP_KO","Triallelic_KO","p300_KO")
levels(metadata_neater$Genotype) <- c("WT_p300flfl","WT_CBPflfl","CBP_KO","Triallelic_KO","p300_KO")

metadata_neater

samples_rnaSeq <- tibble::tibble(`FileName1`= outfiles[seq(1,length(outfiles),2)],
                                      `FileName2`= outfiles[seq(2,length(outfiles),2)],
                                      `SampleName`= paste0(metadata_neater$Accession,
                                                           "_",
                                                           metadata_neater$Genotype)
)
samples_rnaSeq
write_tsv(samples_rnaSeq, "samples_rnaSeq.txt")
sampleFile <- "samples_rnaSeq.txt"


genomeFile <- "Mus_musculus.GRCm39.dna_rm.primary_assembly.fa"
annotFile <- "Mus_musculus.GRCm39.112.gtf"
dir.create(file.path(getwd(), "alignment_out"))

# REMEMBER THAT IF YOU WANT A FRESH ALIGNMENT, NOT ONLY alignment_out, BUT ALSO
#                                      path/to/alignment.Rhisat2 MUST BE DELETED

# This does annoyingly name the alignment files after the first of the 
# pair of read files, although qAlign does in fact take both inputs into
# account when generating its output
alignment <- qAlign(sampleFile,#DANGER: RUNNING THE qAlign HELP PAGE EXAMPLE WILL CREATE DIFFERENT SAMPLE FILE 
                      genome = genomeFile,
                      splicedAlignment = T,
                      aligner = "Rhisat2",
                      geneAnnotation = annotFile,
                      clObj = cl,
                      paired = "fr",
                      alignmentsDir = "alignment_out"
)
bams_dir <- "alignment_out" 
bams <- dir(path = bams_dir, pattern = ".bam$", full.names = T)

gatk <- "path/to/gatk-version/gatk" # modify to fit your installation path
exec_wait(gatk, "--help")
exec_wait(gatk, c("AddOrReplaceReadGroups", "--help"))
# Text Inside ViewBamHeader.sh = samtools view $1 | head -n $2
# Ignore the following R console output:
#samtools view: writing to standard output failed: Broken pipe
#samtools view: error closing standard output: -1
exec_wait("./shell/ViewBamHeader.sh", c(bams[1], 2))

# These publicly available data didn't contain metadata on library and lane,
# so the sequencing libraries are named SampleName-1 and the flowcell lanes
# LibraryName-1 (SampleName-1-1)
for (i in 1:length(bams)) {
  exec_wait("samtools", c("addreplacerg", 
                          c("-r", "PL:ILLUMINA"), # Seq Platform
                          c("-r", 
                            paste0("SM:GSM27058",(i+79))), 
                          c("-r", paste0("LB:GSM27058",(i+79),"-1")), 
                          c("-r", paste0("ID:GSM27058",(i+79),"-1-1")), 
                          c("-m", "overwrite_all"), 
                          c("-o", paste0(td,"/", substr(bams[i],15,24),"_RG.bam")),
                          c("--threads", (detectCores()-1)), # number in addition to main thread
                          bams[i]
  )
  )
}
figlet("done")
alignment
qQCReport(rnaSeq_test, pdfFilename = "quasR_QCreport.pdf", clObj = cl)
figlet("done")
alignmentStats(alignment)

########################## SNP CALLING #########################################

exec_wait(gatk, c("CreateSequenceDictionary", 
                  c("-R", genomeFile)
                  )
          )

td <- tempdir()
RGbams <- dir(path = td, pattern = "*RG.bam$", full.names = T)
sortedbams <- str_replace(basename(RGbams),"RG.bam","RG_sorted.bam")
sortedbams <- file.path(td,sortedbams) 

exec_wait("./shell/ViewBamHeader.sh", c(RGbams[1], 25))

dir.create(file.path(getwd(), "VC_in"))
vc_in <- "VC_in" 

# requires samtools installation
for (i in 1:length(RGbams)) {
  exec_wait("samtools", c("sort", 
                          RGbams[i], 
                          c("-o", sortedbams[i]), 
                          c("--threads", detectCores())
  )
  )
}

# configure multiple inputs
merge_in<-as_tibble(paste0("-I ", sortedbams[1:length(sortedbams)], collapse = " "))
write_file(merge_in[[1]], "VC_in/merge_inputs.txt")
merge_inputs <- "VC_in/merge_inputs.txt"

# GATK can be found at https://github.com/broadinstitute/gatk/releases

exec_wait(gatk, std_in = sortedbams[1:length(sortedbams)],
          c("MergeSamFiles",
            c("--arguments_file", merge_inputs),
            c("-O", paste0(vc_in, "/merged_sorted.bam")),
            c("--SORT_ORDER", "queryname"), # for Spark marking of duplicates
            c("--USE_THREADING", "true"),
            c("-R", genomeFile), 
            c("--TMP_DIR", td)
          )
)

i <- 1
duplicates_out <- file.path(vc_in,"marked_duplicates.bam")

# If you run into chmod exception and can't write your output without sudo rstudio,
# run the line commented out below:
# gatk <- "/path/from/root/to/gatk-4.6.0.0/gatk" # modify to fit your installation version
# Alternatively in Linux cmd line: 
# sudo gatk MarkDuplicatesSpark  --tmp-dir path/to/choice/ -I alignment_out/merged_sorted.bam -O alignment_out/marked_duplicates.bam -M alignment_out/marked_dups_args.txt #Other options
exec_wait(gatk, 
          c("MarkDuplicatesSpark", 
            c("-I", paste0(vc_in, "/merged_sorted.bam")), 
            c("-O", duplicates_out),
            c("-M", paste0(vc_in, "/marked_dups_metrics.txt")),
            c("--tmp-dir", td) 
            
          )
)
gatk <- "path/to/gatk-4.6.0.0/gatk" # modify to fit your installation version


# Using slurm
# Run pipeline SplitNCigarReads, Recalibrate base qualities, analyse covariates
exec_wait("./shell/vCPrep_GRCm39_slurm.sh", c(td))


# Without using slurm
# SplitNCigarReads 
variantFile <- "mus_musculus.vcf.gz"
exec_wait(gatk, c("IndexFeatureFile", 
                  c("-I", variantFile))) # Documentation says this should be .gz

exec_wait("gunzip", c("-vf", "mus_musculus.vcf.gz"))
variantFile <- "mus_musculus.vcf" #Documentation says later steps use decompressed

SplitNCigar_out <- file.path(vc_in, "SplitNCigarReads.bam")
SeqDict <- str_replace(genomeFile, ".fa", ".dict") # Created in line 281 as of this comment

exec_wait(gatk,
          c("SplitNCigarReads",
            c("-R", genomeFile),
            c("-I", duplicates_out),
            c("-O", SplitNCigar_out),
            c("--sequence-dictionary", SeqDict),
            c("--tmp-dir", td)))



exec_wait("./shell/ViewBamHeader.sh", c(SplitNCigar_out, 25))


# Recalibrate Bases, analyse covars
exec_wait(gatk, c("BaseRecalibrator", # Spark version in beta. not working for me
                  c("-I", SplitNCigar_out),
                  c("-R", genomeFile),
                  c("--known-sites", variantFile),
                  c("-O", paste0(vc_in,"/baseRecal.table")),
                  c("--tmp-dir", td)
  )
)

exec_wait(gatk, c("ApplyBQSR", 
                  c("-I", SplitNCigar_out),
                  c("-R", genomeFile),
                  c("-bqsr", paste0(vc_in,"/baseRecal.table")),
                  c("-O", paste0(vc_in,"/baseRecal.bam")),
                  c("--tmp-dir", td)
)
)

exec_wait(gatk, c("AnalyzeCovariates", 
                  c("-bqsr", paste0(vc_in,"/baseRecal.table")),
                  c("-csv", paste0(vc_in,"/AnalyseCovariates.csv"))
))



# REMEMBER TO WRITE IN THE SAMPLE NAME for new projects
# Once you get to second iteration without problems, call it a day/weekend
# /go do some reading/etc.
for (i in 80:94) { # Variable
  exec_wait(gatk, c("HaplotypeCaller", 
                    c("-I", paste0(vc_in, "/baseRecal.bam")),
                    c("--sample-name", paste0("GSM27058", i)), # Variable
                    c("-R", genomeFile),
                    c("-O", paste0("/VC_out/",
                                   "GSM27058",i, # Variable
                                   ".g.vcf.gz")),
                    c("-ERC", "GVCF"),
                    c("--native-pair-hmm-threads", detectCores()),
                    c("--tmp-dir", td)
  )
  )
}

exec_wait("./shell/gunzip_vcf.sh")

# Make vcfR object from any sample
#vcf_file <- read.vcfR("VC_out/GSM2705880.g.vcf", verbose = FALSE)
#vcf_file <- read.vcfR("VC_out/GSM2705881.g.vcf", verbose = FALSE)
#vcf_file <- read.vcfR("VC_out/GSM2705882.g.vcf", verbose = FALSE)
#vcf_file <- read.vcfR("VC_out/GSM2705883.g.vcf", verbose = FALSE)
#vcf_file <- read.vcfR("VC_out/GSM2705884.g.vcf", verbose = FALSE)
#vcf_file <- read.vcfR("VC_out/GSM2705885.g.vcf", verbose = FALSE)
#vcf_file <- read.vcfR("VC_out/GSM2705886.g.vcf", verbose = FALSE)
#vcf_file <- read.vcfR("VC_out/GSM2705887.g.vcf", verbose = FALSE)
#vcf_file <- read.vcfR("VC_out/GSM2705888.g.vcf", verbose = FALSE)
#vcf_file <- read.vcfR("VC_out/GSM2705889.g.vcf", verbose = FALSE)
#vcf_file <- read.vcfR("VC_out/GSM2705890.g.vcf", verbose = FALSE)
#vcf_file <- read.vcfR("VC_out/GSM2705891.g.vcf", verbose = FALSE)
#vcf_file <- read.vcfR("VC_out/GSM2705892.g.vcf", verbose = FALSE)
#vcf_file <- read.vcfR("VC_out/GSM2705893.g.vcf", verbose = FALSE)
#vcf_file <- read.vcfR("VC_out/GSM2705894.g.vcf", verbose = FALSE)

###################featureCounts##############################################


fcSE_test <- featureCounts(bams, 
                           #annot.inbuilt = "mm10", # (Eg)
                           annot.ext = annotFile, isGTFAnnotationFile = TRUE,
                           isPairedEnd = T, nthreads = parallel::detectCores(),
                           verbose=T)

colnames(fcSE_test$counts) <- paste0(metadata_neater$Accession,
                                     "_",
                                     metadata_neater$Genotype)

metadata_neater$Name <- colnames(fcSE_test$counts)

###########################DESeq2############################################
# (Love et al. @ https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html)
# (and, of course, their references)
dds <- DESeqDataSetFromMatrix(countData = fcSE_test$counts,
                              design = ~ Phenotype,
                              colData = metadata_neater
)

nrow(dds)
assays(dds)

smallestGroupSize <- 3 
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
nrow(dds)

vsd <- vst(dds, blind = FALSE) 
head(assay(vsd), 3) # normalised with variance stabilising transformation
colData(vsd)
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3) # normalised with rlog transformation
head(assay(dds),3) # not normalised

dds <- estimateSizeFactors(dds)
df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)) %>% mutate(transformation = "rlog"),
  #as_data_frame(assay(dds)) %>% mutate(transformation = "none") # visualising effect of tranformations. should break plot?
)
lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)
ggplot(df, aes(x = SRR5833394_WT_CBPflfl, y = SRR5833395_p300_KO)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  


sampleDists <- dist(t(assay(vsd)))
sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- vsd@colData@listData$Name
colnames(sampleDistMatrix) <- vsd@colData@listData$Name
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colours)


poisd <- PoissonDistance(t(counts(dds)))

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- dds@colData@listData$Name
colnames(samplePoisDistMatrix) <- dds@colData@listData$Name
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colours)

plotPCA(vsd, intgroup = "Name")
plotPCA(vsd, intgroup="Genotype")
pcaData <- plotPCA(vsd, intgroup = "Genotype", returnData=T)
pcaData <- pcaData[order(pcaData$Genotype),]

#Vector storing percentages of the variance explained by PC1 and then PC2
percentVar <- round(100 * attr(pcaData, "percentVar"))

# (Generalised Linear Model PCA 
gpca <- glmpca(counts(dds), L=2) # notice it uses the raw
gpca.dat <- gpca$factors
gpca.dat$Accession <- dds@colData@listData$Accession
gpca.dat$Genotype <- dds@colData@listData$Genotype


ggplot(gpca.dat, aes(x = dim1, y = dim2, 
                     colour = Genotype #Accession,
                     #shape = Genotype
)) +
  geom_point(size =3) + 
  coord_fixed() + 
  ggtitle("glmpca - Generalized PCA")

# MDS Plot - Useful for multidimensional scaling like PCAs but for
# when you have a matrix of distances rather than data
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`,
                colour = Genotype#Accession,
                #shape = Genotype
)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")

mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, 
                    colour = Genotype#Accession,
                    #shape = Genotype
)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDists")

# Run DiffEx pipeline, changing the dds object
dds <- DESeq(dds)

ds_out <- results(dds)
ds_out 
# Column metadata for results tables
mcols(ds_out, use.names = TRUE)

# You can get specifically chosen contrasts by changing the last 2 arguments
# supplied to contrast - (DATA EXPLORATION, DON'T JUST RUN)
#ds_out_WTvWT <- results(dds, contrast = c("Genotype", "WT_p300flfl", "WT_CBPflfl"))
ds_out_p300vWT <- results(dds, contrast = c("Phenotype", "p300_KO", "WT"))
ds_out_CBPvWT <- results(dds, contrast = c("Phenotype", "CBP_KO", "WT"))
ds_out_TAvWT <- results(dds, contrast = c("Phenotype", "Triallelic_KO", "WT"))


summary(ds_out_p300vWT)
summary(ds_out_CBPvWT)
summary(ds_out_TAvWT)

ds_out.05 <- results(dds, alpha = 0.05)
table(ds_out.05$padj < 0.05)

ds_outLFC1 <- results(dds, lfcThreshold=1)
table(ds_outLFC1$padj < 0.1)

sum(ds_out_p300vWT$pvalue < 0.05, na.rm=TRUE) # p < 0.05...
sum(ds_out_CBPvWT$pvalue < 0.05, na.rm=TRUE) # p < 0.05...
sum(ds_out_TAvWT$pvalue < 0.05, na.rm=TRUE) # p < 0.05...
sum(!is.na(ds_out$pvalue))#           ...out of this many genes

# If a false discovery rate of 10% is acceptable, then a 
# padj < 0.1 is significant
sum(ds_out_p300vWT$padj < 0.05, na.rm=TRUE)
sum(ds_out_CBPvWT$padj < 0.05, na.rm=TRUE)
sum(ds_out_TAvWT$padj < 0.05, na.rm=TRUE)


# Subset your genes that have this significant padj value
dsSig_c1 <- subset(ds_out_p300vWT, padj < 0.05)
head(dsSig_c1[ order(dsSig_c1$log2FoldChange), ]) # most upreg'd
head(dsSig_c1[ order(dsSig_c1$log2FoldChange, decreasing = TRUE), ]) # most downreg'd
# maybe you prefer looking at biggest changes together
dsSig_c1[ order(dsSig_c1$log2FoldChange), ]

# Subset your genes that have this significant padj value
dsSig_c2 <- subset(ds_out_CBPvWT, padj < 0.05)
head(dsSig_c2[ order(dsSig_c2$log2FoldChange), ]) # most upreg'd
head(dsSig_c2[ order(dsSig_c2$log2FoldChange, decreasing = TRUE), ]) # most downreg'd
# maybe you prefer looking at biggest changes together
dsSig_c2[ order(dsSig_c2$log2FoldChange), ]

# Subset your genes that have this significant padj value
dsSig_c3 <- subset(ds_out_TAvWT, padj < 0.05)
head(dsSig_c3[ order(dsSig_c3$log2FoldChange), ]) # most upreg'd
head(dsSig_c3[ order(dsSig_c3$log2FoldChange, decreasing = TRUE), ]) # most downreg'd
# maybe you prefer looking at biggest changes together
dsSig_c3[ order(dsSig_c3$log2FoldChange), ]

###############################PLOTTING#########################################

# You can try the whole plotting section using independent hypothesis weighting
# to filter the dataset after you have looked at the unweighted results
#ds_out <- results(dds, filterFun=ihw)

# What is the most significantly differentially-expressed gene wrt Genotype
# (Change to a contrast of choice if you want, by swapping "ds_out" between options)
topGene <- rownames(ds_out_TAvWT)[which.min(ds_out_TAvWT$padj)]

geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("Genotype","Accession"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = Genotype, y = count, colour = Accession)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)
# (for grouped comparisons like same mouse with and without treatment)
ggplot(geneCounts, aes(x = Genotype, y = count, 
                       colour = Accession,
                       group = Accession
)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()

resultsNames(dds)

ds_out_TAvWT <- lfcShrink(dds, 
                          coef="Phenotype_Triallelic_KO_vs_WT", # or anything from 
                          type="apeglm"                       # running resultsNames()
)
ds_out_CBPvWT <- lfcShrink(dds, 
                           coef="Phenotype_CBP_KO_vs_WT", # or anything from 
                           type="apeglm"                       # running resultsNames()
)
ds_out_p300vWT <- lfcShrink(dds, 
                            coef="Phenotype_p300_KO_vs_WT", # or anything from 
                            type="apeglm"                       # running resultsNames()
)

# Change to whatever comparison you see fit
plotMA(ds_out_TAvWT, ylim = c(-5, 5))
with(ds_out_TAvWT[topGene, ], { 
  points(baseMean, log2FoldChange, col="red", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="red")
})
# to see without shrinkage, enter in console:
# plotMA(results(dds, name = "Phenotype_Triallelic_KO_vs_WT"), ylim = c(-5,5))

hist(ds_out_TAvWT$pvalue[ds_out$baseMean > 1],
     main = "Histogram of p values for genes with mean normalized count larger than 1.",
     breaks = 0:20/20,
     col = "grey50", border = "white")

###################CLUSTERING###########################################
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)

mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("Phenotype",#"Genotype"#(Choice)
                                       "Accession"
)])

pheatmap(mat, annotation_col = anno)

qs <- c(0, quantile(ds_out_TAvWT$baseMean[ds_out_TAvWT$baseMean > 0], 0:6/6))
bins <- cut(ds_out_TAvWT$baseMean, qs)
levels(bins) <- paste0("~", round(signif((qs[-1] + qs[-length(qs)])/2, 2)))
fractionSig <- tapply(ds_out_TAvWT$pvalue, bins, function(p)
  mean(p < .05, na.rm = TRUE))
barplot(fractionSig, xlab = "mean normalized count",
        ylab = "fraction of small p values")

#ds_out.ihw <- results(dds, filterFun=ihw)


################################ANNOTATE/EXPORT################################

# Please refer to Youtube@Bioinformagician "BioMart" video for ins-and-outs
#inputs
ens.ids <- substr(rownames(ds_out), 1, 19) #for mouse, in case of version suffix
ens.ids

listEnsembl()
listEnsemblArchives()
listMarts(host="https://may2024.archive.ensembl.org") # Release 112 (corresponds to downloads)
listMarts()
#ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", host="https://nov2020.archive.ensembl.org")
ensembl <- useEnsembl(biomart = "genes")
ensembl112 <- useEnsembl(biomart = 'genes', 
                         dataset = 'mmusculus_gene_ensembl',
                         version = 112)
ens.sets <- listDatasets(ensembl112) 
ens.conn <- useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl"
                    ,host="https://may2024.archive.ensembl.org"
)

attr <- listAttributes(ens.conn)
filters <- listFilters(ens.conn)

#"go_id"
idVSname <- getBM(attributes = c("ensembl_gene_id","external_gene_name"),
                  filters = "ensembl_gene_id",
                  values = ens.ids,
                  mart = ens.conn,
                  verbose=T)
entrezVSname <- getBM(attributes = c("ensembl_gene_id","entrezgene_id"),
                      filters = "ensembl_gene_id",
                      values = ens.ids,
                      mart = ens.conn,
                      verbose=T)

bm_out <- getBM(attributes = c("ensembl_gene_id",
                               "external_gene_name", 
                               "entrezgene_id"#,
                               #"go_id"
                               ),
                filters = "ensembl_gene_id",
                values = ens.ids,
                mart = ens.conn,
                verbose=T
                )

nrow(ds_out)
bm_out <- bm_out[!duplicated(bm_out$ensembl_gene_id),]

figlet("COMPLETE!")

ds_out_p300vWT$symbol <- bm_out$external_gene_name
ds_out_p300vWT$entrez <- bm_out$entrezgene_id
ds_out_CBPvWT$symbol <- bm_out$external_gene_name
ds_out_CBPvWT$entrez <- bm_out$entrezgene_id
ds_out_TAvWT$symbol <- bm_out$external_gene_name
ds_out_TAvWT$entrez <- bm_out$entrezgene_id


ds_out_p300vWT <- ds_out_p300vWT[order(ds_out_p300vWT$padj),]
head(ds_out_p300vWT)
ds_out_CBPvWT <- ds_out_CBPvWT[order(ds_out_CBPvWT$padj),]
ds_out_TAvWT <- ds_out_TAvWT[order(ds_out_TAvWT$padj),]

dir.create(file.path(getwd(), "DE_out"))
write.csv(ds_out_p300vWT, file = "DE_out/p300vsWT.csv")
write.csv(ds_out_CBPvWT, file = "DE_out/CBPvsWT.csv")
write.csv(ds_out_TAvWT, file = "DE_out/TAvsWT.csv")

#######################ENRICHMENT##########ANALYSIS##########################
# Some documentation at:
# https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/ 

# Which enrichment analyses? GO, KEGG, (more to come)
original_gene_list <- ds_out_TAvWT$log2FoldChange
names(original_gene_list) <- rownames(ds_out_TAvWT)
# omit any NA values 
gene_list<-na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

entrez_list <- ds_out_TAvWT$log2FoldChange
names(entrez_list) <- ds_out_TAvWT$entrez
entrez_list <- na.omit(entrez_list)
entrez_list = sort(entrez_list, decreasing = TRUE)

enrichment_out <- gseGO(geneList=gene_list, 
                        ont ="ALL", 
                        keyType = "ENSEMBL", 
                        #nPerm = 10000, # got warnings saying it was deprecated
                        minGSSize = 3, 
                        maxGSSize = 800, 
                        pvalueCutoff = 0.05, 
                        verbose = TRUE, 
                        OrgDb = org.Mm.eg.db::org.Mm.eg.db, 
                        #pAdjustMethod = "none"
                        eps = 0 # advised to by output of documented version
)


dotplot(enrichment_out, showCategory=10, split=".sign") + facet_grid(.~.sign)

enrich2 <- pairwise_termsim(enrichment_out)
emapplot(enrich2, showCategory = 10)

# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(enrichment_out, categorySize="pvalue", foldChange=gene_list, showCategory = 3)

ridgeplot(enrichment_out) + labs(x = "enrichment distribution")

# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(enrichment_out, by = "all", 
         title = enrichment_out$Description[2], # change this to explore output 
         geneSetID = 1)

terms <- enrichment_out$Description[1:5] # Change the numbers but don't increase
#                                          range length
pmcplot(terms, 2010:substr(paste0(Sys.Date()),1,4), proportion=FALSE)



kk <- gseKEGG(geneList     = entrez_list,
              organism     = 'mmu',
              minGSSize    = 120,
              pvalueCutoff = 0.05)
head(kk)


browseKEGG(kk, 'mmu05166')
