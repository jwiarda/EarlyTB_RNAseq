### Using R version 3.5.0

#############################
# Set up required packages: #
#############################

### Install required packages:
source("https://bioconductor.org/biocLite.R") # will install most updated compatible versions
biocLite("Rsubread")
biocLite("edgeR")

install.packages("dplyr")
install.packages("tidyverse")

### Load required packages: 
library(edgeR) # using version 1.30.0
library(Rsubread) # using version 3.22.0
library(dplyr) # using version ________
library(tidyverse)

####################################################
# Create table of gene counts using featureCounts: #
####################################################

### Create gene count table with all sample bam files:

# Create a vector of all bam files:
bamFiles <- list.files(path = "/Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/STAR_alignments", # only look for files in STAR_alignment directory
                       pattern = "Aligned.out.bam", # look for files containing Aligned.out.bam in name
                       full.names = TRUE) # lists absolute path of each file
length(bamFiles) # output should be equal to number of bam files
head(bamFiles, n = 3) # Check that bam file names are correct

# Create a directory to save featureCount output files in:
dir.create("/Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/featureCounts")

# Run featureCounts on all bam files:
geneCounts <- c(featureCounts(files = bamFiles,
                              annot.ext = "/Users/jayne.wiarda/RNASeq/Bos_taurus_UMD3.1.1/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff", # location of annotation file
                              requireBothEndsMapped = TRUE, # both ends of reads must map
                              countChimericFragments = FALSE, # do not count chimeric fragements
                              GTF.attrType = "Name", # use the Name attribute from GFF file to group features
                              GTF.featureType = "gene", # use gene to define features from GFF file
                              isGTFAnnotationFile = TRUE, # file is in gtf/gff format
                              isPairedEnd = TRUE, # paired-end data
                              strandSpecific = 2, # reversely stranded
                              nthreads = 4, # use 4 threads
                              reportReads = "BAM", # stores featureCounts output in bam file in current directory
                              reportReadsPath = "/Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/featureCounts")) # specifies where to store bam output file

# Check gene count table:
dim(geneCounts$counts) # dimensions should be 26,792 (# of genes in annotation file) and 90 (for number of samples/depths)
head(geneCounts$counts, n=10) # should see list of gene names on left and a column on right with bam file name and gene counts >= 0
head(geneCounts$stat, n = 1) # shows number of reads that were assigned to features (genes)... expected to be ~75-85% of total counts

### Export count data as tab-delimited files:

# Create tab-delimited file of total counts table:
write.table(x=data.frame(geneCounts$annotation[,c("GeneID")],
                         geneCounts$counts,stringsAsFactors=FALSE),
            file="/Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/featureCounts/82_Cattle_gene_counts.txt",
            quote=FALSE,sep="\t",
            row.names=FALSE)

# Create tab-delimited file of total counts statistics:
write.table(x=data.frame(geneCounts$stat),
            file="/Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/featureCounts/82_Cattle_gene_count_stats.txt",
            quote=FALSE,sep="\t",
            row.names=FALSE)
