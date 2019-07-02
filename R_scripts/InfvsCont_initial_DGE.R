#################################################
##                                             ##
## 82-Cattle-16 Project RNA-seq DGE analysis   ##
## TB infected vs uninfected cattle            ##
## Unstimulated whole blood direct collections ##
##                                             ##
#################################################

#############################
# Set up required packages: #
#############################

### Install required packages:
#source("https://bioconductor.org/biocLite.R") # will install most updated compatible versions
#biocLite("Rsubread")
#biocLite("edgeR")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("edgeR")

#install.packages("dplyr")
#install.packages("tidyverse")
#install.packages("VennDiagram")
#install.packages("gplots")
#install.packages("pca3d")
#install.packages("reshape2")

### Load required packages: 
library(edgeR) # using version 1.30.0
library(Rsubread) # using version 3.22.0
library(dplyr) # using version ________
library(tidyverse)
library(VennDiagram)
library(gplots)
library(pca3d)
library(reshape2)

########################################
# Load and organize count data into R: #
########################################

### Upload count data and target file:

# Load target file into R:
target <-read_tsv("/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/82_cattle_target.txt")
# target # view full target file
target <- target[with(target, order(animal_id, weeks_post_infection)),] # reorder target file
# view(target)

# Load counts data back into R:
counts <- read.delim("/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/featureCounts/82_Cattle_genecounts.txt")
colnames(counts) # make sure columns 2-37 are in same order as entries in 'label' column of target
head(target, n = 10)

#Create another entry in target file that combines the trt_rep and weeks_post_infection:
target$label <- paste(target$trt_rep,target$weeks_post_infection, sep="_")
target$label <- paste(target$label,"wpi", sep="")

# Rename columns according to target file:
names(counts) <- c("Gene", target$label[1:39])

# Reorder columns of count table by trt_rep, then by collection timepoint:
target <- target[with(target, order(trt_rep, weeks_post_infection)),]
target$label
sortList <- target$label[1:39]
counts <- counts[c("Gene", sortList)]
colnames(counts)

########################
# EdgeR DGE Analysis:: #
########################

### Samples are from Wk 0, 4, & 10 post-infection timepoints TB-infected & uninfected cattle

# Create additional column of data in target file:
target$treatment <- paste((substr(target$trt_rep, 1, 3)), target$weeks_post_infection, sep="_") # defines the treatment group and date of collection
target$treatment <- paste(target$treatment,"wpi", sep="")
target$treatment

# Define treatment factor:
group <- factor(target$treatment)
group

# Define the experimental design and generalized linear model fit:
design <- model.matrix(~ 0 + group) # defining using group-means parameterization
design

# Create the DGE object:
DGE <- DGEList(counts = counts[,2:40], genes = counts[,1], group = group)
DGEold <- DGE # use this later

# Define the experimental design and generalized linear model fit:
design <- model.matrix(~ 0 + group) # defining using group-means parameterization
design

### Calculate the normalization factors (TMM method):

DGE <- calcNormFactors(DGE)
DGE$samples # shows adjusted normalization factors

### Obtain raw gene count densities for quality check:
rawDGE <- log10(DGE$counts[,1:ncol(DGE$counts)] + 1) # log transform the raw data

### Filter out lowly expressed genes:

# Filter out genes with less than 0.5 counts per million (cpm) in at least 5 samples:
keep <- rowSums(cpm(DGE) > 1) >= 5

table(keep) # shows how many genes are kept vs. filtered out

DGE <- DGE[keep, , keep.lib.sizes = FALSE] # removes lowly expressed genes from DGE object

write.table(DGE$genes, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/infVsCont_filtered_genes.txt", sep = "\t")

head(DGE$counts, n = 5) # check that lowly expressed genes are filtered out

DGE$samples # shows new library sizes

# Plot the old vs new library sizes after filtering out low genes:
mean(DGEold$samples$lib.size)
mean(DGE$samples$lib.size) # average library size remains similar

### Plot raw and post-filtering gene count densities for quality check:

filteredDGE <- log10(DGE$counts[,1:ncol(DGE$counts)] + 1) # log transform the filtered data

plot(density(rawDGE[,1]),
     main = "Gene Count Densities Unfiltered vs Filtered Reads",
     lty = 1,
     lwd = 2,
     xlab = "log10 gene counts",
     ylab = "Density",
     col = "blue4",
     ylim = c(0.0,1.0))
lines(density(filteredDGE[,1]), 
      col = "lightskyblue",
      lty = 1,
      lwd = 2)
legend("topright", 
       legend = c("Unfiltered Reads", "Filtered Reads"), 
       col = c("blue4", "lightskyblue"), 
       lty = 1, 
       lwd =2)

### Create MDS plot of all samples:
col <- c("blue", "blue", "blue", "orange", "orange", "orange")
pch <- c(16, 17, 15, 16, 17, 15)

plotMDS(DGE,
        top = 500,
        col = col [as.factor(group)],
        pch = pch [as.factor(group)],
        cex = 2)
legend("bottomright", 
       legend=c("con0wpi","con4wpi","con10wpi","inf0wpi","inf4wpi","inf10wpi"), 
       pch = c(16,15,17,16,15,17), 
       col = col, 
       cex = 0.8,
       pt.cex = 1.2)

# Estimate common dispersion (required to estimate tagwise dispersion):
DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) # output will provide common dispersion & biological coefficient of variation (BCV)

# Estimate trended dispersion:
DGE <- estimateGLMTrendedDisp(DGE, design)

# Estimate tagwise dispersion (will be used for DE analysis):
DGE <- estimateGLMTagwiseDisp(DGE, design)

### Plot the estimated dispersions:
plotBCV(DGE)

### Test for DGE using GLM:

# Define GLM:
fit <- glmFit(DGE, design)

### Define comparison contrasts & obtain DEGs:
colnames(design) # locate which columns the groups of interest for comparison are in
# for reference:
# genes with positive logFC will be upregulated in groups with contrast = 1 compared to group with -1 contrast
# genes with negative logFC will be downregulated in groups with contrast = 1 compared to group with -1 contrast

### Week 0 infected vs control group (preinfection timepoint):
# Perform likelihood ratio test & identify DEGs:
wk0lrt <- glmLRT(fit, contrast = c(-1, 0, 0, 1, 0, 0)) # perform likelihood ratio test 
summary(wk0DGE <- decideTestsDGE(wk0lrt)) # reports upregulated, downregulated, and not DE genes
wk0topTags <- topTags(wk0lrt, n = Inf) # table of genes sorted by FDR value (FDR < .05 means DE)
wk0DE <- subset(wk0topTags$table, FDR < .05) # create a new table of only DE (FDR < .05) genes
head(wk0DE)
# Create plot smear:
wk0DEnames <- rownames(wk0DE)
plotSmear(wk0lrt, de.tags = wk0DEnames, cex = .75, ylim = c(-6,6))
abline(h = c(-1, 1), col = "blue") 
write.table(wk0DE, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/InfvsCont_0wpi_allDE.txt", sep = "\t")

### Week 4 post-infection infected vs uninfected:
wk4lrt <- glmLRT(fit, contrast = c(0, 0, -1, 0, 0, 1))
summary(wk4DGE <- decideTestsDGE(wk4lrt)) # reports upregulated, downregulated, and not DE genes
wk4topTags <- topTags(wk4lrt, n = Inf) # table of genes sorted by FDR value (FDR < .05 means DE)
wk4DE <- subset(wk4topTags$table, FDR < .05) # create a new table of only DE (FDR < .05) genes
head(wk4DE, n = 10)
# Create plot smear:
wk4DEnames <- rownames(wk4DE)
plotSmear(wk4lrt, de.tags = wk4DEnames, cex = 0.75, ylim = c(-6,6))
abline(h = c(-1, 1), col = "blue") 
write.table(wk4DE, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/InfvsCont_4wpi_allDE.txt", sep = "\t")

### Week 10 post-infection infected vs uninfected:
wk10lrt <- glmLRT(fit, contrast = c(0, -1, 0, 0, 1, 0))
summary(wk10DGE <- decideTestsDGE(wk10lrt)) # reports upregulated, downregulated, and not DE genes
wk10topTags <- topTags(wk10lrt, n = Inf) # table of genes sorted by FDR value (FDR < .05 means DE)
wk10DE <- subset(wk10topTags$table, FDR < .05) # create a new table of only DE (FDR < .05) genes
head(wk10DE, n = 10)
# Create plot smear:
wk10DEnames <- rownames(wk10DE)
plotSmear(wk10lrt, de.tags = wk10DEnames, cex = 0.75, ylim = c(-6,6))
abline(h = c(-1, 1), col = "blue") 
write.table(wk10DE, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/InfvsCont_10wpi_allDE.txt", sep = "\t")

### Heat map of DE genes at weeks 4 and 10:
wk4target <- subset(target, weeks_post_infection == 4) # select for only Wk. 4 samples
# View(wk4target) # confirm only Wk. 4 samples are present
wk4counts <- counts[,wk4target$label] # collect only count data for the subsetted samples
colnames(wk4counts)
wk4counts$genes <- counts[,1]
wk4DEcounts <- merge(wk4counts, as.data.frame(wk4DE), by="genes", type="inner")
colnames(wk4DEcounts)
wk4target$treatment <- paste((substr(wk4target$trt_rep, 1, 3)), wk4target$weeks_post_infection, sep="_") # defines the treatment group and date of collection
wk4group <- factor(wk4target$treatment)
wk4group
colnames(wk4DEcounts)
wk4DGE <- DGEList(counts = wk4DEcounts[,2:14], genes = wk4DEcounts[,1], group = wk4group)
logcounts <- cpm(wk4DGE,log=TRUE)
var_genes <- apply(logcounts, 1, var)
head(var_genes)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:1060] # looking at all the DE genes
head(select_var)
highly_variable_lcpm <- logcounts[select_var,]
colnames(highly_variable_lcpm) <- gsub("_4wpi", "", colnames(highly_variable_lcpm))
dim(highly_variable_lcpm)
head(highly_variable_lcpm)
heatmap.2(highly_variable_lcpm,trace="none",scale="row")

wk10target <- subset(target, weeks_post_infection == 10) # select only wk 10 samples
# View(wk10target) # confirm only Wk. 10 samples are present
wk10counts <- counts[,wk10target$label] # collect only count data for the subsetted samples
colnames(wk10counts)
wk10counts$genes <- counts[,1]
wk10DEcounts <- merge(wk10counts, as.data.frame(wk10DE), by="genes", type="inner")
colnames(wk10DEcounts)
wk10target$treatment <- paste((substr(wk10target$trt_rep, 1, 3)), wk10target$weeks_post_infection, sep="_") # defines the treatment group and date of collection
wk10group <- factor(wk10target$treatment)
wk10group
colnames(wk10DEcounts)
wk10DGE <- DGEList(counts = wk10DEcounts[,2:14], genes = wk10DEcounts[,1], group = wk10group)
logcounts <- cpm(wk10DGE,log=TRUE)
var_genes <- apply(logcounts, 1, var)
head(var_genes)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:162] # looking at all the DE genes
head(select_var)
highly_variable_lcpm <- logcounts[select_var,]
colnames(highly_variable_lcpm) <- gsub("_10wpi", "", colnames(highly_variable_lcpm))
dim(highly_variable_lcpm)
head(highly_variable_lcpm)
heatmap.2(highly_variable_lcpm,trace="none",scale="row")
