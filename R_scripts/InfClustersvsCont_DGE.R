#################################################
##                                             ##
## 82-Cattle-16 Project RNA-seq DGE analysis   ##
## TB infected clusters vs uninfected cattle   ##
## Unstimulated whole blood direct collections ##
##                                             ##
#################################################

## Cluster A = less severely infected
## Cluster B = more severely infected

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
library(dplyr) # using version 0.8.0.1
library(tidyverse) # using version 1.2.1
library(VennDiagram) # using version 1.6.20
library(gplots) # using version 3.0.1.1
library(pca3d)  # using version 0.10
library(reshape2)  # using version 1.4.3

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

# Create MDS of infected animals at individual timepoints to see if they cluster out:
### Subset count data into just samples from wk 0:
wk0target <- subset(target, treatment == "inf_0wpi") # select for only Wk. 0 samples
# View(wk0target) # confirm only Wk. 0 samples are present
wk0counts <- counts[,wk0target$label] # collect only count data for the subsetted samples
colnames(wk0counts) <- wk0target$trt_rep
colnames(wk0counts)

# Define treatment factor:
group <- factor(wk0target$treatment)
group

# Create the DGE object:
DGE <- DGEList(counts = wk0counts, genes = counts[,1], group = group)

### Calculate the normalization factors (TMM method):

DGE <- calcNormFactors(DGE)
DGE$samples # shows adjusted normalization factors

### Examine samples in MDS plot:
col <- c("orange", "green4", "green4", "blue", "blue", "green4", "blue", "blue")
plotMDS(DGE, top = 500, col = c(col))

### Subset count data into just samples from wk 4:
wk4target <- subset(target, treatment == "inf_4wpi") # select for only Wk. 4 samples
# View(wk4target) # confirm only Wk. 4 samples are present
wk4counts <- counts[,wk4target$label] # collect only count data for the subsetted samples
colnames(wk4counts) <- wk4target$trt_rep
colnames(wk4counts)

# Define treatment factor:
group <- factor(wk4target$treatment)
group

# Create the DGE object:
DGE <- DGEList(counts = wk4counts, genes = counts[,1], group = group)

### Calculate the normalization factors (TMM method):

DGE <- calcNormFactors(DGE)
DGE$samples # shows adjusted normalization factors

### Examine samples in MDS plot:
col <- c("orange", "green4", "green4", "blue", "blue", "green4", "blue", "blue")
plotMDS(DGE, top = 500, col = c(col))

### Subset count data into just samples from wk 10:
wk10target <- subset(target, treatment == "inf_10wpi") # select for only Wk. 10 samples
# View(wk10target) # confirm only Wk. 10 samples are present
wk10counts <- counts[,wk10target$label] # collect only count data for the subsetted samples
colnames(wk10counts) <- wk10target$trt_rep
colnames(wk10counts)

# Define treatment factor:
group <- factor(wk10target$treatment)
group

# Create the DGE object:
DGE <- DGEList(counts = wk10counts, genes = counts[,1], group = group)

### Calculate the normalization factors (TMM method):

DGE <- calcNormFactors(DGE)
DGE$samples # shows adjusted normalization factors

### Examine samples in MDS plot:
col <- c("orange", "green4", "green4", "blue", "blue", "green4", "blue", "blue")
plotMDS(DGE, top = 500, col = c(col))

#######################################################################################

### Perform DGE analysis of only infected animals by clustering found in heatmaps:
### Subset count data into just samples from infected animals at wk 0, 4, & 10:
# View(target) # confirm only Wk. 0, 4, & 10 samples are present
target = target[-c(16:18),] # remove inf1 animal samples
# View(target)
### Assign cluster treatments to each sample:
# Inf2, 3, & 6 are cluster A, all other Inf are cluster B
c <- rep("C", 15)
a1 <- rep("A", 6)
b1 <- rep("B", 6)
a2 <- rep("A", 3)
b2 <- rep("B", 6)
clusters <- c(c, a1, b1, a2, b2)

target$cluster = clusters
target$cluster <- paste(target$cluster, target$weeks_post_infection, sep="_") # defines the clustered samples by date
target$cluster <- paste(target$cluster,"wpi", sep="")
target$cluster

### Subset count data to remove inf1 samples that don't cluster consistently:
# View(target) # confirm only Wk. 0, 4, & 10 samples are present
subcounts <- counts[,target$label] # collect only count data for the subsetted samples
colnames(subcounts) # should only have wk 0, 4, & 10 infected samples

# Define treatment factor (uninfected vs infected):
cluster <- factor(target$cluster)
cluster

# Create the DGE object:
DGE <- DGEList(counts = subcounts, genes = counts[,1], group = cluster)

# Define the experimental design and generalized linear model fit:
design <- model.matrix(~ 0 + cluster) # defining using group-means parameterization
design

### Calculate the normalization factors (TMM method):

DGE <- calcNormFactors(DGE)
DGE$samples # shows adjusted normalization factors

### Obtain raw gene count densities for quality check:
rawDGE <- log10(DGE$counts[,1:ncol(DGE$counts)] + 1) # log transform the raw data

### Filter out lowly expressed genes:

# Filter out genes with less than 0.5 counts per million (cpm) in at least 5 samples:
keep <- rowSums(cpm(DGE) > 1) >= 3

table(keep) # shows how many genes are kept vs. filtered out

DGE <- DGE[keep, , keep.lib.sizes = FALSE] # removes lowly expressed genes from DGE object

write.table(DGE$genes, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/AvsBvsC_filtered_genes.txt", sep = "\t")

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
col <- c("green", "green", "green", "red", "red", "red", "blue", "blue", "blue")
pch <- c(16, 17, 15, 16, 17, 15, 16, 17, 15)

plotMDS(DGE,
        top = 500,
        col = col [as.factor(cluster)],
        pch = pch [as.factor(cluster)])
legend("bottomright", legend=c(levels(cluster)), cex = .65, pch = pch, col = col)

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

### Week 0 cluster A vs control group:
# Perform likelihood ratio test & identify DEGs:
wk0Alrt <- glmLRT(fit, contrast = c(1, 0, 0, 0, 0, 0, -1, 0, 0)) # perform likelihood ratio test 
summary(wk0ADGE <- decideTestsDGE(wk0Alrt)) # reports upregulated, downregulated, and not DE genes
wk0AtopTags <- topTags(wk0Alrt, n = Inf) # table of genes sorted by FDR value (FDR < .05 means DE)
wk0ADE <- subset(wk0AtopTags$table, FDR < .05) # create a new table of only DE (FDR < .05) genes
wk0AUp <- subset(wk0ADE, logFC > 0)
wk0ADown <- subset(wk0ADE, logFC < 0)
head(wk0ADE, n = 10)
# Create plot smear:
wk0ADEnames <- rownames(wk0ADE)
plotSmear(wk0Alrt, de.tags = wk0ADEnames, cex = .75, ylim = c(-6,6))
abline(h = c(-1, 1), col = "blue") 

### Week 0 cluster B vs control group:
# Perform likelihood ratio test & identify DEGs:
wk0Blrt <- glmLRT(fit, contrast = c(0, 0, 0, 1, 0, 0, -1, 0, 0)) # perform likelihood ratio test 
summary(wk0BDGE <- decideTestsDGE(wk0Blrt)) # reports upregulated, downregulated, and not DE genes
wk0BtopTags <- topTags(wk0Blrt, n = Inf) # table of genes sorted by FDR value (FDR < .05 means DE)
wk0BDE <- subset(wk0BtopTags$table, FDR < .05) # create a new table of only DE (FDR < .05) genes
wk0BUp <- subset(wk0BDE, logFC > 0)
wk0BDown <- subset(wk0BDE, logFC < 0)
head(wk0BDE, n = 10)
# Create plot smear:
wk0BDEnames <- rownames(wk0BDE)
plotSmear(wk0Blrt, de.tags = wk0BDEnames, cex = .75, ylim = c(-6,6))
abline(h = c(-1, 1), col = "blue") 
write.table(wk0BUp$genes, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/InfBvsCont_wk0_up.txt", sep = "\t")
write.table(wk0BDown$genes, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/InfBvsCont_wk0_down.txt", sep = "\t")
write.table(wk0BDE, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/InfBvsCont_wk0_allDE.txt", sep = "\t")

### Week 0 cluster A vs cluster B group:
wk0ABlrt <- glmLRT(fit, contrast = c(1, 0, 0, -1, 0, 0, 0, 0, 0))
summary(wk0ABDGE <- decideTestsDGE(wk0ABlrt)) # reports upregulated, downregulated, and not DE genes
wk0ABtopTags <- topTags(wk0ABlrt, n = Inf) # table of genes sorted by FDR value (FDR < .05 means DE)
wk0ABDE <- subset(wk0ABtopTags$table, FDR < .05) # create a new table of only DE (FDR < .05) genes
wk0ABUp <- subset(wk0ABDE, logFC > 0)
wk0ABDown <- subset(wk0ABDE, logFC < 0)
head(wk0ABDE, n = 10)
# Create plot smear:
wk0ABDEnames <- rownames(wk0ABDE)
plotSmear(wk0ABlrt, de.tags = wk0ABDEnames, cex = .75, ylim = c(-6,6))
abline(h = c(-1, 1), col = "blue") 
write.table(wk0ABUp$genes, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/InfAvsB_wk0_up.txt", sep = "\t")
write.table(wk0ABDown$genes, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/InfAvsB_wk0_down.txt", sep = "\t")
write.table(wk0ABDE, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/InfAvsB_wk0_allDE.txt", sep = "\t")

### Week 4 cluster A vs control group:
wk4Alrt <- glmLRT(fit, contrast = c(0, 0, 1, 0, 0, 0, 0, 0, -1))
summary(wk4ADGE <- decideTestsDGE(wk4Alrt)) # reports upregulated, downregulated, and not DE genes
wk4AtopTags <- topTags(wk4Alrt, n = Inf) # table of genes sorted by FDR value (FDR < .05 means DE)
wk4ADE <- subset(wk4AtopTags$table, FDR < .05) # create a new table of only DE (FDR < .05) genes
wk4AUp <- subset(wk4ADE, logFC > 0)
wk4ADown <- subset(wk4ADE, logFC < 0)
head(wk4ADE, n = 10)
# Create plot smear:
wk4ADEnames <- rownames(wk4ADE)
plotSmear(wk4Alrt, de.tags = wk4ADEnames, cex = .75, ylim = c(-6,6))
abline(h = c(-1, 1), col = "blue") 
write.table(wk4AUp$genes, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/InfAvsCont_wk4_up.txt", sep = "\t")
write.table(wk4ADown$genes, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/InfAvsCont_wk4_down.txt", sep = "\t")
write.table(wk4ADE, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/InfAvsCont_wk4_allDE.txt", sep = "\t")

### Week 4 cluster B vs control group:
wk4Blrt <- glmLRT(fit, contrast = c(0, 0, 0, 0, 0, 1, 0, 0, -1))
summary(wk4BDGE <- decideTestsDGE(wk4Blrt)) # reports upregulated, downregulated, and not DE genes
wk4BtopTags <- topTags(wk4Blrt, n = Inf) # table of genes sorted by FDR value (FDR < .05 means DE)
wk4BDE <- subset(wk4BtopTags$table, FDR < .05) # create a new table of only DE (FDR < .05) genes
wk4BUp <- subset(wk4BDE, logFC > 0)
wk4BDown <- subset(wk4BDE, logFC < 0)
head(wk4BDE, n = 10)
# Create plot smear:
wk4BDEnames <- rownames(wk4BDE)
plotSmear(wk4Blrt, de.tags = wk4BDEnames, cex = .75, ylim = c(-6,6))
abline(h = c(-1, 1), col = "blue") 
write.table(wk4BUp$genes, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/InfBvsCont_wk4_up.txt", sep = "\t")
write.table(wk4BDown$genes, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/InfBvsCont_wk4_down.txt", sep = "\t")
write.table(wk4BDE, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/InfBvsCont_wk4_allDE.txt", sep = "\t")

### Week 4 cluster A vs cluster B group:
wk4ABlrt <- glmLRT(fit, contrast = c(0, 0, 1, 0, 0, -1, 0, 0, 0))
summary(wk4ABDGE <- decideTestsDGE(wk4ABlrt)) # reports upregulated, downregulated, and not DE genes
wk4ABtopTags <- topTags(wk4ABlrt, n = Inf) # table of genes sorted by FDR value (FDR < .05 means DE)
wk4ABDE <- subset(wk4ABtopTags$table, FDR < .05) # create a new table of only DE (FDR < .05) genes
wk4ABUp <- subset(wk4ABDE, logFC > 0)
wk4ABDown <- subset(wk4ABDE, logFC < 0)
head(wk4ABDE, n = 10)
# Create plot smear:
wk4ABDEnames <- rownames(wk4ABDE)
plotSmear(wk4ABlrt, de.tags = wk4ABDEnames, cex = .75, ylim = c(-6,6))
abline(h = c(-1, 1), col = "blue") 
write.table(wk4ABUp$genes, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/InfAvsB_wk4_up.txt", sep = "\t")
write.table(wk4ABDown$genes, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/InfAvsB_wk4_down.txt", sep = "\t")
write.table(wk4ABDE, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/InfAvsB_wk4_allDE.txt", sep = "\t")

### Week 10 cluster A vs control group:
wk10Alrt <- glmLRT(fit, contrast = c(0, 1, 0, 0, 0, 0, 0, -1, 0))
summary(wk10ADGE <- decideTestsDGE(wk10Alrt)) # reports upregulated, downregulated, and not DE genes
wk10AtopTags <- topTags(wk10Alrt, n = Inf) # table of genes sorted by FDR value (FDR < .05 means DE)
wk10ADE <- subset(wk10AtopTags$table, FDR < .05) # create a new table of only DE (FDR < .05) genes
wk10AUp <- subset(wk10ADE, logFC > 0)
wk10ADown <- subset(wk10ADE, logFC < 0)
head(wk10ADE, n = 10)
# Create plot smear:
wk10ADEnames <- rownames(wk10ADE)
plotSmear(wk10Alrt, de.tags = wk10ADEnames, cex = .75, ylim = c(-6,6))
abline(h = c(-1, 1), col = "blue") 

### Week 10 cluster B vs control group:
wk10Blrt <- glmLRT(fit, contrast = c(0, 0, 0, 0, 1, 0, 0, -1, 0))
summary(wk10BDGE <- decideTestsDGE(wk10Blrt)) # reports upregulated, downregulated, and not DE genes
wk10BtopTags <- topTags(wk10Blrt, n = Inf) # table of genes sorted by FDR value (FDR < .05 means DE)
wk10BDE <- subset(wk10BtopTags$table, FDR < .05) # create a new table of only DE (FDR < .05) genes
wk10BUp <- subset(wk10BDE, logFC > 0)
wk10BDown <- subset(wk10BDE, logFC < 0)
head(wk10BDE, n = 10)
# Create plot smear:
wk10BDEnames <- rownames(wk10BDE)
plotSmear(wk10Blrt, de.tags = wk10BDEnames, cex = .75, ylim = c(-6,6))
abline(h = c(-1, 1), col = "blue") 
write.table(wk10BUp$genes, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/InfBvsCont_wk10_up.txt", sep = "\t")
write.table(wk10BDown$genes, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/InfBvsCont_wk10_down.txt", sep = "\t")
write.table(wk10BDE, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/InfBvsCont_wk10_allDE.txt", sep = "\t")

### Week 10 cluster A vs cluster B group:
wk10ABlrt <- glmLRT(fit, contrast = c(0, 1, 0, 0, -1, 0, 0, 0, 0))
summary(wk10ABDGE <- decideTestsDGE(wk10ABlrt)) # reports upregulated, downregulated, and not DE genes
wk10ABtopTags <- topTags(wk10ABlrt, n = Inf) # table of genes sorted by FDR value (FDR < .05 means DE)
wk10ABDE <- subset(wk10ABtopTags$table, FDR < .05) # create a new table of only DE (FDR < .05) genes
wk10ABUp <- subset(wk10ABDE, logFC > 0)
wk10ABDown <- subset(wk10ABDE, logFC < 0)
head(wk10ABDE, n = 10)
# Create plot smear:
wk10ABDEnames <- rownames(wk10ABDE)
plotSmear(wk10ABlrt, de.tags = wk10ABDEnames, cex = .75, ylim = c(-6,6))
abline(h = c(-1, 1), col = "blue") 
write.table(wk10ABUp$genes, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/InfAvsB_wk10_up.txt", sep = "\t")
write.table(wk10ABDown$genes, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/InfAvsB_wk10_down.txt", sep = "\t")
write.table(wk10ABDE, "/Users/jayne.wiarda/Desktop/TB/Thesis_Project/82_Cattle_RNAseq/geneLists/InfAvsB_wk10_allDE.txt", sep = "\t")

### Venn Diagram comparison of DE genes:
vd <- venn.diagram(x = list("wk10A" = wk10AUp$genes, "wk4A" = wk4AUp$genes, "wk10B" = wk10BUp$genes, "wk4B" = wk4BUp$genes), filename = NULL, fill = c("yellow", "orange", "green", "blue"))
grid.draw(vd, recording = TRUE)
vd <- venn.diagram(x = list("wk10A" = wk10ADown$genes, "wk4A" = wk4ADown$genes, "wk10B" = wk10BDown$genes, "wk4B" = wk4BDown$genes), filename = NULL, fill = c("yellow", "orange", "green", "blue"))
grid.draw(vd, recording = TRUE)
