#!/bin/bash

##############################################################
##	       		       				    ##
## 82-Cattle-16 Project RNA-seq data bioinformatic pipeline ##
## Tuberculosis Infected vs. Uninfected	/ Time Course       ##
## Unstimulated whole blood direct collections 		    ##
## Illumina HiSeq, Paired-end reads, 2 x 100		    ## 
##							    ##
##############################################################	

## Run on Unix shell

## Samples collected pre-infection & post-infection (4 weeks & 10 weeks post-infection)

## Samples (n = 39) sequenced on Illumina HiSeq instrument

## All files downloaded from bioinformatics drive, untarred, and transferred into 
## a new directory folder (/Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/Fastq_raw)

#######################################################
# Do initial quality check of raw reads using FastQC: #
#######################################################

### FastQC v0.11.5 software is required; consult manual at:
### http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create directory for raw read FastQC outputs:

mkdir /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/FastQC_raw

cd !$

### Run FastQC on all raw fastq files:

for file in /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/Fastq_raw/*.fastq.gz
do
	fastqc \
	-f fastq ${file} \
	-o /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/FastQC_raw \
	--noextract \
	--nogroup \
	-t 2 
	done

# Look at FastQC results:

open /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/FastQC_raw/*.html 
		# notice overrepresented sequences of non-discriminatory bases ('NNNNNNNNN') 
		# and TruSeq adapters; files should be trimmed for quality control

###################################################################################
# Trim files with Trimmomatic to remove adapter sequences and poor quality reads: #
###################################################################################

### Trimmomatic 0.36 software is required; consult manual at:
### http://www.usadellab.org/cms/index.php?page=trimmomatic

# Make a list of file name prefixes:
cd /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/Fastq_raw

rm /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/a # make sure that there is not a file a that already exists	

rm /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/samp_name # make sure that there is not a file samp_name already existing

for i in *.fastq.gz
do  
echo $i; echo "$i" | cut -d'_' -f 1,2 >> /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/a   
done
	
uniq /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/a  >> /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/samp_name 
	
rm /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/a

cat /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/samp_name

# Create directory for Trimmomatic outputs:

mkdir /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/Fastq_trimmed

### Run Trimmomatic on all raw fastq files:

cd /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/Fastq_raw

for i in $(cat /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/samp_name) 

do
java -jar /Applications/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
/Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/Fastq_raw/${i}_R1.fastq.gz \
/Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/Fastq_raw/${i}_R2.fastq.gz \
/Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/Fastq_trimmed/${i}_R1_forward_paired.fastq.gz \
/Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/Fastq_trimmed/${i}_R1_forward_unpaired.fastq.gz \
/Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/Fastq_trimmed/${i}_R2_reverse_paired.fastq.gz \
/Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/Fastq_trimmed/${i}_R2_reverse_unpaired.fastq.gz \
ILLUMINACLIP:/Applications/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:36
done

###################################################
# Do quality check of trimmed reads using FastQC: #
###################################################

# Create directory for trimmed read FastQC outputs:

mkdir /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/FastQC_trimmed 

### Run FastQC on all trimmed files:

cd /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/Fastq_trimmed

for i in $(cat /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/samp_name)
do
fastqc \
-o /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/FastQC_trimmed \
--noextract \
--nogroup \
-t 2 \
/Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/Fastq_trimmed/${i}_R1_forward_paired.fastq.gz
done

for i in $(cat /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/samp_name)
do
fastqc \
-o /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/FastQC_trimmed \
--noextract \
--nogroup \
-t 2 \
/Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/Fastq_trimmed/${i}_R2_reverse_paired.fastq.gz
done

open /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/FastQC_trimmed/*.html 
# overrepresented sequences should be gone and files should have good quality reads

#############################################################
# Align trimmed fastq files to Bos taurus genome with STAR: #
#############################################################

### STAR 2.5.3a software is required; consult manual at:
### https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

### Download and store the NCBI genome FASTA file and gff annotation file for 
### Bos_taurus_UMD_3.1.1 assembly (RefSeq accession GCF_00003055.6):

# Download genome file from NCBI
wget -P /Users/jayne.wiarda/RNASeq/Bos_taurus_UMD3.1.1/ \
"ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/055/GCF_000003055.6_Bos_taurus_UMD_3.1.1/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna.gz"

# Download annotation file from NCBI:
wget -P /Users/jayne.wiarda/RNASeq/ \
"ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/055/GCF_000003055.6_Bos_taurus_UMD_3.1.1/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff.gz"

# Uncompress zipped files:
gunzip /Users/jayne.wiarda/RNASeq/Bos_taurus_UMD3.1.1/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna.gz
gunzip /Users/jayne.wiarda/RNASeq/Bos_taurus_UMD3.1.1/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff.gz

### Generate a genomic index using STAR:

# Create a directory to store the index in:
cd /Users/jayne.wiarda/RNASeq
mkdir genomeDir_UMD3.1.1

# Generate the genomic index:
STAR \
--runThreadN 20 \
--runMode genomeGenerate \
--genomeDir /Users/jayne.wiarda/RNASeq/genomeDir_UMD3.1.1 \
--genomeFastaFiles /Users/jayne.wiarda/RNASeq/Bos_taurus_UMD3.1.1/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna \
--sjdbGTFfile /Users/jayne.wiarda/RNASeq/Bos_taurus_UMD3.1.1/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbOverhang 99

### Map and align all trimmed fastq files:

cd /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/Fastq_trimmed

for i in $(cat /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/samp_name) 
do
    fwd=${i}_R1_forward_paired.fastq.gz
    rev=${i}_R2_reverse_paired.fastq.gz
        STAR \
        --genomeDir /Users/jayne.wiarda/RNASeq/genomeDir_UMD3.1.1 \
        --readFilesIn $fwd $rev \
        --readFilesCommand gunzip -c \
        --runMode alignReads \
        --runThreadN 20 \
        --outFilterMultimapNmax 10 \
        --outReadsUnmapped Fastx \
        --outFilterMismatchNmax 10 \
        --outSAMtype BAM Unsorted \
        --outFileNamePrefix "$i"_
done

# Move all STAR output files to a separate directory:

mkdir /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/STAR_alignments

cd /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/Fastq_trimmed

mv *{Aligned.out.bam,Aligned.out.sam,Log.final.out,Log.out,Log.progress.out,SJ.out.tab,Unmapped.out.mate1,Unmapped.out.mate2} \
/Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/STAR_alignments

# Delete files with fastq prefixes:

rm /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/a		
		
rm /Users/jayne.wiarda/RNASeq/82_Cattle_RNAseq/samp_name 

### Continue on to R to count reads and conduct DGE analysis!





