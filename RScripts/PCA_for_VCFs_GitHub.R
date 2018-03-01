#####
#
#This script is to process the data produced from the aln.sh file for the "Pseudomonas" project

#Before this point, we have sequenced 49 genomes (two of which were removed from analysis due to contamination)
#The aln.sh pipeline took the raw sequence data and processed it to the point of creating 
#Variant Call Format (VCF) files following the GATK best processes at the time

#This script largely follows the Bioconductor tools recommendations to go from the VCF format to SNP calling

########################################################################
#Modified from "Tutorials for the R/Bioconductor Package SNPRelate
########################################################################
### Install the Bioconductor R packages gdsfmt and SNPRelate from Bioconductor repository ###

#source("http://bioconductor.org/biocLite.R")
#biocLite("gdsfmt")
#biocLite("SNPRelate")

#Install the development version from Github
#library("devtools")
#install_github("zhengxwen/gdsfmt")
#install_github("zhengxwen/SNPRelate")

########################################################################
#Preparing Data
########################################################################

#Load the R packages: gdsfmt and SNPRelate

#The packages may need to be downloaded from bioconductor if the "library()" command doesn't work

## try http:// if https:// URLs are not supported

#install the latest version of bioconductor 
#source("https://bioconductor.org/biocLite.R")


#biocLite("BiocUpgrade")
#biocLite("gdsfmt")
library("gdsfmt")

## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("SNPRelate")
library("SNPRelate")

## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
biocLite("VariantAnnotation")
library("VariantAnnotation")

library("ggplot2")

#Set working directory to uploaded files following re-analysis of the data
setwd("C:/Users/Lauren/Dropbox/Pseudomonas_fromThunderHorse/VCFs/NoMods")

#Use  your own VCF file

# #Only need to do this once to make the psamples.gds file
# #store the name of the VCF file as vcf.fn for ease later
# vcf.fn <- "psamplestrim.vcf"
# 
# 
# vcf = readVcf(vcf.fn, ignore.chr.prefix = "chr")
# 
# #Reformat from VCF to a GDS file (genomic data structure)
# #GDS is used for storing gentic array-oriented data
# #VCF is the Variant Call Format (VCF) 
# 
# #method - "copy.num.of.ref" to extract and store dosage of the reference allele
# snpgdsVCF2GDS(vcf.fn, "psamples.gds", method = "copy.num.of.ref")
# 

#print the info stores in the gds
snpgdsSummary("psamples.gds")

genofile <-snpgdsOpen("psamples.gds")
genofile

#snpgdsGDS2BED(genofile, bed.fn = "test")

#pca from snps
pca = snpgdsPCA(genofile, autosome.only = FALSE)

#Calculate the percent of variance explained by each principal component
pc.percent <- pca$varprop*100
head(round(pc.percent,2))

#store the eigenvalues for each sample in a data frame
poptab <- data.frame(sample.id = pca$sample.id,
                     pop_code = 0,
                  EV1 = pca$eigenvect[,1],
                  EV2 = pca$eigenvect[,2], 
                  stringsAsFactors = FALSE)

 
poptab

#import metadata
meta = read.csv("C:/Users/Lauren/Dropbox/Pseudomonas/Files/Pseudomonas_metadata.csv", header = TRUE)

#fix some formatting errors in the meta data/sample names
meta$sample.id = gsub("\\-.*","",meta$sample.id)

meta$sample.id[8] = "2864_10_012_CTTGTA"
meta$sample.id[1] = "1005_1_001ATCACG"

#combine metadata with sample ids
poptab_meta = merge(poptab, meta, by = "sample.id")


#recreate plot in ggplot
ggplot(poptab_meta)+
  geom_point(aes(x = EV1, y = EV2, color = poptab_meta$Area), size = 4)
  
#compute a dissimilarity matrix from snps using snpgdsDiss
 dissMatrix = snpgdsDiss(genofile, sample.id = NULL, snp.id = NULL, remove.monosnp = TRUE, 
                         autosome.only = FALSE, maf = NaN, missing.rate = NaN, verbose = TRUE)
 head(dissMatrix)
 
 diss = as.dist(dissMatrix$diss)
 head(diss)

#compute a dissimilarity matrix from the snps using the poppr package



library(vegan)

#do adonis test / permanova test on data to test for differences in groups
adonis(Diss ~ Area, meta, perm = 200)
