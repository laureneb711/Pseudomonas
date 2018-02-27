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
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade")
#biocLite("gdsfmt")
library("gdsfmt")

## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("SNPRelate")
library("SNPRelate")


## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("VariantAnnotation")
library("VariantAnnotation")

#Set working directory to uploaded files following re-analysis of the data
setwd("C:/Users/Lauren/Dropbox/Pseudomonas_fromThunderHorse/VCFs/NoMods")

#Use  your own VCF file
vcf.fn <- "psamplestrim.vcf"

vcf = readVcf("psamplestrim.vcf")

#Reformat to GDS file
snpgdsVCF2GDS("psamplestrim.vcf", "psamples.gds", method = "copy.num.of.ref")

snpgdsSummary("psamples.gds")

genofile <-snpgdsOpen("psamples.gds", allow.duplicate = T)

genofile

#store various info
sample.id = read.gdsn(index.gdsn(genofile, "sample.id"))

pop_code = read.gdsn(index.gdsn(genofile, "sample.id"))

######
#if there is population information

sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

library(ggplot2)

meta = read.csv("C:/Users/Lauren/Dropbox/Pseudomonas/Files/Pseudomonas_metadata.csv", header = TRUE)

poptab_meta = merge(poptab, meta, by = "sample.id")

#dendogram
dissMatrix = snpgdsDiss(genofile, sample.id = NULL, snp.id = NULL, remove.monosnp = TRUE, 
                        autosome.only = FALSE, maf = NaN, missing.rate = NaN, verbose = TRUE)

head(dissMatrix)

snpHCluster = snpgdsHCluster(dissMatrix, sample.id = NULL, need.mat = TRUE, hang = 0.25)
head(snpHCluster)

cutTree = snpgdsCutTree(snpHCluster, z.threshold = 15, outlier.n = 5, 
                        n.perm = 5000, samp.group = NULL, col.outlier = "red", col.list = NULL, 
                        pch.outlier = 4, pch.list = NULL, 
                        label.H = FALSE, label.Z = TRUE, verbose = TRUE)



plot(cutTree$dendrogram)

#snpgdsDrawTree(cutTree)

#make the plot look better

library(ggdendro)

ggdendrogram(cutTree$dendrogram, rotate = T)

#dendr = cutTree$dendrogram

#plot(dendr)

newdend = dendro_data(cutTree$dendrogram)
str(newdend)

#ggdendrogram(newdend, rotate = T)
library(cluster)

clust = cutTree
dendr = dendro_data(clust$dendrogram, type = "rectangle")

clust.df = data.frame(label = poptab_meta$sample.id, cluster = factor(clust))

dendr$labels$label

id = data.frame(sample.id = dendr$labels$label)

id$sample.id
id$Area = NULL

library(plyr)

tail(meta)
tail(id)
head(meta_label)

meta$sample.id


meta$sample.id = gsub("\\-.*","",meta$sample.id)

meta$sample.id[8] = "2864_10_012_CTTGTA"
meta$sample.id[1] = "1005_1_001ATCACG"

  meta_label = join(id, meta, by = "sample.id")
tail(meta)
tail(id)
head(meta_label)

library(dplyr)
df <- mutate_all(meta_label, funs(toupper))

df <- mutate_all(meta_label, funs(toupper))

dendr$leaf_labels = df$Area

dendr$leaf_labels = gsub( " ", "", dendr$leaf_labels) 


ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(dendr), aes(x, y, label=label, hjust=0, color= dendr$leaf_labels), size=5)+
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        legend.text=element_text(size=16),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())+
  scale_color_discrete(name = "Area")

ggsave("C:/Users/Lauren/Dropbox/Pseudomonas/Draft_tree.png", width = 25, height = 13, units = "in")

