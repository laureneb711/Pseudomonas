d#####
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
#source("https://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade")
#biocLite("gdsfmt")
library("gdsfmt")

## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("SNPRelate")
library("SNPRelate")

#Set working directory to uploaded files following re-analysis of the data
setwd("C:/Users/Lauren/Dropbox/Pseudomonas_fromThunderHorse/VCFs/NoMods")

#Use  your own VCF file
vcf.fn <- "psamples_trim.vcf"

vcf = readVcf("psamplestrim.vcf")

#Reformat
snpgdsVCF2GDS(vcf.fn, "psamples.gds", method = "copy.num.of.ref")

snpgdsSummary("psamples.gds")

genofile <-snpgdsOpen("psamples.gds")

genofile


sample.id = read.gdsn(index.gdsn(genofile, "sample.id"))

pop_code = read.gdsn(index.gdsn(genofile, "sample.id"))

pca = snpgdsPCA(genofile, autosome.only = FALSE)

# tab <- data.frame(sample.id = pca$sample.id,pop =
#   factor(pop_code)[match(pca$sample.id, sample.id)],EV1 =
#   pca$eigenvect[,1],EV2 = pca$eigenvect[,2],stringsAsFactors = FALSE)
# 
# plot(tab$EV2, tab$EV1, col=as.integer(tab$pop),xlab="eigenvector 2",
#      ylab="eigenvector 1") 
# 
# 
# legend("topleft", legend=levels(tab$pop),
#                                   pch="o", col=1:nlevels(tab$pop))

pc.percent <- pca$varprop*100
head(round(pc.percent,2))

tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],
                  EV2 = pca$eigenvect[,2], 
                  stringsAsFactors = FALSE)

head(tab)

plot(tab$EV2, tab$EV1, xlab = "eigenvector 2", ylab = "eigenvector 1")

######
#if there is population information

sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# Get population information
#   or pop_code <- scan("pop.txt", what=character())
#   if it is stored in a text file "pop.txt"

poptab <- data.frame(sample.id = pca$sample.id,
                     pop_code = 0,
                  EV1 = pca$eigenvect[,1],
                  EV2 = pca$eigenvect[,2], 
                  stringsAsFactors = FALSE)

  for (i in 1:nrow(poptab)){
    print(i)
    if (grepl("bath", poptab$sample.id[i], ignore.case = TRUE) == TRUE) {
      poptab$pop_code[i] = "Bath"
  #      print("TRUE")
    }
    
    if (grepl("CF", poptab$sample.id[i], ignore.case = TRUE) == TRUE) {
      poptab$pop_code[i] = "CF"
 #     print("TRUE")
    }
    
    if (grepl("Kitchen", poptab$sample.id[i], ignore.case = TRUE) == TRUE) {
      poptab$pop_code[i] = "Kitchen"
#      print("TRUE")
    }
     
  }
    
poptab


library(ggplot2)

meta = read.csv("C:/Users/Lauren/Dropbox/Pseudomonas/Files/Pseudomonas_metadata.csv", header = TRUE)

poptab_meta = merge(poptab, meta, by = "sample.id")



ggplot(poptab_meta)+
  geom_point(aes(x = EV1, y = EV2, color = poptab_meta$Area), size = 4)#+
#  geom_text(aes(x = EV1+.05, y = EV2+.05, label = sample.id))
  

#parallel coordinates plot for the top principal components
library(MASS)

datpop = factor(pop_code)[match(pca$sample.id, sample.id)]
parcoord(pca$eigenvect, col=datpop)


#Calculate SNP correlations betwen eigenvetors and SNP genotypes
chr = read.gdsn(index.gdsn(genofile, "snp.chromosome"))

CORR = snpgdsPCACorr(pca, genofile, eig.which = 1:4)

savepar = par(mfrow=c(3,1), mai=c(0.3, ))

ibs = snpgdsIBS(genofile, autosome.only = FALSE)

pop.idx = order(poptab_meta$Strain.type)

image(ibs$ibs[pop.idx,pop.idx], col = terrain.colors(16))


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

snpgdsDrawTree(cutTree)

library(ggdendro)

ggdendrogram(cutTree$dendrogram, rotate = T)

dendr = cutTree$dendrogram

plot(dendr)
gg<-rect.hclust(dendr,k=0)

newdend = dendro_data(cutTree$dendrogram)
str(newdend)

ggdendrogram(newdend, rotate = T)
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

getwd()
#library(ggplot2)

ggplot(segment(newdend))+
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend))+
  geom_text(data=newdend$labels, aes(x=x, y=y, label=label))+
  coord_flip()+scale_y_reverse(expand=c(0.2,0))+
  theme_dendro()


diss = as.dist(dissMatrix$diss)

head(meta)

library(vegan)
adonis(diss ~ Area, meta, perm = 200)
