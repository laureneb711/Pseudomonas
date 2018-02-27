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


## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("VariantAnnotation")
library("VariantAnnotation")
#gds file
#example <- snpgdsOpen(snpgdsExampleFileName())


setwd("C:/Users/Lauren/Dropbox/Pseudomonas/NoMod")

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

meta = read.csv("C:/Users/Lauren/Dropbox/Pseudomonas/Pseudomonas_metadata.csv", header = TRUE)

poptab_meta = merge(poptab, meta, by = "sample.id")



ggplot(poptab_meta)+
  geom_point(aes(x = EV1, y = EV2, color = poptab_meta$), size = 4)#+
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

ggdendrogram(cutTree$dendrogram)

newdend = dendro_data(cutTree$dendrogram)
str(newdend)

library(ggplot2)

ggplot(segment(newdend))+
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend))+
  geom_text(data=newdend$labels, aes(x=x, y=y, label=label))+
  coord_flip()+scale_y_reverse(expand=c(0.2,0))+
  theme_dendro()


diss = as.dist(dissMatrix$diss)

head(meta)

library(vegan)
adonis(diss ~ Area, meta, perm = 200)
