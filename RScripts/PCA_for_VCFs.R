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
library("gdsfmt", lib.loc="~/R/win-library/3.2")
library("SNPRelate", lib.loc="~/R/win-library/3.2")

#gds file
#example <- snpgdsOpen(snpgdsExampleFileName())

setwd("C:/Users/Lauren/Dropbox/Pseudomonas/NoMod")

#Use  your own VCF file
vcf.fn <- "psamplestrim.vcf"

#Reformat
snpgdsVCF2GDS(vcf.fn, "psamples.gds", method="biallelic.only")

snpgdsSummary("psamples.gds")

genofile <-snpgdsOpen("psamples.gds")

genofile


#example
genotype = openfn.gds("psamples.gds", readonly=TRUE, allow.duplicate = TRUE)
#exampe = openfn.gds(snpgdsExampleFileName(), readonly=TRUE, allow.duplicate = TRUE)


#example <- snpgdsPCA(exampe)

pca <- snpgdsPCA(genotype, autosome.only=FALSE)

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
  geom_point(aes(x = EV1, y = EV2), size = 4)+
  geom_text(aes(x = EV1+.05, y = EV2+.05, label = sample.id))
  

##################
#Kmeans from PCA#
##################
pca

#use elbow method to determine best vectors to determine variance
plot(pca$varprop)

pca$loadings = pca$eigenvect*sqrt(pca$eigenval)

ggplot()+
  geom_point(aes(x = pca$loadings[,1], y = pca$loadings[,2]), size = 8, shape = 5)

k4 = kmeans(pca$loadings, 4)

plot(comp, col=k4$clust, pch=16)




out=princomp(matrix(c(1,1,-1,-1),byrow=T,nrow=2,ncol=2))
out$score



#############
#Install pophelper package

