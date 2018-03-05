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
#biocLite("gdsfmt, "C:/Users/Lauren/R/win-library/3.4" ")
library("gdsfmt")

## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("SNPRelate", "C:/Users/Lauren/R/win-library/3.4" )
library("SNPRelate")

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
#biocLite("VariantAnnotation", lib = "C:/Users/Lauren/R/win-library/3.4" )
#library(VariantAnnotation, lib.loc = "C:/Users/Lauren/R/win-library/3.4")

library(vegan)

#Set working directory to uploaded files following re-analysis of the data
setwd("C:/Users/Lauren/Dropbox/Pseudomonas_fromThunderHorse/VCFs/FromCluster")

#Use  your own VCF file

 # #Only need to do this once to make the psamples.gds file
 # #store the name of the VCF file as vcf.fn for ease later
 # vcf.fn <- "NoMods_psamplestrim.vcf"
 # 
 # 
 # vcf = readVcf(vcf.fn)
 # 
 # #Reformat from VCF to a GDS file (genomic data structure)
 # #GDS is used for storing gentic array-oriented data
 # #VCF is the Variant Call Format (VCF) 
 # 
 # #method - "copy.num.of.ref" to extract and store dosage of the reference allele
 # snpgdsVCF2GDS(vcf.fn, "NoMods_psamples.gds", method = "copy.num.of.ref")
 # 

#print the info stores in the gds
snpgdsSummary("NoMods_psamples.gds")

genofile <-snpgdsOpen("NoMods_psamples.gds")
genofile

#import metadata
meta = read.csv("C:/Users/Lauren/Dropbox/Pseudomonas/Files/Pseudomonas_metadata.csv", header = TRUE)

meta = subset(meta, meta$Area != "Test")

meta$sample.id
#fix some formatting errors in the meta data/sample names
meta$sample.id = gsub("\\-.*","",meta$sample.id)

meta$sample.id[8] = "2864_10_012_CTTGTA"
meta$sample.id[1] = "1005_1_001ATCACG"


#compute a dissimilarity matrix from snps using snpgdsDiss
 dissMatrix = snpgdsDiss(genofile, sample.id = NULL, snp.id = NULL, remove.monosnp = TRUE, 
                         autosome.only = FALSE, maf = NaN, missing.rate = NaN, verbose = TRUE)
 head(dissMatrix)
 
 str(dissMatrix)
 dist= as.dist(dissMatrix$diss)
 
 dm = as.matrix(dissMatrix$diss)
 
  ############
 #looking for difference without any specific hypothesis
 ad_nohyp = adonis(dist ~ Strain.type, meta, perm = 200, method = "binary")
 
 NMDS=metaMDS(dm,k=2,trymax=10000, distance = "binary", wascores = T)
 
 stressplot(NMDS)
 plot(NMDS)
 
 #plot in ggplot
 library(ggplot2)
 
 NMDS
 
 #adapt from this tutorial to make a nice looking plot
 # ============================================================
 # Tutorial on drawing an NMDS plot using ggplot2
 # by Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz)
 # =============================================================
 
 #Get MDS stats
 
 #Make a new data frame, and put country, latrine, and depth information there, to be useful for coloring, and shape of points
 NMDS_df=data.frame(x=NMDS$point[,1],y=NMDS$point[,2],
                    Strain.type=as.factor(meta$Strain.type),
                    Area=as.factor(meta$Area),
                    House = as.factor(meta$House)
 )
 
 #Get spread of points based on countries
 plot.new()
 ord<-ordiellipse(NMDS, as.factor(meta$Strain.type) ,display = "sites", kind ="sd", conf = 0.95, label = T)
 dev.off()
 
 #Reference: http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
 #Data frame df_ell contains values to show ellipses. It is calculated with function veganCovEllipse which is hidden in vegan package. This function is applied to each level of NMDS (group) and it uses also function cov.wt to calculate covariance matrix.
 veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
 {
   theta <- (0:npoints) * 2 * pi/npoints
   Circle <- cbind(cos(theta), sin(theta))
   t(center + scale * t(Circle %*% chol(cov)))
 }
 
 #Generate ellipse points
 df_ell <- data.frame()
 
 for(g in levels(NMDS_df$Strain.type)){
   if(g!="" && (g %in% names(ord))){
     
     df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS_df[NMDS$Strain.type==g,],
                                                      veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                   ,Strain.type=g))
   }
 }
 
 head(df_ell)
 
 #Generate mean values from NMDS plot grouped on Countries
 NMDS.mean=aggregate(NMDS_df[,1:2],list(group=NMDS_df$Strain.type),mean)
 
 NMDS.mean
 
 #Now do the actual plotting
 
 shape_values<-seq(1,16)
 
 p = ggplot(data=NMDS_df,aes(x,y,colour=Strain.type))+
   annotate("text",x=NMDS.mean$x,y=NMDS.mean$y,label=NMDS.mean$group,size=4)+
   geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2), size=1, linetype=2)+
   geom_point(aes(shape=House))+scale_shape_manual(values=shape_values)+theme_bw() 
 
 p
 #pdf("NMDS.pdf")
 
 ad_nohyp
 
 ############
 #Repeat the process, but this time address hyp 1
 #hyp 1: Comparative genomics of CF patient strains vs household drain strains; 
 #do these sequences cluster as two major clades in the analysis?
 ############
 
 #first, remove 3 healthy strains from the analysis
 meta_hyp1 = subset(meta, meta$Strain.type != "Non-CF UR")
 
  # #Only need to do this once to make the psamples.gds file
  # #store the name of the VCF file as vcf.fn for ease later
  # vcf.fn <- "Hyp1_psamplestrim.vcf"
  # 
  # 
  # vcf = readVcf(vcf.fn)
  # 
  # #Reformat from VCF to a GDS file (genomic data structure)
  # #GDS is used for storing gentic array-oriented data
  # #VCF is the Variant Call Format (VCF) 
  # 
  # #method - "copy.num.of.ref" to extract and store dosage of the reference allele
  # snpgdsVCF2GDS(vcf.fn, "Hyp1_psamples.gds", method = "copy.num.of.ref")
  # 
 #print the info stores in the gds
 snpgdsSummary("Hyp1_psamples.gds")
 
 Hyp1_genofile <-snpgdsOpen("Hyp1_psamples.gds")
 Hyp1_genofile
 
  #compute a dissimilarity matrix from snps using snpgdsDiss
 Hyp1_dissMatrix = snpgdsDiss(Hyp1_genofile, sample.id = NULL, snp.id = NULL, remove.monosnp = TRUE, 
                         autosome.only = FALSE, maf = NaN, missing.rate = NaN, verbose = TRUE)
 head(Hyp1_dissMatrix)
 
 str(Hyp1_dissMatrix)
 dist_hyp1= as.dist(Hyp1_dissMatrix$diss)
 
 dm_hyp1 = as.matrix(Hyp1_dissMatrix$diss)
 
 ad_hyp1 = adonis(dist_hyp1 ~ Strain.type, meta_hyp1, perm = 200, method = "binary")
 
 NMDS_hyp1=metaMDS(dm_hyp1,k=2,trymax=10000, distance = "binary", wascores = T)
 
 stressplot(NMDS_hyp1)
 plot(NMDS_hyp1)
 
 #plot in ggplot
 library(ggplot2)
 
 NMDS_hyp1
 
 #adapt from this tutorial to make a nice looking plot
 # ============================================================
 # Tutorial on drawing an NMDS plot using ggplot2
 # by Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz)
 # =============================================================
 
 #Get MDS stats
 
 #Make a new data frame, and put country, latrine, and depth information there, to be useful for coloring, and shape of points
 NMDS_hyp1_df=data.frame(x=NMDS_hyp1$point[,1],y=NMDS_hyp1$point[,2],
                         Strain.type=as.factor(meta_hyp1$Strain.type),
                         Area=as.factor(meta_hyp1$Area),
                         House = as.factor(meta_hyp1$House)
 )
 
 plot.new()
 ord<-ordiellipse(NMDS_hyp1, as.factor(meta_hyp1$Strain.type) ,display = "sites", kind ="sd", conf = 0.95, label = T)
 dev.off()
 
 #Reference: http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
 #Data frame df_ell contains values to show ellipses. It is calculated with function veganCovEllipse which is hidden in vegan package. This function is applied to each level of NMDS (group) and it uses also function cov.wt to calculate covariance matrix.
 veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
 {
   theta <- (0:npoints) * 2 * pi/npoints
   Circle <- cbind(cos(theta), sin(theta))
   t(center + scale * t(Circle %*% chol(cov)))
 }
 
 #Generate ellipse points
 df_ell <- data.frame()
 
 for(g in levels(NMDS_hyp1_df$Strain.type)){
   if(g!="" && (g %in% names(ord))){
     
     df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS_hyp1_df[NMDS_hyp1$Strain.type==g,],
                                                      veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                   ,Strain.type=g))
   }
 }
 
 head(df_ell)
 
 #Generate mean values from NMDS plot grouped on strain type
 NMDS.mean=aggregate(NMDS_hyp1_df[,1:2],list(group=NMDS_hyp1_df$Strain.type),mean)
 
 NMDS.mean
 
 #Now do the actual plotting
 
 shape_values<-seq(1,16)
 
 p = ggplot(data=NMDS_hyp1_df,aes(x,y,colour=Strain.type))+
   annotate("text",x=NMDS.mean$x,y=NMDS.mean$y,label=NMDS.mean$group,size=4)+
   geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2), size=1, linetype=2)+
   geom_point(aes(shape=House))+scale_shape_manual(values=shape_values)+theme_bw() 
 
 p
 
 ggsave("C:/Users/Lauren/Dropbox/Pseudomonas/NMDS_hyp1_SNP.png")
 #pdf("NMDS.pdf")
 
 ad_hyp1
 
 
 ############
 #Repeat the process, but this time address hyp 2
 #hyp 2:Comparative genomics of strains isolated from within CF patients vs location(s) within their houses; 
 #do these pairs cluster according to household environment? 
 #(need sufficient power in this analysis, and target is to have on order of 10 or more pairs for the analysis)
 ############
 #Note: there are 4 homes that fit this analysis.  
 #For 2 homes have bath drain /kitchen drain /human strains.
 #For 2 homes we have 1 bath drain & 1 sputum sample. 
 #The 5th house that had a individual with CF living there, did not have any drain strains. 
 #I don’t know if this number of strains is sufficient for this analysis.
 
 #first, remove 3 healthy strains from the analysis
 table(meta$House)
 
 meta_hyp2 = meta[meta$House %in% names(which(table(meta$House) > 2)), ]
 meta_hyp2 = subset(meta_hyp2, meta_hyp2$House != ".")
 meta_hyp2 = subset(meta_hyp2, meta_hyp2$House != "")
 
 # #Only need to do this once to make the psamples.gds file
 # #store the name of the VCF file as vcf.fn for ease later
# vcf.fn <- "Hyp2_psamplestrim.vcf"
 # 
 #  vcf = readVcf(vcf.fn)
 # 
 # #Reformat from VCF to a GDS file (genomic data structure)
 # #GDS is used for storing gentic array-oriented data
 # #VCF is the Variant Call Format (VCF) 
 # 
 # #method - "copy.num.of.ref" to extract and store dosage of the reference allele
 #  snpgdsVCF2GDS(vcf.fn, "Hyp2_psamples.gds", method = "copy.num.of.ref")
 # 
 #print the info stores in the gds
 snpgdsSummary("Hyp2_psamples.gds")
 
 Hyp2_genofile <-snpgdsOpen("Hyp2_psamples.gds")
 Hyp2_genofile
 
 #compute a dissimilarity matrix from snps using snpgdsDiss
 Hyp2_dissMatrix = snpgdsDiss(Hyp2_genofile, sample.id = NULL, snp.id = NULL, remove.monosnp = TRUE, 
                              autosome.only = FALSE, maf = NaN, missing.rate = NaN, verbose = TRUE)
 head(Hyp2_dissMatrix)
 
 str(Hyp2_dissMatrix)
 dist_hyp2= as.dist(Hyp2_dissMatrix$diss)
 
 dm_hyp2 = as.matrix(Hyp2_dissMatrix$diss) 

 ad_hyp2 = adonis(dist_hyp2 ~ House, meta_hyp2, perm = 200, method = "binary")
 
 ad_hyp2
 
 NMDS_hyp2=metaMDS(dm_hyp2,k=2,trymax=1000, distance = "binary", wascores = T)
 
 stressplot(NMDS_hyp2)
 plot(NMDS_hyp2)
 
 #plot in ggplot
 library(ggplot2)
 
 NMDS_hyp2
 
 #adapt from this tutorial to make a nice looking plot
 # ============================================================
 # Tutorial on drawing an NMDS plot using ggplot2
 # by Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz)
 # =============================================================
 
 #Get MDS stats
 
 #Make a new data frame, and put country, latrine, and depth information there, to be useful for coloring, and shape of points
 NMDS_hyp2_df=data.frame(x=NMDS_hyp2$point[,1],y=NMDS_hyp2$point[,2],
                         Strain.type=as.factor(meta_hyp2$Strain.type),
                         Area=as.factor(meta_hyp2$Area),
                         House = as.factor(meta_hyp2$House)
 )
 
 plot.new()
 ord<-ordiellipse(NMDS_hyp2, as.factor(meta_hyp2$House) ,display = "sites", kind ="sd", conf = 0.95, label = T)
 dev.off()
 #Reference: http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
 #Data frame df_ell contains values to show ellipses. It is calculated with function veganCovEllipse which is hidden in vegan package. This function is applied to each level of NMDS (group) and it uses also function cov.wt to calculate covariance matrix.
 veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
 {
   theta <- (0:npoints) * 2 * pi/npoints
   Circle <- cbind(cos(theta), sin(theta))
   t(center + scale * t(Circle %*% chol(cov)))
 }
 
 #Generate ellipse points
 df_ell <- data.frame()
 
 for(g in levels(NMDS_hyp2_df$House)){
   if(g!="" && (g %in% names(ord))){
     
     df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS_hyp2_df[NMDS_hyp2$House==g,],
                                                      veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                   ,House=g))
   }
 }
 
 head(df_ell)
 
 #Generate mean values from NMDS plot grouped on strain type
 NMDS.mean=aggregate(NMDS_hyp2_df[,1:2],list(group=NMDS_hyp2_df$House),mean)
 
 NMDS.mean
 
 #Now do the actual plotting
 
 shape_values<-seq(1,16)
 
 p = ggplot(data=NMDS_hyp2_df,aes(x,y,colour=House))+
   annotate("text",x=NMDS.mean$x,y=NMDS.mean$y,label=NMDS.mean$group,size=4)+
   geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2), size=1, linetype=2)+
   geom_point(aes(shape=Strain.type))+scale_shape_manual(values=shape_values)+theme_bw() 
 
 p
 
 ggsave("C:/Users/Lauren/Dropbox/Pseudomonas/NMDS_hyp2_SNP.png")
 
 ad_hyp2
 
 ############
 #Repeat the process, but this time address hyp 3
 #hyp 3: Comparative genomics of household strains in light of 3 kinds of drains from which they were isolated 
 #(kitchen, bathroom, sink); 
 #is there an association between sampling location and pattern of relatedness?
 ############
 
 #first, remove 3 healthy strains from the analysis
 meta_hyp3 = subset(meta, meta$Strain.type == "drain")
 
 # #Only need to do this once to make the psamples.gds file
 # #store the name of the VCF file as vcf.fn for ease later
 # vcf.fn <- "Hyp3_psamplestrim.vcf"
 # 
 # vcf = readVcf(vcf.fn)
 # 
 # #Reformat from VCF to a GDS file (genomic data structure)
 # #GDS is used for storing gentic array-oriented data
 # #VCF is the Variant Call Format (VCF) 
 # 
 # #method - "copy.num.of.ref" to extract and store dosage of the reference allele
 # snpgdsVCF2GDS(vcf.fn, "Hyp3_psamples.gds", method = "copy.num.of.ref")
 # 
 #print the info stores in the gds
 snpgdsSummary("Hyp3_psamples.gds")
 
 Hyp3_genofile <-snpgdsOpen("Hyp3_psamples.gds")
 Hyp3_genofile
 
 #compute a dissimilarity matrix from snps using snpgdsDiss
 Hyp3_dissMatrix = snpgdsDiss(Hyp3_genofile, sample.id = NULL, snp.id = NULL, remove.monosnp = TRUE, 
                              autosome.only = FALSE, maf = NaN, missing.rate = NaN, verbose = TRUE)
 head(Hyp3_dissMatrix)
 
 str(Hyp3_dissMatrix)
 dist_hyp3= as.dist(Hyp3_dissMatrix$diss)
 
 dm_hyp3 = as.matrix(Hyp3_dissMatrix$diss)  
 ad_hyp3 = adonis(dist_hyp3 ~ Area, meta_hyp3, perm = 200, method = "binary")
 
 ad_hyp3
 
 NMDS_hyp3=metaMDS(dm_hyp3,k=2,trymax=1000, distance = "binary", wascores = T)
 
 stressplot(NMDS_hyp3)
 plot(NMDS_hyp3)
 
 #plot in ggplot
 library(ggplot2)
 
 NMDS_hyp3
 
 #adapt from this tutorial to make a nice looking plot
 # ============================================================
 # Tutorial on drawing an NMDS plot using ggplot2
 # by Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz)
 # =============================================================
 
 #Get MDS stats
 
 #Make a new data frame, and put country, latrine, and depth information there, to be useful for coloring, and shape of points
 NMDS_hyp3_df=data.frame(x=NMDS_hyp3$point[,1],y=NMDS_hyp3$point[,2],
                         Strain.type=as.factor(meta_hyp3$Strain.type),
                         Area=as.factor(meta_hyp3$Area),
                         House = as.factor(meta_hyp3$House)
 )
 
 
 plot.new()
 ord<-ordiellipse(NMDS_hyp3, as.factor(meta_hyp3$Area) ,display = "sites", kind ="sd", conf = 0.95, label = T)
 dev.off()
 
 #Reference: http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
 #Data frame df_ell contains values to show ellipses. It is calculated with function veganCovEllipse which is hidden in vegan package. This function is applied to each level of NMDS (group) and it uses also function cov.wt to calculate covariance matrix.
 veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
 {
   theta <- (0:npoints) * 2 * pi/npoints
   Circle <- cbind(cos(theta), sin(theta))
   t(center + scale * t(Circle %*% chol(cov)))
 }
 
 #Generate ellipse points
 df_ell <- data.frame()
 
 for(g in levels(NMDS_hyp3_df$Area)){
   if(g!="" && (g %in% names(ord))){
     
     df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS_hyp3_df[NMDS_hyp3$Area==g,],
                                                      veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                   ,Area=g))
   }
 }
 
 head(df_ell)
 
 #Generate mean values from NMDS plot grouped on strain type
 NMDS.mean=aggregate(NMDS_hyp3_df[,1:2],list(group=NMDS_hyp3_df$Area),mean)
 
 NMDS.mean
 
 #Now do the actual plotting
 
 shape_values<-seq(1,16)
 
 p = ggplot(data=NMDS_hyp3_df,aes(x,y,colour=Area))+
   annotate("text",x=NMDS.mean$x,y=NMDS.mean$y,label=NMDS.mean$group,size=4)+
   geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2), size=1, linetype=2)+
   geom_point(aes(shape=House))+scale_shape_manual(values=shape_values)+theme_bw() 
 
 p
 
 ggsave("C:/Users/Lauren/Dropbox/Pseudomonas/NMDS_hyp3_SNP.png")
 #pdf("NMDS.pdf")
 
 ad_hyp3
 