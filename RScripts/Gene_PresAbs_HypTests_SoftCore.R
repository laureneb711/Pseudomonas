#import the gene presence/absence table from roary
gap = read.csv("C:/Users/Lauren/Dropbox/Pseudomonas_fromThunderHorse/gene_presence_absence.csv", header = T, stringsAsFactors = F)

#Roary defines "SoftCore" as genes found in >44 genomes ad <47 
gap = subset(gap, gap$No..isolates >= 44 &
               gap$No..isolates < 47)

#rt = read.table("C:/Users/Lauren/Dropbox/Pseudomonas_fromThunderHorse/gene_presence_absence.rtab")

#remove hypothetical proteins from the file
no_hp = subset(gap, gap$Annotation != "hypothetical protein")

no_hp_table = as.data.frame(table(no_hp$Annotation), stringsAsFactors = F)

colnames(no_hp)

presabs = no_hp[,15:61]

presabs = replace(presabs, presabs !="", 1)
presabs = replace(presabs, presabs =="", 0)

colnames(presabs)
rownames(presabs)

presabs_t = as.data.frame(t(presabs))

colnames(presabs_t)
rownames(presabs_t)

tail(presabs_t[,1:5])

dm = as.matrix(dist(presabs_t, method = "binary"))

dist = dist(presabs_t, method = "binary")
head(dist)

library(vegan)
?adonis

meta = read.csv("C:/Users/Lauren/Dropbox/Pseudomonas/Files/Pseudomonas_metadata.csv", header = TRUE)

meta = subset(meta, meta$Area != "Test")
meta$Area
meta$Strain.type


############
#looking for difference without any specific hypothesis
ad_nohyp = adonis(dist ~ Strain.type, meta, perm = 200, method = "binary")

NMDS=metaMDS(dm,k=2,trymax=1000, distance = "binary", wascores = T)

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

#find which ones were removed
meta$sample.id[which(is.na(match(meta$sample.id, meta_hyp1$sample.id)))]

#remove same ones from rownames of presabs_t      
presabs_t_hyp1 = presabs_t[-which(is.na(match(meta$sample.id, meta_hyp1$sample.id))),]
head(presabs_t_hyp1[1:4])

rownames(presabs_t_hyp1)

dm_hyp1 = as.matrix(dist(presabs_t_hyp1, method = "binary"))

dist_hyp1 = dist(presabs_t_hyp1, method = "binary")
head(dist_hyp1)

ad_hyp1 = adonis(dist_hyp1 ~ Strain.type, meta_hyp1, perm = 200, method = "binary")

NMDS_hyp1=metaMDS(dm_hyp1,k=2,trymax=1000, distance = "binary", wascores = T)

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

ggsave("C:/Users/Lauren/Dropbox/Pseudomonas/NMDS_hyp1_SoftCore.png")
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

meta_hyp2$House

#find which ones were removed
#remove same ones from rownames of presabs_t      
presabs_t_hyp2 = presabs_t[-which(is.na(match(meta$sample.id, meta_hyp2$sample.id))),]
head(presabs_t_hyp2[1:4])

rownames(presabs_t_hyp2)

dm_hyp2 = as.matrix(dist(presabs_t_hyp2, method = "binary"))

dist_hyp2 = dist(presabs_t_hyp2, method = "binary")
head(dist_hyp2)

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

ggsave("C:/Users/Lauren/Dropbox/Pseudomonas/NMDS_hyp2_SoftCore.png")

ad_hyp2

############
#Repeat the process, but this time address hyp 3
#hyp 3: Comparative genomics of household strains in light of 3 kinds of drains from which they were isolated 
#(kitchen, bathroom, sink); 
#is there an association between sampling location and pattern of relatedness?
############

#first, remove 3 healthy strains from the analysis
meta_hyp3 = subset(meta, meta$Strain.type == "drain")

#find which ones were removed
meta$sample.id[which(is.na(match(meta$sample.id, meta_hyp3$sample.id)))]

#remove same ones from rownames of presabs_t      
presabs_t_hyp3 = presabs_t[-which(is.na(match(meta$sample.id, meta_hyp3$sample.id))),]
head(presabs_t_hyp3[1:4])

rownames(presabs_t_hyp3)

dm_hyp3 = as.matrix(dist(presabs_t_hyp3, method = "binary"))

dist_hyp3 = dist(presabs_t_hyp3, method = "binary")
head(dist_hyp3)

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

ggsave("C:/Users/Lauren/Dropbox/Pseudomonas/NMDS_hyp3_SoftCore.png")
#pdf("NMDS.pdf")

ad_hyp3
