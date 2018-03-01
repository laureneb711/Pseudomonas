#import the gene presence/absence table from roary
gap = read.csv("C:/Users/Lauren/Dropbox/Pseudomonas_fromThunderHorse/gene_presence_absence.csv", header = T, stringsAsFactors = F)

rt = read.table("C:/Users/Lauren/Dropbox/Pseudomonas_fromThunderHorse/gene_presence_absence.rtab")
head(rt)

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

ad_nohyp = adonis(dist ~ Strain.type, meta, perm = 200, method = "binary")

NMDS=metaMDS(dm,k=2,trymax=1000, distance = "binary", wascores = T)

stressplot(NMDS)
plot(NMDS)

NMDS

data.scores <- as.data.frame(scores(NMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$grp <- meta$Area  #  add the grp variable created earlier
head(data.scores)  #look at the data

species.scores <- as.data.frame(scores(NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data
#start hypothesis testing

source("http://bioconductor.org/biocLite.R")
biocLite("topGO")
