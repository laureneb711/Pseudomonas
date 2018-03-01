library(seqminer)

readVCFToMatrixByGene("C:/Users/Lauren/Dropbox/Pseudomonas_fromThunderHorse/VCFs/NoMods/psamplestrim.vcf")

#import the gene presence/absence table from roary
gap = read.csv("C:/Users/Lauren/Dropbox/Pseudomonas_fromThunderHorse/gene_presence_absence.csv", header = T, stringsAsFactors = F)

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

ad_nohyp = adonis(dist ~ Strain.type, meta, perm = 200)

NMDS=metaMDS(dm,k=2,trymax=1000)

stressplot(NMDS)
plot(NMDS)


#start hypothesis testing