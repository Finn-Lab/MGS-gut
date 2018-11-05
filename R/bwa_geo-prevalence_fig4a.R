# load libraries
library(readxl)
library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(grid)

# load files
bwa.prev = read.csv("bwa_presence-absence.csv", row.names=1) # csv file with a presence/absence binary matrix
bwa.counts = read.csv("bwa_counts-unique.csv", row.names=1) # csv file with bwa counts
metadata.raw = read_excel("metadata.xlsx") # file with metadata
metadata = as.data.frame(metadata.raw[,c("pub_state", "pub_disease", "pub_disease_secondary", "pub_agestrat", 
                                         "pub_antibio", "country", "continent", "read_count_total")])
rownames(metadata) = metadata.raw$run_accession
colnames(metadata) = c("disease_state", "disease_name", "disease_secondary", "age", 
                       "antibio", "country", "continent", "read_count")

# create dataset
dset = data.frame(Species=rownames(bwa.prev), Prev_Past=0, Prop=0, Abund=0, Continent=NA)
dset = dset[rep(1:nrow(dset),each=length(unique(metadata$continent))),]
dset$Continent = rep(unique(metadata$continent), nrow(bwa.prev))
rel.ab = t(t(bwa.counts)/metadata$read_count*100)
rel.ab[bwa.prev == 0] = 0
rel.ab.thresh = 0.01
for (n in 1:nrow(dset)){
  cat("Running",n,"\n")
  species = as.vector(dset[n,1])
  total.samples = rownames(metadata[which(metadata$continent == dset[n,5]),])
  present.past = which(rel.ab[species,total.samples] > rel.ab.thresh)
  dset[n,2] = length(present.past)
  dset[n,3] = dset[n,2]/length(total.samples)*100
  dset[n,4] = mean(rel.ab[species,names(present.past)]) 
}
dset$Genome = NA
dset$Genome[grep("bin", dset$Species)] = "UMGS"
dset$Genome[grep("bin", dset$Species, invert=TRUE)] = "HGR"
dset = dset[dset$Continent != "NA",]

# plot samples with species > 0.01 abundance
print(ggplot(dset, aes(x=Continent, y=log10(Prev_Past), fill=Genome))
      + geom_boxplot(outlier.size = 0.1, alpha=0.6, position=position_dodge(0.6), width=0.5)
      + theme_bw()
      + xlab("")
      + ylab("Number of samples (log10)")
      + scale_x_discrete(limits=c("North America", "Europe", "Asia", "Oceania", "South America", "Africa"),
                         labels=c("North\nAmerica", "Europe", "Asia", "Oceania", "South\nAmerica", "Africa"))
      + scale_fill_manual(values=c("steelblue", "darkgreen"))
      + theme(axis.text.y = element_text(size=11))
      + theme(axis.text.x = element_text(size=11)))

# calculate number of species in proportion of samples
dset.quant = data.frame(Genome=c(rep("UMGS", 7),rep("HGR",7)), Number=0, Prop=0, Continent=NA)
dset.quant$Continent = rep(unique(metadata$continent), 2)
percent = 20
for (n in 1:nrow(dset.quant)){
  genome = as.vector(dset.quant[n,1])
  dset.quant[n,2] = length(which(dset$Genome == genome & 
                                   dset$Continent == dset.quant[n,4] & 
                                   dset$Prop > percent))
}

#  plot species present in > 20% samples
print(ggplot(dset.quant, aes(x=Continent, y=Number, fill=Genome))
      + geom_bar(stat="identity", alpha=0.7, colour="black", size=0.2)
      + theme_bw()
      + xlab("")
      + ylab("Number of species in > 20% samples")
      + scale_x_discrete(limits=c("North America", "Europe", "Asia", "Oceania", "South America", "Africa"),
                         labels=c("North\nAmerica", "Europe", "Asia", "Oceania", "South\nAmerica", "Africa"))
      + scale_fill_manual(values=c("steelblue", "darkgreen"))
      + theme(axis.text.y = element_text(size=11))
      + theme(axis.text.x = element_text(size=11)))
