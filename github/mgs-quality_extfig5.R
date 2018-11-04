# load libraries
library(ggplot2)
library(ggrastr)
library(VennDiagram)
library(reshape2)

# setup input files and datasets
setwd("~/Documents/ESPOD/Analyses/Project_UMGS/MetaSpecies_revision/quality_metrics/")
options(scipen=999)

# load checkm results and rRNA counts
dset = read.delim("checkm_final.tab", row.names=1)
rRNAs.raw = read.delim("rRNAs_id_mgs.tab", header=FALSE)
rRNAs = acast(rRNAs.raw, V1 ~ V2)
colnames(rRNAs) = c("5S", "23S", "16S")
tRNAs = read.delim("tRNAs_id_mgs.tab", header=FALSE, row.names=1)
colnames(tRNAs) = "tRNAs"
RNAs = merge(rRNAs, tRNAs, by="row.names")
rownames(RNAs) = RNAs$Row.names
RNAs = RNAs[,-1]
mgs = scan("../taxonomy/mgs_genomes.txt", what="")
dset = dset[mgs,]
dset = merge(dset, RNAs, by="row.names")
rownames(dset) = dset$Row.names
dset = dset[,-1]
markers = read.delim("checkm_marker-genes.tab", row.names=1)
dset = merge(dset, markers, by="row.names")
dset$Marker_cont = dset$Markers.over.1/dset$Total.Markers*100
#genomes = scan("drep_mgs-hq.txt", what="")
#genomes = scan("drep_mgs-mq.txt", what="")
#dset = dset[which(dset$Row.names %in% genomes),]

# table formatting
dset.boxplot = dset[,c("complet", "cont", "Marker_cont")]
dset.boxplot = melt(dset.boxplot)

# plot boxplot of completeness and purity (100-contamination)
plot = print(ggplot(dset.boxplot, aes(x=variable, y=value, colour=variable, fill=variable)) 
             #+ geom_violin(alpha=0.5)
             + geom_point(position = position_jitter(width=0.2, height=0), alpha=0.5, size=0.3)
             + geom_boxplot(alpha=0.5, width=0.5, outlier.colour=NA, colour="black")
             #+ geom_density(position = "identity", alpha = 0.6)
             + ylab("%")
             + theme_bw()
             + guides(fill=FALSE)
             + scale_y_continuous(breaks = c(0,20,40,60,80,100))
             + ylim(0,100)
             + scale_x_discrete(labels=c("Completeness", "Contamination", "Multiple MGs"))
             + theme(axis.title.y = element_text(size=14))
             + theme(axis.text.y = element_text(size=12))
             + theme(axis.title.x = element_blank())
             + scale_colour_manual(values=c("darkgreen", "red3", "steelblue"))
             + scale_fill_manual(values=c("darkgreen", "red3", "steelblue"))
             + guides(colour=FALSE)
             + theme(axis.text.x = element_text(size=12)))

# plot histogram for one variable
dset.hist = data.frame(tRNAs=dset[,"tRNAs"], row.names=dset$Row.names)
print(ggplot(dset.hist, aes(x=tRNAs)) 
      + geom_histogram(colour="white", fill="darkgrey", alpha=0.8, bins=18)
      + theme_bw()
      + ylab("Frequency")
      + xlab("Number of tRNAs")
      + guides(fill=FALSE)
      #+ scale_fill_manual(name="",values=c("grey", "#66CCFF", "#000033"))
      #+ xlim(60,100) # ANI
      + scale_fill_manual(name="",
                          values=c("gray", "darkgreen"),
                          limits=c("Unclassified", "Reference match"))
      + theme(axis.title.y = element_text(size=14))
      + theme(axis.text.y = element_text(size=12))
      + theme(axis.title.x = element_text(size=14))
      + theme(axis.text.x = element_text(size=12)))

# get boxplot stats
ggplot_build(plot)$data