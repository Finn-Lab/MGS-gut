# load libraries
library(ggplot2)
library(ggrastr)
library(VennDiagram)
library(reshape2)

# setup input files and datasets
setwd("~/Documents/ESPOD/Analyses/Project_UMGS/MetaSpecies_revision//quality_metrics/")
options(scipen=999)

# load checkm results and rRNA counts
dset = read.delim("checkm_final.tab", row.names=1)
rRNAs.raw = read.delim("rRNAs_id_mags-high.tab", header=FALSE)
rRNAs = acast(rRNAs.raw, V1 ~ V2)
colnames(rRNAs) = c("5S", "23S", "16S")
tRNAs = read.delim("tRNAs_id_mags-high.tab", header=FALSE, row.names=1)
colnames(tRNAs) = "tRNAs"
RNAs = merge(rRNAs, tRNAs, by="row.names")
rownames(RNAs) = RNAs$Row.names
RNAs = RNAs[,-1]
dset = merge(dset, RNAs, by="row.names", all.x=TRUE)

# add quality score
dset$checkm_qs = dset$complet-(5*dset$cont)
#dset = dset[grep("k__Bacteria", dset$taxon),] # bacteria

# classification for colouring
dset$class_qs = rep(0,nrow(dset))
dset$class_qs[which(dset$complet >= 50 & dset$cont < 10)] = "Medium quality (QS <= 50)"
dset$class_qs[which(dset$checkm_qs > 50)] = "Medium quality (QS > 50)"
dset$class_qs[which(dset$complet > 90 & dset$cont < 5)] = "Near complete"
dset$class_qs[which(dset$complet > 90 & dset$cont < 5 &
        dset$`5S` > 80 & dset$`16S` > 80 & dset$`23S` > 80 &
        dset$tRNAs >= 18)] = "High quality"
dset$class_qs[which(dset$class_qs == 0)] = "Low quality"

# plot results (coloured by checkm classification)
print(ggplot(dset, aes(x=complet, y=cont, colour=class_qs)) 
      + geom_point_rast(size=0.2, alpha=0.2, dpi=500)
      + scale_color_manual(name="",limits=c("Low quality", "Medium quality (QS <= 50)", 
                                            "Medium quality (QS > 50)", "Near complete", "High quality"),
                           values=c("grey", "steelblue", "#000033", "darkolivegreen3", "darkgreen"))
      + guides(colour = guide_legend(override.aes = list(size=5, alpha=0.8)))
      + theme_bw()
      + ylab("Contamination (%)")
      + xlab("Completeness (%)")
      + ylim(0,100)
      + theme(axis.title.y = element_text(size=14))
      + theme(axis.text.y = element_text(size=12))
      + theme(axis.title.x = element_text(size=14))
      + theme(axis.text.x = element_text(size=12)))

# calculate total counts per quality category
low = length(which(dset$class_qs=="Low quality"))
med = length(which(dset$class_qs=="Medium quality (QS <= 50)"))
med_qs50 = length(which(dset$class_qs == "Medium quality (QS > 50)"))
med_near = length(which(dset$class_qs == "Near complete"))
high_qual = length(which(dset$class_qs == "High quality"))
df_qual = data.frame(matrix(c(low,med,med_qs50,med_near,high_qual),1))
colnames(df_qual) = c("Low quality", "Medium quality (QS <= 50)", "Medium quality (QS > 50)", 
                      "Near complete", "High quality")
df_stack = melt(df_qual)
df_stack$quality = factor(c("Low", "Medium", "Medium", "High", "High"))

# plot bargraph with counts
print(ggplot(df_stack, aes(x=quality, y=as.numeric(as.character(value)), fill=variable)) 
      + geom_bar(stat="identity", width = 0.6, alpha=0.6)
      + theme_bw()
      + scale_fill_manual(name="",
                          limits=c("Low quality", "Medium quality (QS <= 50)", 
                                   "Medium quality (QS > 50)", "Near complete", "High quality"),
                          values=c("grey", "steelblue", "#000033", "darkgreen", "darkgreen"))
      + scale_x_discrete(limits=c("Low","Medium", "High"),
                                  labels=c("< 50% complet.\n> 10% cont.", 
                                  "> 50% complet.\n< 10% cont.",
                                  "> 90% complet.\n< 5% cont."))
      + ylab("Number of MAGs")
      + theme(axis.title.y = element_text(size=14))
      + theme(axis.text.y = element_text(size=12))
      + theme(axis.title.x = element_blank())
      + theme(axis.text.x = element_text(size=12)))

# saving source data
#source = dset[,1:3]
#colnames(source) = c("Bin", "Completeness", "Contamination")
#write.csv(source, file="../figures/Extended_Data_Figure_2_source.csv", row.names=FALSE, quote=FALSE)
