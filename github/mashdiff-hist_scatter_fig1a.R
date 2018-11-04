# load libraries
library(ggplot2)
library(ggrastr)
library(reshape2)

# setup input files and datasets
setwd("~/Documents/ESPOD/Analyses/Project_UMGS//MetaSpecies_revision/dnadiff/")
col.template = c("_hit", "_lenref", "_aliref", "_lenquer", "_aliquer", "_ani", "_mash_dist")
hgr = read.delim("dnadiff_hgr_final.tab", row.names=1)
colnames(hgr) = paste("hgr", col.template, sep="")
refseq = read.delim("dnadiff_refseq_final.tab", row.names=1)
colnames(refseq) = paste("refseq", col.template, sep="")
checkm = read.delim("../quality_metrics/checkm_final.tab", row.names=1)
dset = merge(hgr, refseq, by="row.names")
rownames(dset) = dset$Row.names
dset = dset[,-1]
dset = merge(dset, checkm, by="row.names")
colnames(dset) = c("bins", colnames(dset)[-1])

# add quality score
dset$checkm_qs = dset$complet-(5*dset$cont)

# filter data
#dset = dset[which(dset$checkm_qs > 50 & dset$cont < 5),] # QS > 50
dset = dset[which(dset$complet > 90 & dset$cont < 5),]
dset = dset[grep("k__Bacteria", dset$taxon),] # bacteria

# classify as matched
dset$hgr_match = rep(0,nrow(dset))
dset$hgr_match[which(dset$hgr_aliquer >= 60 & dset$hgr_ani >= 95)] = "Reference match"
dset$hgr_match[which(dset$hgr_match == 0)] = "Unclassified"

dset$hgr_clani = rep(0,nrow(dset))
dset$hgr_clani[which(dset$hgr_ani >= 95)] = ">= 95% ANI"
dset$hgr_clani[which(dset$hgr_ani >= 99)] = ">= 99% ANI"
dset$hgr_clani[which(dset$hgr_clani == 0)] = "< 95% ANI"

# plot results
print(ggplot(dset, aes(x=hgr_aliquer, y=hgr_aliref, colour=hgr_match)) 
      + geom_point_rast(size=0.2, alpha=0.1, dpi = 500)
      #+ scale_color_manual(name="",values=c("grey", "#66CCFF", "#000033")) # ANI
      + scale_color_manual(name="",
                           values=c("gray", "darkgreen"),
                           limits=c("Unclassified", "Reference match"))
      + guides(colour = guide_legend(override.aes = list(size=5, alpha=0.8)))
      + theme_bw()
      + ylab("% Reference aligned")
      + xlab("% MAG aligned")
      + xlim(0,100)
      #+ xlim(60,100) # ANI
      #+ ylim(60,100) # ANI
      + ylim(0,100)
      + theme(axis.title.y = element_text(size=14))
      + theme(axis.text.y = element_text(size=12))
      + theme(axis.title.x = element_text(size=14))
      + theme(axis.text.x = element_text(size=12)))

# plot histogram for one variable
plot_dset = melt(dset)
plot_dset = plot_dset[which(plot_dset$variable=="hgr_aliquer"),] # subset by variable
print(ggplot(plot_dset, aes(x=value, fill=hgr_match)) 
      + geom_histogram(colour="black", alpha=0.5, bins=50, size=0.2)
      + theme_bw()
      + ylab("Number of MAGs")
      + xlab("% MAG aligned")
      # guides(fill=FALSE)
      #+ scale_fill_manual(name="",values=c("grey", "#66CCFF", "#000033"))
      #+ xlim(60,100) # ANI
      + scale_fill_manual(name="",
                          values=c("gray", "darkgreen"),
                          limits=c("Unclassified", "Reference match"))
      + theme(axis.title.y = element_text(size=14))
      + theme(axis.text.y = element_text(size=10))
      + theme(axis.title.x = element_text(size=14))
      + theme(axis.text.x = element_text(size=10))
      + theme(legend.position = "top"))

# save supplementary information
supp.dset = dset[,c("bins", "complet", "cont", "strain", "checkm_qs", "taxon")]
colnames(supp.dset) = c("Bin", "Completeness", "Contamination", "Strain_Heterogeneity",
                        "Quality Score (QS)", "CheckM lineage")
write.table(supp.dset, file="../tables/40K_mags.tab", sep= "\t", quote=FALSE, row.names=FALSE)

# save source file (Fig. 1a)
fig1.dset = dset[,c("bins", "hgr_aliref", "hgr_aliquer", "hgr_ani")]
colnames(fig1.dset) = c("Bin", "% Reference aligned", "% MAG aligned", "ANI (%)")
write.csv(fig1.dset, file="../figures/Figure_1a_source.csv", quote=FALSE, row.names=FALSE)