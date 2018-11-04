# load libraries
library(ggplot2)
library(VennDiagram)

# setup input files and datasets
setwd("~/Documents/ESPOD/Analyses/Assemb_Binning/MetaSpecies_revision/dnadiff/")
col.template = c("_hit", "_lenref", "_aliref", "_lenquer", "_aliquer", "_ani", "_mash_dist")
hgr = read.delim("dnadiff_hgr_final.tab", row.names=1)
colnames(hgr) = paste("hgr", col.template, sep="")
refseq = read.delim("dnadiff_refseq_final.tab", row.names=1)
colnames(refseq) = paste("refseq", col.template, sep="")
checkm = read.delim("../quality_metrics/checkm_final.tab", row.names=1)
dset = merge(hgr, refseq, by="row.names")
rownames(dset) = dset$Row.names
dset = dset[,-1]
mgs = read.delim("dnadiff_mgs_final.tab", row.names=1) # MGS
colnames(mgs) = paste("mgs", col.template, sep="") # MGS
dset = merge(dset, mgs, by="row.names") # MGS
rownames(dset) = dset$Row.names # MGS
dset = dset[,-1] # MGS
dset = merge(dset, checkm, by="row.names")
colnames(dset) = c("bins", colnames(dset)[-1])

# add quality score
dset$checkm_qs = dset$complet-(5*dset$cont)

# filter data
dset = dset[which(dset$checkm_qs > 50),] # QS > 50
dset = dset[which(dset$complet > 90 & dset$cont < 5),]
dset = dset[grep("k__Bacteria", dset$taxon),] # bacteria

# counts
#vars = c("HGR", "RefSeq", "Unknown")
vars = c("HGR", "RefSeq", "MGS", "Unknown") # MGS
counts = rep(0,length(vars))
counts[1] = length(which(dset$hgr_aliquer >= 60 & dset$hgr_ani >= 95))
counts[2] = length(which(dset$refseq_aliquer >= 60 & dset$refseq_ani >= 95))
counts[3] = length(which(dset$mgs_aliquer >= 60 & dset$mgs_ani >= 95)) # MGS
counts[4] = nrow(dset) - length(which(dset$refseq_aliquer >= 60 & dset$refseq_ani >= 95 # MGS
                                      | dset$hgr_aliquer >= 60 & dset$hgr_ani >= 95
                                      | dset$mgs_aliquer >= 60 & dset$mgs_ani >= 95)) # MGS
#counts[3] = nrow(dset) - length(which(dset$refseq_aliquer >= 60 & dset$refseq_ani >= 95 
#                         | dset$hgr_aliquer >= 60 & dset$hgr_ani >= 95))
df = data.frame(vars, counts)

# plot bargraph of match counts
print(ggplot(df, aes(x=vars, y=counts))
      + geom_bar(stat="identity", fill="#000033", alpha=0.6)
      + scale_x_discrete(limits=df$vars)
      + theme_bw()
      + ylab("Number of MAGs")
      + theme(axis.title.y = element_text(size=14))
      + theme(axis.text.y = element_text(size=12))
      + theme(axis.title.x = element_blank())
      + theme(axis.text.x = element_text(size=12)))

# get counts for venn diagram
hgr_match = length(which(dset$hgr_aliquer >= 60 & dset$hgr_ani >= 95))
refseq_match = length(which(dset$refseq_aliquer >= 60 & dset$refseq_ani >= 95))
hgr_refseq = length(which(dset$hgr_aliquer >= 60 & dset$hgr_ani >= 95
                          & dset$refseq_aliquer >= 60 & dset$refseq_ani >= 95))
mgs_match = length(which(dset$mgs_aliquer >= 60 & dset$mgs_ani >= 95))

# double venn with HGR and RefSeq
dev.off()
draw.pairwise.venn(area1 = hgr_match, area2 = refseq_match, cross.area=hgr_refseq,
                 category = c("HR", "RefSeq"),
                 lty = rep("blank", 2), fill = c("steelblue", "red3"), 
                 alpha = rep(0.3, 2),
                 cex=1.5, fontfamily = rep("sans", 3),
                 cat.cex = 1.7, cat.fontfamily = rep("sans", 2), 
                 euler.d = TRUE, scaled = TRUE)

# write unknown bins to file
known = which(dset$refseq_aliquer >= 60 & dset$refseq_ani >= 95 
              | dset$hgr_aliquer >= 60 & dset$hgr_ani >= 95
              | dset$mgs_aliquer >= 60 & dset$mgs_ani >= 95) # MGS
unkn = dset[-known,]
#write.table(unkn$bins, file="unclass_mags_hq.txt", row.names=FALSE, 
#            col.names=FALSE, quote=FALSE)
#write.table(unkn$bins, file="unclass_mags_mq.txt", row.names=FALSE, 
#            col.names=FALSE, quote=FALSE)