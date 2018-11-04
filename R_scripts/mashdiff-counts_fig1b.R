# load libraries
library(ggplot2)
library(VennDiagram)

# setup input files and datasets
dset = read.delim("mashdiff_results.tab") # load file with mashdiff results with HR and RefSeq

# add quality score
dset$checkm_qs = dset$complet-(5*dset$cont)

# filter data
dset = dset[which(dset$checkm_qs > 50),] # QS > 50
dset = dset[grep("k__Bacteria", dset$taxon),] # bacteria

# counts
vars = c("HR", "RefSeq", "Unknown")
counts = rep(0,length(vars))
counts[1] = length(which(dset$hr_aliquer >= 60 & dset$hr_ani >= 95))
counts[2] = length(which(dset$refseq_aliquer >= 60 & dset$refseq_ani >= 95))
counts[3] = nrow(dset) - length(which(dset$refseq_aliquer >= 60 & dset$refseq_ani >= 95 
                         | dset$hr_aliquer >= 60 & dset$hr_ani >= 95))
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
hr_match = length(which(dset$hr_aliquer >= 60 & dset$hr_ani >= 95))
refseq_match = length(which(dset$refseq_aliquer >= 60 & dset$refseq_ani >= 95))
hr_refseq = length(which(dset$hr_aliquer >= 60 & dset$hr_ani >= 95
                          & dset$refseq_aliquer >= 60 & dset$refseq_ani >= 95))

# Venn with HR and RefSeq
dev.off()
draw.pairwise.venn(area1 = hr_match, area2 = refseq_match, cross.area=hr_refseq,
                 category = c("HR", "RefSeq"),
                 lty = rep("blank", 2), fill = c("steelblue", "red3"), 
                 alpha = rep(0.3, 2),
                 cex=1.5, fontfamily = rep("sans", 3),
                 cat.cex = 1.7, cat.fontfamily = rep("sans", 2), 
                 euler.d = TRUE, scaled = TRUE)
