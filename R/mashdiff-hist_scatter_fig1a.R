# load libraries
library(ggplot2)
library(ggrastr)
library(reshape2)

# setup input files and datasets
dset = read.delim("mashdiff_results.tab") # load file with mashdiff results with HR and RefSeq

# add quality score
dset$checkm_qs = dset$complet-(5*dset$cont)

# filter data (qs and taxon)
dset = dset[which(dset$checkm_qs > 50),]
dset = dset[grep("k__Bacteria", dset$taxon),]

# classify as matched
dset$hr_match = rep(0,nrow(dset))
dset$hr_match[which(dset$hr_aliquer >= 60 & dset$hr_ani >= 95)] = "Reference match"
dset$hr_match[which(dset$hr_match == 0)] = "Unclassified"

dset$hr_clani = rep(0,nrow(dset))
dset$hr_clani[which(dset$hr_ani >= 95)] = ">= 95% ANI"
dset$hr_clani[which(dset$hr_ani >= 99)] = ">= 99% ANI"
dset$hr_clani[which(dset$hr_clani == 0)] = "< 95% ANI"

# plot results
print(ggplot(dset, aes(x=hr_aliquer, y=hr_aliref, colour=hr_match)) 
      + geom_point_rast(size=0.2, alpha=0.1, dpi = 500)
      + scale_color_manual(name="",
                           values=c("gray", "darkgreen"),
                           limits=c("Unclassified", "Reference match"))
      + guides(colour = guide_legend(override.aes = list(size=5, alpha=0.8)))
      + theme_bw()
      + ylab("% Reference aligned")
      + xlab("% MAG aligned")
      + xlim(0,100)
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
      + scale_fill_manual(name="",
                          values=c("gray", "darkgreen"),
                          limits=c("Unclassified", "Reference match"))
      + theme(axis.title.y = element_text(size=14))
      + theme(axis.text.y = element_text(size=10))
      + theme(axis.title.x = element_text(size=14))
      + theme(axis.text.x = element_text(size=10))
      + theme(legend.position = "top"))
