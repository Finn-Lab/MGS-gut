# load libraries
library(ggplot2)
library(ggrastr)
library(VennDiagram)
library(reshape2)

# load input
dset = read.delim("mgs-quality.tab", row.names=1) # load file with checkm quality and RNA counts

# table formatting
dset.boxplot = dset[,c("complet", "cont")]
dset.boxplot = melt(dset.boxplot)

# plot boxplot of completeness and contamination
plot = print(ggplot(dset.boxplot, aes(x=variable, y=value, colour=variable, fill=variable)) 
             + geom_point(position = position_jitter(width=0.2, height=0), alpha=0.5, size=0.3)
             + geom_boxplot(alpha=0.5, width=0.5, outlier.colour=NA, colour="black")
             + ylab("%")
             + theme_bw()
             + guides(fill=FALSE)
             + scale_y_continuous(breaks = c(0,20,40,60,80,100))
             + ylim(0,100)
             + scale_x_discrete(labels=c("Completeness", "Contamination"))
             + theme(axis.title.y = element_text(size=14))
             + theme(axis.text.y = element_text(size=12))
             + theme(axis.title.x = element_blank())
             + scale_colour_manual(values=c("darkgreen", "red3", "steelblue"))
             + scale_fill_manual(values=c("darkgreen", "red3", "steelblue"))
             + guides(colour=FALSE)
             + theme(axis.text.x = element_text(size=12)))

# plot histogram of tRNA counts
dset.hist = data.frame(tRNAs=dset[,"tRNAs"], row.names=dset$Row.names)
print(ggplot(dset.hist, aes(x=tRNAs)) 
      + geom_histogram(colour="white", fill="darkgrey", alpha=0.8, bins=18)
      + theme_bw()
      + ylab("Frequency")
      + xlab("Number of tRNAs")
      + guides(fill=FALSE)
      + scale_fill_manual(name="",
                          values=c("gray", "darkgreen"),
                          limits=c("Unclassified", "Reference match"))
      + theme(axis.title.y = element_text(size=14))
      + theme(axis.text.y = element_text(size=12))
      + theme(axis.title.x = element_text(size=14))
      + theme(axis.text.x = element_text(size=12)))

# get boxplot stats
ggplot_build(plot)$data
