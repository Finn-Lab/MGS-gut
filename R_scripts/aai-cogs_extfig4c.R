# load libraries
library(ggplot2)
setwd("~/Documents/ESPOD/Analyses/Project_UMGS/MetaSpecies_revision/pairwise_aai/")
options(scipen=999)

# load input
aai = data.frame(AAI=as.numeric(scan("mean_aais.txt", what="")))

# plots
print(ggplot(aai, aes(x=AAI))
      + geom_histogram(bins=50, fill="grey", colour="black", alpha=0.5, size=0.2)
      + theme_bw()
      + xlab("Average Amino Acid Identity (%)")
      + ylab("Number of pairwise comparisons")
      + theme(axis.text.x = element_text(size=10))
      + theme(axis.text.y = element_text(size=10))
      + theme(axis.title.x = element_text(size=14))
      + theme(axis.title.y = element_text(size=14)))