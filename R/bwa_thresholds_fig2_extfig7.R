# merge multiple csv files
library(dplyr)
library(readr)
library(ggplot2)
library(ggrastr)
library(rpgm)

# load files
bwa.depth = read.csv("bwa_depth.csv", row.names=1)
bwa.cov = read.csv("bwa_coverage.csv", row.names=1)
bwa.coeff = read.csv("bwa_coeffvar.csv", row.names=1)
bwa.covscore = (100-bwa.cov)*log10(bwa.depth)
bwa.coeffscore = (100-bwa.cov)*bwa.coeff

# calculate parameters
cov = as.vector(as.matrix(bwa.cov))[which(bwa.cov >= 60)]
coeff = as.vector(as.matrix(bwa.coeff))[which(bwa.cov >= 60)]
depth = as.vector(as.matrix(bwa.depth))[which(bwa.cov >= 60)]
coeffscore = as.vector(as.matrix(bwa.coeffscore))[which(bwa.cov >= 60)]
covscore = as.vector(as.matrix(bwa.covscore))[which(bwa.cov >= 60)]
dset = data.frame(Cov=cov, Depth=log10(depth), CovScore=covscore, Coeff=coeff, CoeffScore=coeffscore)

# visualize scatterplot
print(ggplot(dset, aes(x=Cov, y=CovScore))
      + geom_point_rast(size=0.2, alpha=0.2, colour="dodgerblue4", dpi=500)
      + geom_hline(yintercept = quantile(covscore, 0.99), lty = "dotted", size=0.8, colour="red")
      #+ geom_hline(yintercept = quantile(coeffscore, 0.99), lty = "dotted", size=0.8, colour="red")
      + theme_bw()
      + ylab("Coverage penalty score")
      #+ ylab("Variation penalty score")
      + xlab("Coverage (%)")
      + theme(axis.title.y = element_text(size=14))
      + theme(axis.text.y = element_text(size=12))
      + theme(axis.title.x = element_text(size=14))
      + theme(axis.text.x = element_text(size=12)))

# define presence/absence and save matrix to csv file
bwa.thresh = bwa.cov
bwa.thresh[bwa.cov >= 60 & 
             bwa.coeffscore < quantile(coeffscore, 0.99) & 
             bwa.covscore < quantile(covscore, 0.99)] = 1
bwa.thresh[bwa.thresh != 1] = 0
bwa.thresh = cbind(rownames(bwa.thresh), bwa.thresh)
colnames(bwa.thresh) = c("Genome", colnames(bwa.thresh)[-1])
write.csv(bwa.thresh, file="bwa_presence-absence.csv", quote=FALSE, row.names=FALSE)