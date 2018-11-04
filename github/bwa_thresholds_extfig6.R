# merge multiple csv files
library(dplyr)
library(readr)
library(ggplot2)
library(ggrastr)
library(rpgm)

# set working directory
setwd("~/Documents/ESPOD/Analyses/Project_UMGS/MetaSpecies_revision/read_mapping/bwa")

# load files
bwa.depth = read.csv("bwa_depth.csv", row.names=1)
bwa.cov = read.csv("bwa_coverage.csv", row.names=1)
bwa.coeff = read.csv("bwa_coeffvar.csv", row.names=1)
bwa.covscore = (100-bwa.cov)*log10(bwa.depth)
bwa.coeffscore = (100-bwa.cov)*bwa.coeff

cov = as.vector(as.matrix(bwa.cov))[which(bwa.cov >= 60)]
coeff = as.vector(as.matrix(bwa.coeff))[which(bwa.cov >= 60)]
depth = as.vector(as.matrix(bwa.depth))[which(bwa.cov >= 60)]
coeffscore = as.vector(as.matrix(bwa.coeffscore))[which(bwa.cov >= 60)]
covscore = as.vector(as.matrix(bwa.covscore))[which(bwa.cov >= 60)]
dset = data.frame(Cov=cov, Depth=log10(depth), CovScore=covscore, Coeff=coeff, CoeffScore=coeffscore)

# plot results (thresholds: cov >= 60 and coeff < 12)
print(ggplot(dset, aes(x=Cov, y=CovScore))
      + geom_point_rast(size=0.2, alpha=0.2, colour="dodgerblue4", dpi=500)
      + geom_hline(yintercept = quantile(covscore, 0.99), lty = "dotted", size=0.8, colour="red")
      #+ geom_hline(yintercept = quantile(coeffscore, 0.99), lty = "dotted", size=0.8, colour="red")
      #+ geom_smooth(method="lm", size=1)
      + theme_bw()
      #+ ylab("Missing coverage x coefficient of variation")
      + ylab("Missing coverage x log10(mean depth)")
      + xlab("Coverage (%)")
      + theme(axis.title.y = element_text(size=14))
      + theme(axis.text.y = element_text(size=12))
      + theme(axis.title.x = element_text(size=14))
      + theme(axis.text.x = element_text(size=12)))

# define presence/absence
bwa.thresh = bwa.cov
bwa.thresh[bwa.cov >= 60 & 
             bwa.coeffscore < quantile(coeffscore, 0.99) & 
             bwa.covscore < quantile(covscore, 0.99)] = 1
bwa.thresh[bwa.thresh != 1] = 0
length(grep("bin", names(which(rowSums(bwa.thresh) > 0)))) # how many UMGS are detected
bwa.thresh = cbind(rownames(bwa.thresh), bwa.thresh)
colnames(bwa.thresh) = c("Genome", colnames(bwa.thresh)[-1])
#write.csv(bwa.thresh, file="bwa_presence-absence.csv", quote=FALSE, row.names=FALSE)

# check concordance between bins and samples of origin
umgs.hq = scan("../../taxonomy/umgs-hq_genomes.txt", what="")
umgs.mq = scan("../../taxonomy/umgs-mq_genomes.txt", what="")
umgs = c(umgs.hq, umgs.mq)
bwa.thresh = bwa.cov[umgs,]
bwa.thresh.cv = bwa.coeff[umgs,]
control.samp = data.frame(Cov=rep(0,nrow(bwa.thresh)), CV=rep(0,nrow(bwa.thresh)), row.names=rownames(bwa.thresh))
for (n in 1:nrow(control.samp)){
  samp.name = gsub("_.*","",rownames(control.samp)[n])
  control.samp[n,1] = bwa.thresh[n,samp.name]
  control.samp[n,2] = bwa.thresh.cv[n,samp.name]
}