#!/usr/bin/env Rscript

# load libraries
library(dplyr)
library(readr)
library(optparse)

# prepare arguments
option_list = list(
  make_option(c("-t", "--total"), type="character", default=NULL,
              help="total counts directory", metavar="total"),
  make_option(c("-u", "--unique"), type="character", default=NULL,
              help="unique counts directory", metavar="unique"))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$total) || is.null(opt$unique)){
  print_help(opt_parser)
  stop("Provide input directories with -t and -u")
}

# prepare function
merge.sm=function(ext) {
  samp_list = Sys.glob(ext)
  datalist= lapply(samp_list, function (x) read.delim(x, row.names=1))
  do.call(cbind, datalist)}

# load files
cat("Reading BWA results (total) ...\n")
bwa.data = merge.sm(paste(opt$total,"/*total.tab", sep=""))

cat("Reading BWA results (unique) ...\n")
bwa.unique = merge.sm(paste(opt$unique,"/*unique.tab", sep=""))

# parse and save dsets
cat("Splitting files...\n")
bwa.counts = bwa.data[,grep("Counts", colnames(bwa.data))]
colnames(bwa.counts) = gsub("_.*", "", colnames(bwa.counts))

bwa.depth = bwa.data[,grep("MeanDepth", colnames(bwa.data))]
colnames(bwa.depth) = gsub("_.*", "", colnames(bwa.depth))

bwa.cov = bwa.data[,grep("Coverage", colnames(bwa.data))]
colnames(bwa.cov) = gsub("_.*", "", colnames(bwa.cov))

bwa.coeff = bwa.data[,grep("CoeffVar", colnames(bwa.data))]
colnames(bwa.coeff) = gsub("_.*", "", colnames(bwa.coeff))

bwa.counts.unique = bwa.unique[,grep("Counts", colnames(bwa.unique))]
colnames(bwa.counts.unique) = gsub("_.*", "", colnames(bwa.counts.unique))

# define thresholds for presence/absence
cat("Defining thresholds for presence/absence...\n")
bwa.covscore = (100-bwa.cov)*log10(bwa.depth)
bwa.coeffscore = (100-bwa.cov)*bwa.coeff
coeffscore = as.vector(as.matrix(bwa.coeffscore))[which(bwa.cov >= 60)]
covscore = as.vector(as.matrix(bwa.covscore))[which(bwa.cov >= 60)]

bwa.thresh = data.frame(matrix(0,nrow(bwa.cov),ncol(bwa.cov)))
rownames(bwa.thresh) = rownames(bwa.cov)
colnames(bwa.thresh) = colnames(bwa.cov)
bwa.thresh[bwa.cov >= 60 & 
             bwa.coeffscore < quantile(coeffscore, 0.99) & 
             bwa.covscore < quantile(covscore, 0.99)] = 1

# save output files
cat("Saving files...\n")
bwa.thresh = cbind(rownames(bwa.thresh), bwa.thresh)
colnames(bwa.thresh) = c("Genome", colnames(bwa.thresh)[-1])
write.csv(bwa.thresh, file="bwa_presence-absence.csv", quote=FALSE, row.names=FALSE)

bwa.counts.unique = cbind(rownames(bwa.counts.unique), bwa.counts.unique)
colnames(bwa.counts.unique) = c("Genome", colnames(bwa.counts.unique)[-1])
write.csv(bwa.counts.unique, file="bwa_counts-unique.csv", quote=FALSE, row.names=FALSE)
