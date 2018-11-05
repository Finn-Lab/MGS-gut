#!/usr/bin/env Rscript

# load libraries
library(VirFinder)
library(optparse)

# prepare arguments
option_list = list(
  make_option(c("-f", "--fasta"), type="character", default=NULL, 
              help="fasta file", metavar="fasta")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$fasta)){
  print_help(opt_parser)
  stop("Provide input FASTA file with -f")
}

# prepare input files
path = normalizePath(opt$fasta)

# run VirFinder
cat("Running VirFinder on:",opt$fasta,"\n")
predResult = VF.pred(path)

# correct for multiple testing
predResult$fdr = p.adjust(predResult$pvalue, method="BH")
predResult$qvalue = VF.qvalue(predResult$pvalue)

# write to output file
cat("Saving output files","\n")
result_table = paste(sub("\\.fa.*", "", opt$fasta), "_VirFinder_results.tab", sep="")
write.table(predResult, file=result_table, sep="\t", quote=FALSE, row.names=FALSE)
