# load libraries
library(reshape2)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(readxl)

# load input
dset = read.table("unclass_mags_dist.tab") # load mash distance tabular file
meta = read_excel("SuppInfo_metadata.xlsx") # load excel file with metadata

# carry hclust analysis
dist.dset = acast(dset, V1~V2, value.var="V3")
hc = hclust(as.dist(dist.dset))
memb = cutree(hc, h=0.2) # cut tree at 0.2 Mash distance
clusters = unique(memb)
ngroups = length(clusters)
cat("Number of groups:\n")
print(ngroups)

# write names of bins per cluster into files
for (i in clusters){
  if (length(names(which(memb==i))) > 1){
    write.table(names(which(memb==i)), file=paste("clusters/", "cluster_", i, ".txt", sep=""), 
                quote=FALSE, row.names=FALSE, col.names=FALSE)
  } else {
    write.table(names(which(memb==i)), file=paste("clusters/", "single_", i, ".txt", sep=""),
                quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
}

# create dset with counts per cluster
dset_ref = data.frame(matrix(0,ngroups,4))
dset_ref[,1] = clusters
colnames(dset_ref) = c("Cluster", "MAGs", "Samples", "Projects")
for (i in 1:nrow(dset_ref)){
  dset_ref[i,2] = length(which(memb == dset_ref[i,1]))
  runs = strsplit(names(which(memb == dset_ref[i,1])),"_")
  runs = sapply(runs, "[", 2)
  dset_ref[i,3] = nrow(unique(meta[which(meta$run_accession %in% runs),"sample_accession"]))
  dset_ref[i,4] = length(unique(gsub("_.*","",names(which(memb == dset_ref[i,1])))))
}

# check correlation between MAG counts and projects (dset_ref$Projects) or samples (dset$Samples)
#res = cor.test(log10(dset_ref$MAGs), log10(dset_ref$Projects), method = "pearson")
res = cor.test(log10(dset_ref$MAGs), log10(dset_ref$Samples), method = "pearson")
rsqd = res$estimate**2
#print(ggplot(dset_ref, aes(x=log10(dset_ref$MAGs), y=log10(dset_ref$Projects)))
print(ggplot(dset_ref, aes(x=log10(dset_ref$MAGs), y=log10(dset_ref$Samples)))
      + geom_point(size=1) 
      + theme_bw()
      #+ labs(x = "log10 (MAG counts)", y = "log10 (Project counts)")
      + labs(x = "log10 (MAG counts)", y = "log10 (Sample counts)")
      + ylim(-0.1,3)
      + xlim(-0.1,3)
      + geom_smooth(method="lm", color="steelblue")
      + theme(axis.text.x = element_text(face="plain", size=12))
      + theme(axis.text.y = element_text(face="plain", size=12))
      + theme(axis.title.y = element_text(size=14))
      + theme(axis.title.x = element_text(size=14)))
