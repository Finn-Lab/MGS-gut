# load libraries
library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(grid)

# load input
as.umgs = read.delim("antismash_bgc_umgs.tab", row.names=1, header=FALSE) # tabular file with antismash counts
as.hgr = read.delim("antismash_bgc_hgr.tab", row.names=1, header=FALSE) # tabular file with antismash counts
colnames(as.umgs) = c("Total_BGC", "Known_BGC")
colnames(as.hgr) = c("Total_BGC", "Known_BGC")

# 1. calculate totals
total.umgs = sum(as.umgs[,1])
total.hgr = sum(as.hgr[,1])
novel.umgs = (total.umgs-sum(as.umgs[,2]))/total.umgs*100
novel.hgr = (total.hgr-sum(as.hgr[,2]))/total.hgr*100

novel = data.frame(rbind(novel.umgs, novel.hgr))
rownames(novel) = c("UMGS", "HGR")
colnames(novel) = "Percent"

# plot novel BGC
print(ggplot(novel, aes(x=rownames(novel), y=Percent, fill=rownames(novel))) 
      + geom_bar(stat="identity", width = 0.6, alpha=0.5)
      + theme_bw()
      + coord_flip()
      + ylab("% Novel BGC")
      + ylim(0,100)
      + scale_fill_manual(name="Genome", values=c("steelblue", "red3"))
      + theme(axis.title.x = element_text(size=14))
      + theme(axis.text.y = element_text(size=12))
      + theme(axis.title.y = element_blank())
      + theme(axis.text.x = element_text(size=12)))

# 2. calculate BGC per category
cats.umgs = read.delim("antismash_cats_umgscounts.tab", header=FALSE)
colnames(cats.umgs) = c("Name", "Category", "Count")
cats.hgr = read.delim("antismash_cats_hgrcounts.tab", header=FALSE)
colnames(cats.hgr) = c("Name", "Category", "Count")
cats = as.vector(unique(cats.umgs$Category))
umgs.tax = read.delim("taxumgs.tab", header=FALSE)
colnames(umgs.tax) = c("Name", "Tax", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "UMGS")
hgr.tax = read.delim("taxhgr.tab", header=FALSE)
colnames(hgr.tax) = c("Name", "Tax", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# UMGS
cats.umgs.dset = data.frame(cbind(cats, rep(0, length(cats))), stringsAsFactors = FALSE)
colnames(cats.umgs.dset) = c("Category", "Count")
cats.umgs.dset$Genome = "UMGS"

for (i in 1:nrow(cats.umgs.dset)){
  cats.umgs.dset[i,2] = sum(cats.umgs[which(cats.umgs$Category == cats.umgs.dset[i,1]),3])
}

# HGR
cats.hgr.dset = data.frame(cbind(cats, rep(0, length(cats))), stringsAsFactors = FALSE)
colnames(cats.hgr.dset) = c("Category", "Count")
cats.hgr.dset$Genome = "HGR"

for (i in 1:nrow(cats.hgr.dset)){
  cats.hgr.dset[i,2] = sum(cats.hgr[which(cats.hgr$Category == cats.hgr.dset[i,1]),3])
}

cats.all = rbind(cats.umgs.dset, cats.hgr.dset)
cats.all$Count = as.numeric(cats.all$Count)

# define cat order
cats.order = data.frame(cbind(cats, rep(0, length(cats))), stringsAsFactors = FALSE)
colnames(cats.order) = c("Category", "Count")

for (i in 1:nrow(cats.order)){
  cats.order[i,2] = sum(cats.all[which(cats.all$Category == cats.order[i,1]),2])
}
cats.order$Count = as.numeric(cats.order$Count)
order = cats.order[order(cats.order$Count, decreasing=TRUE)[1:25],"Category"]

# plot BGC categories
print(ggplot(cats.all, aes(x=Category, y=as.numeric(Count), fill=Genome)) 
      + geom_bar(stat="identity", width = 0.6, alpha=0.5)
      + theme_bw()
      + ylab("Number of BGCs")
      + scale_fill_manual(name="Genome", values=c("steelblue", "red3"))
      + scale_x_discrete(limits=order)
      #+ ggtitle("Chloroflexi")
      + theme(plot.title = element_text(hjust = 0.5))
      + theme(axis.title.x = element_blank())
      + theme(axis.text.y = element_text(size=12))
      + theme(axis.title.y = element_text(size=14))
      + theme(axis.text.x = element_text(face="plain", size=12, angle=45, hjust = 1, vjust=1)))
