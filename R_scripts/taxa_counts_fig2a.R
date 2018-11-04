# load libraries
library(reshape2)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gdata)
library(Rmisc)

# set working directory and load large data
dset = read.table("taxacounts_umgs-all.tab", sep="\t") # load UMGS counts per taxon
colnames(dset) = c("Rank", "Taxon", "Counts")
dset$Percentage = dset$Counts/dset[dset$Rank=="Kingdom","Counts"]*100

# order datasets
rank_order = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
dset$Rank = reorder.factor(dset$Rank, new.order=rank_order)
dset = dset[order(dset$Rank,dset$Counts),]
dset$Taxon = factor(dset$Taxon, levels=dset$Taxon)

# label other taxa
ranks_selected = c("Phylum", "Class", "Order", "Family", "Genus")
rm(new.dset)
for (r in 1:length(ranks_selected)){
  prefix = tolower(substr(ranks_selected[r], 1, 1))
  dset.rank = dset[dset$Rank == ranks_selected[r],]
  dset.rank.others = dset.rank[1:(nrow(dset.rank)-5),]
  dset.newrow = data.frame(Rank=ranks_selected[r], Taxon=paste(prefix,"__Other", sep=""), Counts=sum(dset.rank.others$Counts), 
                            Percentage=sum(dset.rank.others$Percentage), Colour="grey")
  dset.rank.top = dset.rank[(nrow(dset.rank)-4):nrow(dset.rank),]
  dset.toprow = data.frame(Rank=ranks_selected[r], Taxon=dset.rank.top$Taxon, Counts=dset.rank.top$Counts, 
                           Percentage=dset.rank.top$Percentage, Colour=rev(brewer.pal(5,"Set1"))) 
  if (!exists("new.dset")) {
    new.dset = rbind(dset.newrow, dset.toprow)
  } else {
    new.dset = rbind(new.dset, dset.newrow, dset.toprow)
  }
}
new.dset$Taxon = factor(new.dset$Taxon, levels=new.dset$Taxon)

# Phylum counts
print(dset[which(dset$Rank == "Phylum"),])

# plot stacked plot taxa counts
print(ggplot(new.dset, aes(x=Rank, y=Percentage, fill=Taxon)) 
      + geom_bar(stat="identity", colour="black", alpha=0.5, size=0.2)
      + theme_bw()
      + ylab("Proportion (%)")
      + scale_fill_manual(values=as.vector(new.dset$Colour))
      + scale_x_discrete(limits=ranks_selected)
      + ylim(0,100)
      + theme(axis.title.y = element_text(size=18))
      + theme(axis.text.y = element_text(size=15))
      + theme(axis.title.x = element_blank())
      + theme(axis.text.x = element_text(size=15, angle=90, hjust=1, vjust=0.3)))
