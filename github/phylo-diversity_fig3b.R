# load libraries
library(phytools)
library(ggplot2)

# load files
setwd("~/Documents/ESPOD/Analyses/Project_UMGS/MetaSpecies_revision/phylogeny/")
tre = read.tree("hgr-umgs_phylo/raxml_hgr-umgs_phylogeny.nwk") # load tree
tax.hgr = read.delim("../taxonomy/taxonomy_hgr.tab", header=FALSE) # load taxonomy file
tax.umgs = read.delim("../taxonomy/taxonomy_umgs-all.tab", header=FALSE)[1:8] # load taxonomy file
colnames(tax.hgr) = c("Genome", "Name", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
colnames(tax.umgs) = c("Genome", "Name", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
tax = rbind(tax.hgr, tax.umgs)

# prepare file
all.phyla = c(as.vector(tax.hgr$Phylum), as.vector(tax.umgs$Phylum))
phyla = unique(all.phyla[which(!is.na(all.phyla))])
pd = data.frame(matrix(0, ncol=3, nrow=length(phyla)+1))
colnames(pd) = c("Phylum", "HGR_PD", "HGR+UMGS_PD")

# calculate phylogenetic diversity (PD)
for (n in 1:length(phyla)){
  phylum = as.character(phyla[n])
  species.hgr = as.vector(tax.hgr[which(tax.hgr$Phylum == phylum),"Genome"])
  if (length(species.hgr) > 0){
    subs.tre.hgr = drop.tip(tre, tre$tip.label[-match(species.hgr, tre$tip.label)])
    pd[n,2] = sum(subs.tre.hgr$edge.length)
  } 
  species.all = as.vector(tax[which(tax$Phylum == phylum),"Genome"])
  subs.tre.all = drop.tip(tre, tre$tip.label[-match(species.all, tre$tip.label)])
  pd[n,1] = phylum
  pd[n,3] = sum(subs.tre.all$edge.length)
}
tre.hgr = drop.tip(tre, tre$tip.label[-match(tax.hgr$Genome, tre$tip.label)])
pd[length(phyla)+1,] = c("Total", sum(tre.hgr$edge.length), sum(tre$edge.length))
pd$Proportion = ((as.numeric(pd$`HGR+UMGS_PD`)-as.numeric(pd$HGR_PD))
                  /as.numeric(pd$`HGR+UMGS_PD`)*100) # proportion
pd$Improvement = (as.numeric(pd$`HGR+UMGS_PD`)-as.numeric(pd$HGR_PD)) # difference

# reorder phylum
pd.fi = pd[order(pd$Improvement[-nrow(pd)], decreasing=TRUE),]
pd.fi$Phylum = factor(pd.fi$Phylum, levels = pd.fi$Phylum)

# prepare categorical variable
col_type = c(Actinobacteria="#004cce", Bacteroidetes="#059633", Cyanobacteria="#bdf2c1", 
                  Firmicutes="#CB1414", Fusobacteria="#C3C00E", Proteobacteria="#9c5999", 
                  Saccharibacteria="#ffc49c", Spirochaetes="#7a4b03", Synergistetes="#95cde8", 
                  Tenericutes="#e59735", Verrucomicrobia="#ababab") # phylum

# plot phylogenetic diversity
print(ggplot(pd.fi, aes(x=Phylum, y=Improvement, fill=Phylum))       
      + geom_bar(stat="identity", alpha=0.7, size=0.2)
      + theme_bw()
      + scale_fill_manual(values=col_type)
      #+ scale_x_discrete(limits=rev(levels(prop.freq$Var1)))
      + ylab("Increase in phylogenetic diversity\nprovided by the UMGS (total branch length)")
      #+ ylab("Proportion of the total phylogenetic diversity\nprovided by the UMGS (%)")
      + guides(fill=FALSE)
      + coord_flip()
      #+ scale_y_reverse()
      + theme(axis.title.y = element_blank())
      + theme(axis.text.y = element_text(size=12))
      + theme(axis.title.x = element_text(size=14))
      #+ theme(axis.text.x = element_text(size=15))
      #+ theme(axis.title.x = element_text(size=18))
      + theme(axis.text.x = element_text(size=12)))