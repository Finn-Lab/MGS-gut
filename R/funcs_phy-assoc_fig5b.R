# load library
library(ALDEx2)
library(ggplot2)
library(RColorBrewer)

# load input
go.slim = read.delim("go-slim_summary.tab", sep="\t", row.names=1, check.names=FALSE) # GO slim summary counts
go.names = go.slim[,1:2]
go.dset = go.slim[,3:ncol(go.slim)]
tax.hgr = read.delim("taxonomy_hgr.tab", row.names=1) # taxonomy file of HGR species
tax.umgs = read.delim("taxonomy_umgs.tab", row.names=1) # taxonomy file of UMGS species
rownames(tax.umgs) = tax.umgs$MAG_ID

# define genomes to analyse
phylum = "Firmicutes"
hgr.genomes = rownames(tax.hgr)[which(tax.hgr$Phylum == phylum)]
umgs.genomes = rownames(tax.umgs)[which(tax.umgs$Phylum == phylum)]

# prepare dataset
analy.dset = go.dset[,c(hgr.genomes, umgs.genomes)]
analy.dset[is.na(analy.dset)] = 0 # replace NAs with 0
analy.conds = c(rep("HGR",length(hgr.genomes)),rep("UMGS", length(umgs.genomes)))

# perform aldex analysis
aldex.analy = aldex.clr(reads=analy.dset, conds=analy.conds, mc.samples=128)
aldex.eff = aldex.effect(aldex.analy, analy.conds, useMC=TRUE)
aldex.res = aldex.ttest(aldex.analy)
res.all = data.frame(rownames(aldex.eff), aldex.eff,aldex.res)

# summarize results and get most significant
effect_thresh = 0.2
res.sig = res.all[which(res.all$wi.eBH<0.05),]
res.eff = res.sig[which(abs(res.sig$effect) > effect_thresh),]
res.plot = merge(res.eff, go.names, by="row.names") # go terms

# plot barcharts with effect size
res.top = res.plot[order(res.plot$effect, decreasing=TRUE)[1:5],]
res.bot = res.plot[order(res.plot$effect, decreasing=FALSE)[1:5],]
res.best = rbind(res.top, res.bot)
res.fi = res.best[order(res.best$effect, decreasing=TRUE),]
print(ggplot(res.fi, aes(x=Function_specific, y=effect, fill=Function_general))
      + geom_bar(stat="identity", alpha=0.8, size=0.4, alpha=0.8)
      + theme_bw()
      + coord_flip()
      + theme(plot.title = element_text(hjust = 0.5))
      + ylab("Effect size")
      + scale_x_discrete(limits=res.fi$Function_specific)
      + scale_fill_manual(name="", values=c("red3", "darkgreen", "steelblue"),
                          labels=c("Biological process", "Cellular component", "Molecular function"))
      + theme(axis.title.y = element_blank()))
