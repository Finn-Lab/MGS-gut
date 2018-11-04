# load library
library(ALDEx2)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(FactoMineR)
library(CoDaSeq)

# set working directory and read input files
setwd("~/Documents/ESPOD/Analyses/Assemb_Binning/MetaSpecies_revision/functions/")
#go.slim = read.delim("go-terms/go-slim_summary.tab", sep="\t", row.names=1, check.names=FALSE)
#go.full = read.delim("go-terms/go-full_summary.tab", sep="\t", row.names=1, check.names=FALSE)
#pfam = read.delim("pfam/pfam_summary.tab", sep="\t", row.names=1, check.names=FALSE)
kegg.dset = read.csv("kegg/kegg_summary.csv", row.names=1, check.names=FALSE)
kegg.names = read.delim("kegg/KEGG_orthology_simple.tab", sep="\t", row.names=1, check.names=FALSE, quote="")

# parse dset
#go.names = go.slim[,1:2]
#go.dset = go.slim[,3:ncol(go.slim)]

#go.names = go.full[,1:2]
#go.dset = go.full[,3:ncol(go.full)]

#pfam.names = data.frame(pfam[,1], row.names=rownames(pfam))
#colnames(pfam.names) = "Annotation"
#pfam.dset = pfam[,2:ncol(pfam)]

# load and prepare taxonomy file
tax.hgr = read.delim("../taxonomy/taxonomy_hgr.tab", header = FALSE, row.names=1)
tax.hgr$Genome = "HGR" 
tax.mgs = read.delim("../taxonomy/taxonomy_umgs-hq.tab", header = FALSE, row.names=1)
tax.mgs = tax.mgs[,1:7]
tax.qs50 = read.delim("../taxonomy/taxonomy_umgs-mq.tab", header = FALSE, row.names=1)
tax.qs50 = tax.qs50[,1:7]
tax.umgs = rbind(tax.mgs, tax.qs50)
tax.umgs$Genome = "UMGS"
tax = rbind(tax.hgr, tax.umgs)
colnames(tax) = c("Name", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Genome")

# define  genomes to analyse
phylum = "Tenericutes"
phy.selected = tax[which(tax$Phylum == phylum),]
hgr.genomes = rownames(phy.selected)[which(phy.selected$Genome == "HGR")]
umgs.genomes = rownames(phy.selected)[which(phy.selected$Genome == "UMGS")]

# prepare dset
#analy.dset = go.dset[,c(hgr.genomes, umgs.genomes)] # go
#analy.dset = pfam.dset[,c(hgr.genomes, umgs.genomes)] # pfam
analy.dset = kegg.dset[,c(hgr.genomes, umgs.genomes)] # kegg
analy.dset[is.na(analy.dset)] = 0 # replace NAs with 0
analy.conds = c(rep("HGR",length(hgr.genomes)),rep("UMGS", length(umgs.genomes)))

# perform aldex analysis
aldex.analy = aldex.clr(reads=analy.dset, conds=analy.conds, mc.samples=128)
aldex.eff = aldex.effect(aldex.analy, analy.conds, useMC=TRUE)
aldex.res = aldex.ttest(aldex.analy, analy.conds)
res.all = data.frame(rownames(aldex.eff), aldex.eff,aldex.res)
write.csv(res.all, file="kegg/teneri_kegg_results.csv", row.names=FALSE, quote=FALSE)
#write.csv(res.all, file="go-terms/teneri_go-slim_results.csv", row.names=FALSE, quote=FALSE)
#res.all = read.csv("go-terms/teneri_go-slim_results.csv", row.names=1)

# summarize results and get most significant
effect_thresh = 0.2
res.sig = res.all[which(res.all$wi.eBH<0.05),]
res.eff = res.sig[which(abs(res.sig$effect) > effect_thresh),]
res.plot = merge(res.eff, kegg.names, by="row.names") # kegg
#res.plot = merge(res.eff, go.names, by="row.names") # go terms
#res.plot = merge(res.eff, pfam.names, by="row.names") # pfam

# plot barplot with effect size
res.top = res.plot[order(res.plot$effect, decreasing=TRUE)[1:5],]
res.bot = res.plot[order(res.plot$effect, decreasing=FALSE)[1:5],]
res.best = rbind(res.top, res.bot)
#res.best = res.bot
res.fi = res.best[order(res.best$effect, decreasing=TRUE),]
#print(ggplot(res.fi, aes(x=Product, y=effect, fill=`Functional category`)) # kegg
#print(ggplot(res.fi, aes(x=Annotation, y=effect, fill=effect)) # pfam
print(ggplot(res.fi, aes(x=Function_specific, y=effect, fill=Function_general)) # go
      + geom_bar(stat="identity", alpha=0.8, size=0.4, alpha=0.8)
      + theme_bw()
      + coord_flip()
      #+ ylim(-0.9, 0.9)
      + theme(plot.title = element_text(hjust = 0.5))
      + ylab("Effect size")
      #+ scale_fill_manual(name="", values=brewer.pal(12,"Set3")) # kegg
      #+ scale_x_discrete(limits=res.fi$Product) # kegg
      #+ scale_x_discrete(limits=res.fi$Annotation) # pfam
      + scale_x_discrete(limits=res.fi$Function_specific) # go
      + scale_fill_manual(name="", values=c("red3", "darkgreen", "steelblue"),
                          labels=c("Biological process", "Cellular component", "Molecular function")) # go
      + theme(axis.title.y = element_blank()))

# plot effect size scatter plots
res.all$label = rep(0,nrow(res.all))
res.all[which(res.all$wi.eBH < 0.05),"label"] = "Significant"
res.all[which(res.all$wi.eBH < 0.05 & abs(res.all$effect) > 0.5),"label"] = "Strong effect"
res.all[which(res.all$label == 0), "label"] = "Not significant"
print(ggplot(res.all, aes(x=diff.win, y=diff.btw, colour=label))
      + geom_point(alpha=0.5, size=0.4)
      + geom_abline(intercept = 0, slope = effect_thresh, linetype=3, colour="steelblue", size=1)
      + geom_abline(intercept = 0, slope = -effect_thresh, linetype=3, colour="steelblue", size=1)
      + theme_bw()
      + guides(colour=FALSE)
      #+ ylim(-25,20)
      + scale_colour_manual(values=c("black", "red", "steelblue"))
      + ylab("Median Log2 Difference")
      + xlab("Median Log2 Dispersion"))