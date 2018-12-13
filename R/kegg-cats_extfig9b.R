# load library
library(ALDEx2)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(FactoMineR)
library(CoDaSeq)

# set working directory and read input files
kegg.dset = read.csv("kegg/kegg_summary.csv", row.names=1, check.names=FALSE)
kegg.names = read.delim("kegg/KEGG_orthology_simple.tab", sep="\t", row.names=1, check.names=FALSE, quote="")

# load kegg results per phylum
actino = read.csv("kegg/actino_kegg_results.csv", row.names=1)
actino = merge(actino, kegg.names, by="row.names")
firm = read.csv("kegg/firm_kegg_results.csv", row.names=1)
firm = merge(firm, kegg.names, by="row.names")
bacteroi = read.csv("kegg/bacteroi_kegg_results.csv", row.names=1)
bacteroi = merge(bacteroi, kegg.names, by="row.names")
teneri = read.csv("kegg/teneri_kegg_results.csv", row.names=1)
teneri = merge(teneri, kegg.names, by="row.names")
proteo = read.csv("kegg/proteo_kegg_results.csv", row.names=1)
proteo = merge(proteo, kegg.names, by="row.names")
res.all = rbind(actino, firm, bacteroi, teneri, proteo)

# summarize results and get most significant
effect_thresh = 0.2
res.sig = res.all[which(res.all$wi.eBH<0.05),]
res.eff = res.sig[which(abs(res.sig$effect) > effect_thresh),]
res.plot = res.eff
res.hgr = res.plot[which(res.plot$effect < 0),]
res.umgs = res.plot[which(res.plot$effect > 0),]
res.hgr = data.frame(sort(table(droplevels(res.hgr$`Functional category`))))
res.hgr$Prop = res.hgr$Freq/sum(res.hgr$Freq)*100
res.hgr$Type = "HGR"
res.umgs = data.frame(sort(table(droplevels(res.umgs$`Functional category`))))
res.umgs$Prop = res.umgs$Freq/sum(res.umgs$Freq)*100
res.umgs$Type = "UMGS"
res.cats = rbind(res.hgr, res.umgs)

# plot results
print(ggplot(res.cats, aes(x=Var1, y=Freq, fill=Type))
      + geom_bar(stat="identity", alpha=0.7, size=0.4, alpha=0.8)
      + theme_bw()
      + coord_flip()
      + ylab("Number of enriched genes")
      + theme(axis.title.x = element_text(size=12))
      + theme(axis.text.y = element_text(size=12))
      + theme(axis.text.x = element_text(size=12))
      + theme(axis.title.y = element_blank())
      + scale_fill_manual(name="", values=c("steelblue", "tomato3"))
      + scale_x_discrete(limits=rev(res.umgs$Var1[order(res.umgs$Freq)])[1:10]))
