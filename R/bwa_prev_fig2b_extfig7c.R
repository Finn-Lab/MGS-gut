# load libraries
library(readxl)
library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(grid)

# load files
bwa.prev = read.csv("bwa_presence-absence.csv", row.names=1) # csv file with presence/absence binary matrix
tax.umgs.hq = read.delim("../../taxonomy/taxonomy_umgs-hq.tab", header=FALSE) # tabular file with taxonomic info
tax.umgs.mq = read.delim("../../taxonomy/taxonomy_umgs-mq.tab", header=FALSE) # tabular file with taxonomic info
tax.id = data.frame(rbind(tax.umgs.hq, tax.umgs.mq))
colnames(tax.id) = c("Genome", "Name", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "MGS")
umgs.genomes = as.vector(tax.id$Genome)

# counts per MGS
prev = as.data.frame(rowSums(bwa.prev))
prev.total = data.frame(Genome=rownames(prev), Counts=as.numeric(prev[,1]))
prev.umgs = prev.total[which(prev.total$Genome %in% umgs.genomes),]
prev.umgs.top = prev.umgs[order(prev.umgs[,2], decreasing=TRUE),]
prev.umgs.top = prev.umgs.top[1:20,]
prev.umgs.top = merge(prev.umgs.top, tax.id, by="Genome")
prev.umgs.top = prev.umgs.top[order(prev.umgs.top$Counts, decreasing=TRUE),]

# plot bargraph of counts per ref genome for iTOL
color_class = c(Actinobacteria= "#09e516", Alphaproteobacteria= "#1d7b37", Bacilli= "#CB1414", Bacteroidia= "#CB7014",
                Betaproteobacteria= "#ffc1f6", Clostridia= "#6992d8", Coriobacteriia= "#C30EAD", Deltaproteobacteria= "#7a4e82",
                Epsilonproteobacteria= "#003ea7", Erysipelotrichia= "#6b552b", Flavobacteriia= "#0ba306", Fusobacteriia= "#C3C00E",
                Gammaproteobacteria= "#9c5999", Mollicutes= "#2100f9", Negativicutes= "#ffc096", Opitutae= "#0b0a0a", Sphingobacteriia= "#909090",
                Spirochaetia= "#775b7b", Synergistia= "#95cde8", Tissierellia= "#FDB103", Verrucomicrobiae= "#fd4b1e")

print(ggplot(prev.umgs.top, aes(x=Genome, y=Counts, fill=Class))
      + geom_bar(stat="identity", color="black", size=0.4, alpha=0.6)
      + scale_x_discrete(limits=prev.umgs.top$Genome, labels=paste(prev.umgs.top$MGS, " (",prev.umgs.top$Name,")",sep=""))
      + scale_fill_manual(values=color_class)
      + theme_bw()
      + ylab("Frequency")
      + theme(axis.text.x = element_text(face="plain", size=10, angle=45, hjust = 1, vjust=1))
      + theme(axis.text.y = element_text(face="plain", size=11))
      + theme(axis.title.y = element_text(size=14))
      + theme(axis.title.x = element_blank()))

# analyse distribution
prev.total = rbind(prev.umgs, prev.hgr)
print(ggplot(prev.total, aes(x=Counts, fill=Type)) 
      + geom_histogram(colour="black", alpha=0.5, size = 0.1, bins=200)
      + geom_vline(xintercept = 31, lty = "dashed", size=0.5)
      + theme_bw()
      + ylab("Number of UMGS")
      + xlab("Number of samples")
      + scale_fill_manual(values=c("darkgreen", "steelblue"))
      + xlim(0,100)
      + theme(axis.title.y = element_text(size=12))
      + theme(axis.text.y = element_text(size=12))
      + theme(axis.title.x = element_text(size=12))
      + theme(axis.text.x = element_text(size=12)))
