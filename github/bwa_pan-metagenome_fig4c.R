# load libraries
library(readxl)
library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(FactoMineR)
library(pheatmap)
library(grid)
library(stats)

# load workspace
setwd("~/Documents/ESPOD/Analyses/Assemb_Binning/MetaSpecies_revision/read_mapping/bwa/")

# load files
bwa.prev = read.csv("bwa_presence-absence.csv", row.names=1)
metadata.raw = read_excel("../../tables/SuppInfo_metadata.xlsx")
metadata = as.data.frame(metadata.raw[,c("pub_state", "pub_disease", "pub_disease_secondary", "pub_agestrat", 
                                         "pub_antibio", "country", "continent")])
rownames(metadata) = metadata.raw$run_accession
colnames(metadata) = c("disease_state", "disease_name", "disease_secondary", "age", 
                       "antibio", "country", "continent")
umgs.hq = scan("../../taxonomy/umgs-hq_genomes.txt", what="")
umgs.mq = scan("../../taxonomy/umgs-mq_genomes.txt", what="")
umgs = c(umgs.hq, umgs.mq)
select = umgs

# prepare dataset
bwa.subs = bwa.prev[select,] # select species

locations = c("Asia", "Africa", "Europe", "North America", "South America", "Oceania")
for (c in locations) {
  samp.conds = rownames(metadata[which(metadata$continent == c),])
  samp.cond.analys = bwa.subs[,samp.conds] # select samples
  # subsample in X increments (50 for Asia/North America and Europe, 10 for the rest)
  if (c == "Asia" | c == "North America" | c == "Europe"){
    subsampling = seq(round(ncol(samp.cond.analys)/50),ncol(samp.cond.analys),
                      round(ncol(samp.cond.analys)/50))
  } else {
    subsampling = seq(round(ncol(samp.cond.analys)/10),ncol(samp.cond.analys),
                      round(ncol(samp.cond.analys)/10))
  }
  dset.subsamp = data.frame(matrix(0,length(subsampling),4))
  colnames(dset.subsamp) = c("Samples", "Species", "SD", "Origin")
  dset.subsamp$Origin = c
  
  reps = 10
  for (i in 1:length(subsampling)){
    rep = c(rep(0,reps))
    for (n in 1:reps){
      cols = sample(colnames(samp.cond.analys),subsampling[i])
      rep[n] = length(which(rowSums(samp.cond.analys[,cols]) > 0))
    }
    dset.subsamp[i,1] = subsampling[i]
    dset.subsamp[i,2] = mean(rep)
    dset.subsamp[i,3] = sd(rep)
  }
  if (c == "Asia") { dset.asia = dset.subsamp }
  else if (c == "Europe") { dset.europe = dset.subsamp }
  else if (c == "North America") { dset.namerica = dset.subsamp }
  else if (c == "South America") { dset.samerica = dset.subsamp }
  else if (c == "Africa") { dset.africa = dset.subsamp }
  else if (c == "Oceania") { dset.oceani = dset.subsamp }
}

dset.line = rbind(dset.asia, dset.europe, dset.namerica, dset.africa, dset.samerica, dset.oceani) # all
dset.nea = dset.line[which(dset.line$Origin != "North America" & dset.line$Origin != "Europe" 
                           & dset.line$Origin != "Asia"),] # just NEA

# plot linegraphs with error bars
print(ggplot(dset.nea, aes(x=Samples, y=Species, colour=Origin))
      #+ geom_errorbar(aes(ymin=Species-SD, ymax=Species+SD), size=0.1, width=0.5)
      #+ geom_line(linetype="dotted", size=1, alpha = 0.5)
      + geom_smooth(method="nls", formula=y~SSasymp(x, Asym, R0, lrc), se=F, size=0.35)
      + geom_point(size=0.2)
      + theme_bw()
      + scale_color_manual(values=c("red3", "salmon", "steelblue", "green4", "orchid4", "darkorange"), 
                           limits=c("Asia", "Africa", "Europe", "North America", "South America",
                                    "Oceania"))
      + ylab("Number of UMGS detected")
      + xlab("Number of samples")
      #+ scale_x_continuous(breaks=c(0,200,400,600,800,1000,1200,1400))
      #+ scale_y_continuous(breaks=c(0,5000,10000,15000,20000,25000,30000))
      + theme(legend.title = element_blank())
      + theme(legend.text = element_text(size = 10))
      + theme(axis.title.y = element_text(size=14))
      + theme(axis.text.y = element_text(size=12))
      + theme(axis.title.x = element_text(size=14))
      + theme(axis.text.x = element_text(size=12)))

# non-linear regression (asymptotic)
nls.dset = dset.namerica
print(summary(nls(nls.dset$Species ~ SSasymp(nls.dset$Samples, Asym, R0, lrc), data = nls.dset)))
