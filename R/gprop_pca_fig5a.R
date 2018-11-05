# load libraries
library(RColorBrewer)
library(FactoMineR)
library(ggplot2)
library(vegan)

# load summary files
file_list = Sys.glob("*gp.tab")
gprop.annot = read.delim(file_list[1], sep="\t", header=FALSE)[,1:2]

# load and prepare taxonomy file
tax = read.delim("taxonomy_all.tab", header = FALSE, row.names=1)
colnames(tax) = c("Name", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Genome")

# make sure gprop object is not in environment
if (exists("gprop")){
  rm(gprop)
}

# create dataset
n = 0
for (file in file_list){
  n = n + 1
  time = strptime(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(time$h,":",time$min,":",time$s," Loading: sample ",n, " (",round(n/length(file_list)*100,1),"%)\n",sep="")
  # if the merged dataset doesn't exist, create it
  if (!exists("gprop")){
    gprop = read.delim(file, header=FALSE, stringsAsFactors=FALSE)
    gprop$V3[gprop$V3 == "YES"] = 2
    gprop$V3[gprop$V3 == "PARTIAL"] = 1
    gprop$V3[gprop$V3 == "NO"] = 0
    str = gsub("_gp.tab","",file)
    colnames(gprop) = c("GP", "Function", str)
  }
  # if the merged dataset does exist, append to it
  else {
    gprop_tmp = read.delim(file, header=FALSE, stringsAsFactors=FALSE)
    gprop_tmp$V3[gprop_tmp$V3 == "YES"] = 2
    gprop_tmp$V3[gprop_tmp$V3 == "PARTIAL"] = 1
    gprop_tmp$V3[gprop_tmp$V3 == "NO"] = 0
    gprop = cbind(gprop, gprop_tmp[,3])
    str = gsub("_gp.tab","",file)
    colnames(gprop) = c(colnames(gprop)[-ncol(gprop)], str)
    rm(gprop_tmp)
  }
}

# format dataset
gprop_trim = gprop[,-c(1,2)]
gprop_plot = sapply(gprop_trim, function(x) as.numeric(as.character(x)))
rownames(gprop_plot) = gprop$GP

# prepare dataset for pca plot
dset.plots = t(gprop_plot)
dset.tax = merge(dset.plots, tax, by="row.names")
dset.tax = dset.tax[order(dset.tax$Genome),]
dset.tax$Phylum = factor(dset.tax$Phylum, levels = c("Actinobacteria", "Bacteroidetes", "Cyanobacteria", 
                                                     "Firmicutes", "Fusobacteria", "Proteobacteria", 
                                                     "Saccharibacteria", "Spirochaetes", "Synergistetes", 
                                                     "Tenericutes", "Verrucomicrobia"))

# colours for categorical variables
phy_type = c("#004cce", "#059633", "#bdf2c1", "#CB1414", "#C3C00E", 
             "#9c5999", "#ffc49c", "#7a4b03", "#95cde8", "#e59735", "#ababab") # phylum
gen_type = c("steelblue", "red3", "orange") # genome type

# principal component analysis
dset.plot = dset.tax
pca = PCA(dset.plot[,2:1200], scale.unit = FALSE, graph=FALSE)
pcs = as.data.frame(pca$ind$coord)
pc1_var = round(pca$eig[1,2],1)
pc2_var = round(pca$eig[2,2],1)
gg.dset = data.frame(Sample=dset.plot$Row.names, PC1=pcs$Dim.1, PC2=pcs$Dim.2, 
                     Genome=dset.plot$Genome, Phylum=dset.plot$Phylum)
print(ggplot(gg.dset, aes(x=PC1, y=PC2, colour=Phylum, shape=Genome))
      + geom_point(size=1, alpha=0.8)
      + scale_color_manual(name="Phylum", values=phy_type, na.value="black")
      + scale_shape_manual(name="Genome", values=c(0,4,16))
      + guides(colour = guide_legend(override.aes = list(size=5, alpha=0.8)))
      + theme_bw()
      + guides(shape=FALSE, colour=FALSE)
      + ylab(paste("PC2"," (",pc2_var,"%)",sep=""))
      + xlab(paste("PC1"," (",pc1_var,"%)",sep=""))
      + theme(axis.title.y = element_text(size=14))
      + theme(axis.text.y = element_text(size=12))
      + theme(axis.title.x = element_text(size=14))
      + theme(axis.text.x = element_text(size=12)))

# perform anosim test by phylum
dset.anosim = dset.tax[which(!is.na(dset.tax$Phylum)),]
rownames(dset.anosim) = dset.anosim$Row.names
dset.anosim = dset.anosim[,2:1200]
res.phyl = anosim(dset.anosim, dset.tax$Phylum[which(!is.na(dset.tax$Phylum))], distance="gower")
print(res.phyl)
