# load libraries
library(readxl)
library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# load data
setwd("~/Documents/ESPOD/Analyses/Assemb_Binning/MetaSpecies_revision/read_mapping/sourmash/")
sm.dset = read.csv("sourmash_classification.csv", check.names=FALSE, row.names=1)
metadata.raw = read_excel("../../tables/SuppInfo_metadata.xlsx")
metadata = as.data.frame(metadata.raw[,c("pub_state", "pub_disease", "pub_disease_secondary", "pub_agestrat", 
                                         "pub_antibio", "country", "continent")])
rownames(metadata) = metadata.raw$run_accession
colnames(metadata) = c("Disease_state", "Disease_name", "Disease_secondary", "Age", 
                       "Antibio", "Country", "Continent")

# format datasets
dset = merge(sm.dset, metadata, by="row.names")
rownames(dset) = dset$Row.names
dset = dset[,c("HR", "+ RefSeq", "+ UMGS", "Continent")]
dset$Improv = (dset$`+ UMGS`-dset$`+ RefSeq`)/dset$`+ RefSeq`*100

# plot boxplot of read assignment %
dset.box.all = melt(dset)

plot = print(ggplot(dset.box.all, aes(x=variable, y=value, fill=variable)) # overall
#plot = print(ggplot(dset, aes(x=Continent, y=Improv, fill=Continent)) # by continent
             + geom_boxplot(alpha=0.5, width=0.4, outlier.colour=NA, colour="black")
             #+ geom_point(position = position_jitter(width=0.2, height=0), alpha=0.5, size=0.3)
             + theme_bw()
             #+ ylim(0,100)
             #+ ylim(0,750)
             + scale_x_discrete(limits=rev(c("HR", "+ RefSeq", "+ UMGS")))
             + scale_fill_manual(values=c("skyblue1", "skyblue2", "skyblue4"))
             #+ scale_x_discrete(limits=rev(c("North America","Asia", "Europe", "Oceania", "South America", "Africa")))
             #+ scale_fill_manual(limits=rev(c("North America","Asia", "Europe", "Oceania", "South America", "Africa")),
            #                    values=c("salmon", "orchid4", "darkorange", "steelblue", "red3", "green4"))
             + coord_flip()
             + guides(fill=FALSE)
             + ylab("Read classification (%)")
             + ylab("Read classification increase (%)")
             + theme(axis.title.y = element_blank())
             + theme(axis.text.y = element_text(size=12))
             + theme(axis.title.x = element_text(size=14))
             + theme(axis.text.x = element_text(size=12)))

# boxplot data
ggplot_build(plot)$data
