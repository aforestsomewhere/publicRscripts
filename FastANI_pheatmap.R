#20.04.22 FastANI graph
#This script generates a basic plot of FastANI values (pairwise comparison within dataset) using the pheatmap R package
#It also includes section to annotate by isolate source (or other available metadata)
#Main input is the fastANI output file (.txt)
#Additionally you can provide a metadata file ('stats.csv') and pull out metadata of interest with which to annotate the heatmap.

library(pheatmap)
library(dplyr)
library(reshape2)
library(gplots)

#get fastani, convert to matrix
x <- read.table("fastANIoutput_pared.txt")
matrix <- acast(x, V1~V2, value.var="V3")
matrix[is.na(matrix)] <- 70

#remove .fna from rownames and colnames, if present
row.names(matrix) <- gsub('.fna', '', row.names(matrix))
colnames(matrix) <- gsub('.fna', '', colnames(matrix))

### define the colors within 2 zones
#breaks = seq(min(matrix), max(100), length.out=100)
#gradient1 = colorpanel( sum( breaks[-1]<=95 ), "red", "white" )
#gradient2 = colorpanel( sum( breaks[-1]>95 & breaks[-1]<=100), "white", "blue" )

## Define colour breaks, 98% and 95% cutoffs
breaks = seq(min(matrix), max(100), length.out=100)
gradient1 = colorpanel( sum( breaks[-1]<=95), "red", "white" )
gradient2 = colorpanel( sum( breaks[-1]>95 & breaks[-1]<=98), "yellow", "green" )
gradient3 = colorpanel( sum( breaks[-1]>98 & breaks[-1]<=100), "white", "blue" )

#Generate colour palette
hm.colors = c(gradient1, gradient2)

#Plot the base heatmap
pheatmap(matrix, scale = "none", cexRow=.30, cexCol=.30, fontsize=5)

#add basic annotations
metadata <- read.csv("stats.csv")
my_metadata <- cbind(metadata$Assembly_name, metadata$Source)
rownames(my_metadata) <- my_metadata[,1]
my_metadata <- my_metadata[,-1]
df<- as.data.frame(my_metadata)
#If NA present (in source attribution) modify to 'unknown'
df[is.na(df)] <- 'unknown'

#reconcile matrix/df names - change all to _ from -
row.names(df) <- gsub('-', '_', row.names(df))
colnames(df) <- gsub('-', '_', colnames(df))
row.names(matrix) <- gsub('-', '_', row.names(matrix))
colnames(matrix) <- gsub('-', '_', colnames(matrix))


# add group annotation
anno<-data.frame(row.names=rownames(df), Group=df$df)
annoCol<-list(Group=c(Clinical = "#000000", Environmental="#009E73", Food="#F0E442", unknown="Gray"))

#aka3 = list( c(Clinical = "#a65852", Environmental="#343c24", Food="#c4aa23", Unknown="Gray"))
pheatmap(matrix, scale = "none", annotation_col = anno, annotation_colors = annoCol, cexRow=.30, cexCol=.30, fontsize=3)
