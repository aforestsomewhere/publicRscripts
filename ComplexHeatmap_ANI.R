#01.05.22
#This script uses the ComplexHeatmap package to generate a series of heatmaps using multiple datasets concerning the same subjects
#Individual heatmaps are constructed for each dataset, they are then joined together and annotated with a column annotation

library(reshape2)
library(ComplexHeatmap)
library(gplots)
library(dplyr)
library(circlize)

######## Generate ht1, which will display the ANI heatmap alone ########

#Read the fastani, convert to matrix
x <- read.table("fastANIoutput_Pared.txt")
ani <- acast(x, V1~V2, value.var="V3")
ani[is.na(ani)] <- 70
#QC on matching names between ani and cardmatrix
#remove .fna from rownames and colnames
row.names(ani) <- gsub('.fna', '', row.names(ani))
colnames(ani) <- gsub('.fna', '', colnames(ani))
#reconcile matrix/df names - change all to _ from -
row.names(ani) <- gsub('-', '_', row.names(ani))
colnames(ani) <- gsub('-', '_', colnames(ani))
#Remove X from colnames
row.names(ani) <- gsub('X', '', row.names(ani))
colnames(ani) <- gsub('X', '', colnames(ani))

#Convert ani dataframe to numeric matrix
ani<- as.matrix(ani)
class(ani) <- "numeric"

#Define colour palette for the ANI heatmap
col_ht1 = colorRamp2(c(100,95,80), c("#F21A00", "#EBCC2A","#78B7C5"))

#Finally, plot the ANI heatmap
ht1 <-Heatmap(ani,
  name = "ANI (%ID)",
  col=col_ht1,
  width = unit(15, "cm"), height = unit(15, "cm"),
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize=5),
  cluster_rows = TRUE,
  column_title = "Average Nucleotide Identity (%)")
draw(ht1)

######## Now generate ht2, containing the CARD results ########

#read card data, adjust rownames, convert NA to 0
y <- read.csv("card_processed2.csv", stringsAsFactors = FALSE)
rownames(y) <- y[,1]
y <- y[,-1]
y[is.na(y)] <- 0

#QC
row.names(y) <- gsub('.fna', '', row.names(y))
colnames(y) <- gsub('.fna', '', colnames(y))
row.names(y) <- gsub('-', '_', row.names(y))
colnames(y) <- gsub('-', '_', colnames(y))
#remove X in genome names
row.names(y) <- gsub('X', '', row.names(y))
colnames(y) <- gsub('X', '', colnames(y))

#plot
y<-t(y)
y<- as.matrix(y)
col_ht2 = colorRamp2(c(0,1,2), c("white", "#F2AD00", "#F98400"))

ht2 <- Heatmap(y, 
               name = "CARD predicted AMR genes", 
               col=col_ht2, 
               row_names_gp = gpar(fontsize = 7),
               column_names_gp = gpar(fontsize=7), 
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               row_dend_reorder = FALSE,
               rect_gp = gpar(col = "white", lwd = 2),
               column_title = "AMR genes (CARD)")
draw(ht2)


######Now generate ht3, containing the PlasmidFinder results ########

p <- read.csv("plasmid_matrix.csv", stringsAsFactors = FALSE)
#fix rownames
rownames(p) <- p[,1]
p <- p[,-1]

#QC on matching names between ani and p
#remove .fna from rownames and colnames
row.names(p) <- gsub('.fna', '', row.names(p))
colnames(p) <- gsub('.fna', '', colnames(p))
#reconcile matrix/df names - change all to _ from -
row.names(p) <- gsub('-', '_', row.names(p))
colnames(p) <- gsub('-', '_', colnames(p))
#Remove X from colnames
row.names(p) <- gsub('X', '', row.names(p))
colnames(p) <- gsub('X', '', colnames(p))

#reorder alphabetically (allows it to cluster properly in next step) and convert to numeric
p<- p[order(rownames(p)) , ,drop=F]
lapply(p,as.numeric)
p<- as.matrix(p)

#Set colours
col_ht3 = colorRamp2(c(0,2,4), c("white", "#00A08A", "#5BBCD6"))

#Plot heatmap
ht3 <- Heatmap(p, 
               col=col_ht3,  
               name="PlasmidFinder predicted replicons", 
               cluster_rows = FALSE, 
               cluster_columns = FALSE,
               row_dend_reorder = FALSE,
               row_names_gp = gpar(fontsize = 7),
               column_names_gp = gpar(fontsize=7),
               rect_gp = gpar(col = "white", lwd = 2),
               column_title = "Replicons (PlasmidFinder)")
draw(ht3)

####### Now generate ht3, containing the VFDB results ########

#Read input file, fix rownames, convert NA to 0
v <- read.csv("vfdb.csv", stringsAsFactors = FALSE)
rownames(v) <- v[,1]
v <- v[,-1]
v[is.na(v)] <- 0

#QC on matching names between ani and p
#remove .fna from rownames and colnames
row.names(v) <- gsub('.fna', '', row.names(v))
colnames(v) <- gsub('.fna', '', colnames(v))
#reconcile matrix/df names - change all to _ from -
row.names(v) <- gsub('-', '_', row.names(v))
colnames(v) <- gsub('-', '_', colnames(v))
#Remove X from colnames
row.names(v) <- gsub('X', '', row.names(v))
colnames(v) <- gsub('X', '', colnames(v))

#alphabeticaly and numeric
v<- v[order(rownames(v)) , ,drop=F]
lapply(v,as.numeric)

#plot
v<- as.matrix(v)
col_ht7 = colorRamp2(c(0,1,2), c("white", "#9986A5", "#79402E"))

ht4 <- Heatmap(v, 
               name = "VFDB", 
               col=col_ht7, 
               row_names_gp = gpar(fontsize = 7),
               column_names_gp = gpar(fontsize=7), 
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               row_dend_reorder = FALSE,
               rect_gp = gpar(col = "white", lwd = 2),
               column_title = "Virulence (VFDB)")
ht4

#Join the individual heatmaps together
ht5 = ht1 + ht2 + ht3 +ht4
draw(ht5)

#####Now to generate row annotations ***************************************************

tax <- read.csv("species_tax.csv", header=FALSE)
#keep just relevant columns on taxonomy
my_metadata <- cbind(tax$V1, tax$V2)
rownames(my_metadata) <- my_metadata[,1]
my_metadata <- my_metadata[,-1]
t <- as.data.frame(my_metadata)

#QC on t
#remove .fna from rownames and colnames
row.names(t) <- gsub('.fna', '', row.names(t))
#reconcile matrix/df names - change all to _ from -
row.names(t) <- gsub('-', '_', row.names(t))
#Remove X from colnames
row.names(t) <- gsub('X', '', row.names(t))

#Generate Species Row annotation
RowAnn <- data.frame(t$my_metadata)
colnames(RowAnn) <- c("Species")

#Read metadata, extract columns of interest, format
o <- read.csv("stats.csv", stringsAsFactors = FALSE)
source <- cbind(o$Assembly_name, o$Source)
rownames(source) <- source[,1]
source <- source[,-1]
# s <- as.data.frame(source)

#QC
#remove .fna from rownames and colnames
row.names(s) <- gsub('.fna', '', row.names(s))
#reconcile matrix/df names - change all to _ from -
row.names(s) <- gsub('-', '_', row.names(s))
#Remove X from colnames
row.names(s) <- gsub('X', '', row.names(s))
s <- as.data.frame(s)

#order both alphabetically by rowname, then merge to new dataframe u, change NA to 'Unknown'
t<- t[order(rownames(t)) , ,drop=F]
s<- s[order(rownames(s)) , ,drop=F]
u<- cbind(t, s$source)
colnames(u) <- c('Species', 'Source')
u[is.na(u)] <- "Unknown"

#Define colour palette
colours <- list("Species"=c("paralvei"="#D8B70A","alvei"="#02401B","americana"="#C7B19C"),
                'Source' = c('Clinical' = "#B40F20", 'Environmental' = "#D69C4E", 'Food'="#ABDDDE", 'Unknown'='gray'))
#Generate the column annotation
colAnno <- HeatmapAnnotation(df = u,
                              which = 'col',
                              col = colours,
                              annotation_width = unit(c(1, 4), 'cm'),
                              gap = unit(1, 'mm'))

#Plot a heatmap with the column annotation
ht6 <-Heatmap (ani,
                 name = "ANI (%ID)",
                 col=col_ht1,
                 width = unit(15, "cm"), height = unit(15, "cm"),
                 row_names_gp = gpar(fontsize = 3),
                 column_names_gp = gpar(fontsize=5),
                 cluster_rows = TRUE,
                 top_annotation = colAnno,
                 column_title = "Average Nucleotide Identity (%)")

#Lastly, 
draw(ht5,heatmap_legend_side="left", annotation_legend_side="left",
     legend_grouping = "original")

draw(ht6+ht2+ht3+ht7, heatmap_legend_side="left", annotation_legend_side="left",
     legend_grouping = "original")
