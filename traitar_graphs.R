library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)

#Read the Traitar output file
data <- read.csv("predictions_majority-vote_combined.csv", header=TRUE)
#Make genome names rownames
rownames(data) <- data[,2]
data <- data[,-2]
#remove last row (error)
data <- data[-83,]

#remove source info for now
myMat <- data[,-1]
myMat <- as.matrix(sapply(myMat, as.numeric))
#plot plain heatmap
pheatmap(myMat)

#add basic annotations from data

my_metadata <- cbind(rownames(data), data$X)
rownames(my_metadata) <- my_metadata[,1]
my_metadata <- my_metadata[,-1]
my_metadata <- as.data.frame(my_metadata)
my_metadata$my_metadata <- as.factor(my_metadata$my_metadata)

#Make source a factor
my_metadata$X <- as.factor(my_metadata$my_metadata)

# add group annotation
anno<-data.frame(row.names=rownames(my_metadata), Group=my_metadata$my_metadata)
annoCol<-list(Group=c(Clinical = "#000000", Environmental="#009E73", Food="#F0E442", Undocumented="Gray"))

#aka3 = list( c(Clinical = "#a65852", Environmental="#343c24", Food="#c4aa23", Undocumented="Gray"))
#transpose myMat
myMat<-t(myMat)
myMat<- as.matrix(myMat)
colnames(myMat) <- rownames(data)

pheatmap(myMat, scale = "none", annotation_col = anno, annotation_colors = annoCol,fontsize=3)

#############################################################################################

#Collapse all into same category
# #Read the Traitar output file
data <- read.csv("predictions_majority-vote_combined.csv", header=TRUE)
# #Make genome names rownames
rownames(data) <- data[,2]
data <- data[,-2]
# #remove last row (error)
data <- data[-83,]


#Make Source (X) a factor, then recursively split into individual tables
data$X <- as.factor(data$X)
for (i in levels(data$X)) {
  command <- paste0(i, "<-subset(data, X=='", i, "')")
  eval(parse(text=command))
}

################################################################################
#Total for Clinical
#Remove factor column
Clinical <- Clinical[,-1]

df <- mutate_all(Clinical, function(x) as.numeric(as.character(x)))

# Clinicalm <- as.data.frame.numeric(Clinical)
# Clinicalm <- as.matrix(Clinical)
# Clinucaln <- as.data.frame.numeric(Clinicalm)
# Clinicaln <- lapply(Clinicalm,as.numeric)
# remove(Clinicaln)
#Clinicalm <- as.numeric(Clinicalm)
df["Total" ,] <- colSums(df)

#Summary for Clinical
Clinical %>%
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.character), ~"Total")))

################################################################################
#Total for Environmental
Environmental <- Environmental[,-1]
df2 <- mutate_all(Environmental,function(x) as.numeric(as.character(x)))
df2["Total" ,] <- colSums(df2)


################################################################################
#Total for Food
Food <- Food[,-1]
df3 <- mutate_all(Food,function(x) as.numeric(as.character(x)))
df3["Total" ,] <- colSums(df3)

###############################################################################
#Combine the totals rows
totals <- rbind(df["Total",], df2["Total",], df3["Total",])
rownames(totals) <- c("Clinical", "Environmental", "Food")

######Vertical heatmap, scale within each column to normalise across sources

totals<-t(totals)
totals<- as.matrix(totals)
pheatmap(totals, scale="column")

######todo, correct for number of observations
