#21.04.22
#This script reads a .csv file formatted as follows:
#Column1: Assembly_name, list of assembly names
#Column2: X.GC, list of %GC values associated with each assembly
#Column3: GTDB.tk, contains species assigned to each assembly as per GTDB.tk
#and generates basic graphs to visualise the data, as well as computing basic stats

library(ggplot2)
library(dplyr)

#Read the csv file, convert GTDB.tk to factor
data <- read.csv("GC_content.csv")
data$GTDB.tk <- as.factor(data$GTDB.tk)

#Boxplots
graph1 <- ggplot (data, aes(x=data$GTDB.tk, y=data$X.GC,color=data$GTDB.tk))+
  geom_boxplot()+
  ggtitle("%GC content of assemblies assigned to each species")+
  xlab("Species")+
  ylab("%GC content")
graph1

#Scatterplots
graph2<- ggplot(data, aes(x=`GTDB.tk`, y=`X.GC`, color=`GTDB.tk`))+
  geom_point()+
  ggtitle("%GC content of assemblies assigned to each species",)+
  xlab("Species")+
  ylab("%GC content")+
  scale_x_discrete(labels = NULL, breaks = NULL)+
  labs(x = "")
graph2

#Stats-lite
#Split the dataframe into a list named l of dataframes, one for each species
#Then recover as a vector 'alvei' the %GC values for each species as a vector
l <- split(data, data$GTDB.tk)
species1 <- l$`species1`$X.GC
species2 <- l$`species2`$X.GC
#etc.....

#Check assumptions of t-test via summary stats
group_by(data, GTDB.tk) %>%
  summarise(
    count = n(),
    mean = mean(X.GC, na.rm = TRUE),
    sd = sd(X.GC, na.rm = TRUE)
  )

t.test(species1, species2, alternative = "two.sided", var.equal = TRUE)
