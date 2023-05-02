#Katie Falà
#Script to produce boxplots for analysis of Micromatrix plate count data (for Giulia B.)
#02/05/23

library(openxlsx)
library(janitor)
library(dplyr)
library(ggplot2)
library(devEMF)

#################
# read n' clean #
#################
df <- read.xlsx("data\\box_plot_bixon.xlsx")

#split df into 3 separate dataframes - coliforms, LAB, Bifs.
#insert extra column to denote replicate experiment number (fsi1, fsi2, fsi3, fsi4)
df$exp <- c("0", rep(1:4,7))
df$sample <- paste(df$X1,"_",df$exp)
rownames(df) <- df$sample

###############
# coliforms   #
###############
coli <- df %>% select(2:4)
colnames(coli) <- coli[1,]
coli<-coli[-1,]
coli <- t(coli)
coli<-clean_names(coli)
coli<-t(coli)

#reshape to long format tibble
coli |> 
  as.data.frame() |> 
  tibble::rownames_to_column("letter") |> 
  tidyr::pivot_longer(colnames(coli)) |> 
  subset(value != 0)-> colil
colil$value <- as.numeric(colil$value) #ensure counts are numeric
colil$letter <- gsub('.{2}$', '', colil$letter)

#plot
ggplot(colil, aes(x=value, y=letter, fill=letter)) + 
  geom_boxplot() +
  geom_jitter() + 
  theme(axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        legend.position = "none") + 
  xlab("log CFU/mL") +
  ylab("") +
  facet_wrap(~name) + ggtitle("Coliforms") +
  theme(plot.title = element_text(size=20, face="bold"))

ggsave("figures\\coliform.emf", width = 10, height = 10, dpi = 300, device = {function(filename, ...) devEMF::emf(file = filename, ...)})

###############
# LAB   #
###############
lab <- df %>% select(5:7)
colnames(lab) <- lab[1,]
lab<-lab[-1,]
lab <- t(lab)
lab<-clean_names(lab)
lab<-t(lab)

#reshape to long format tibble
lab |> 
  as.data.frame() |> 
  tibble::rownames_to_column("letter") |> 
  tidyr::pivot_longer(colnames(lab)) |> 
  subset(value != 0)-> labl
labl$value <- as.numeric(labl$value) #ensure counts are numeric
labl$letter <- gsub('.{2}$', '', labl$letter) #remove replicate number (fsi_1 etc)

ggplot(labl, aes(x=value, y=letter, fill=letter)) + 
  geom_boxplot() +
  geom_jitter() + 
  theme(axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        legend.position = "none") + 
  xlab("log CFU/mL") +
  ylab("") +
  facet_wrap(~name) + ggtitle("LAB") +
  theme(plot.title = element_text(size=20, face="bold"))

ggsave("figures\\lab.emf", width = 10, height = 10, dpi = 300, device = {function(filename, ...) devEMF::emf(file = filename, ...)})


###############
# bif   #
###############
bif <- df %>% select(8:10)
colnames(bif) <- bif[1,]
bif<-bif[-1,]
bif <- t(bif)
bif<-clean_names(bif)
bif<-t(bif)

#reshape to long format tibble
bif |> 
  as.data.frame() |> 
  tibble::rownames_to_column("letter") |> 
  tidyr::pivot_longer(colnames(bif)) |> 
  subset(value != 0)-> bifl
bifl$value <- as.numeric(bifl$value) #ensure counts are numeric
bifl$letter <- gsub('.{2}$', '', bifl$letter) #remove replicate number (fsi_1 etc)

ggplot(bifl, aes(x=value, y=letter, fill=letter)) + 
  geom_boxplot() +
  geom_jitter() + 
  theme(axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        legend.position = "none") + 
  xlab("log CFU/mL") +
  ylab("") +
  facet_wrap(~name) + ggtitle("Bifidobacteria") +
  theme(plot.title = element_text(size=20, face="bold"))

ggsave("figures\\bif.emf", width = 10, height = 10, dpi = 300, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
