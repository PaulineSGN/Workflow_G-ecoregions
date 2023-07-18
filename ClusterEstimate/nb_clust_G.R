#Date : 21/03/2023
# Script permettant de determiner le nombres de clusters optimal grâce à l'optimisation de l'indice SIH 
# et de produire les fichiers dont on a besoin à l'étape suivante la clusterisation

#load packages
library(cluster)
library(dplyr)
library(tidyverse)

#load arguments
args = commandArgs(trailingOnly=TRUE) 
if (length(args)==0)
{
    stop("This tool needs at least one argument")
}else{
    enviro <- args[1]
    taxa_list <- args[2]
    preds <- args[3]
    max_k <- as.numeric(args[4])
    metric <- args[5]
    sample <- as.numeric(args[6])
}

#load data 

env.data <- read.table(enviro, header = TRUE, dec = ".", na.strings = "-9999.00") 

## liste de taxons modélisés retenus pour la clusterisation
tv <- read.table(taxa_list, dec=".", sep=" ", header=F, na.strings = "NA") 
names(tv) <- c("a")

################Groupement des taxons si plusieurs fichier de prediction entrés ################

data_split = str_split(preds,",")
data.bio = NULL

for (i in 1:length(data_split[[1]])) {
data.bio1 <- read.table(data_split[[1]][i], dec=".", sep=" ", header=T, na.strings = "NA")
data.bio <- rbind(data.bio,data.bio1)
remove(data.bio1)
}

names(data.bio) <- c("lat", "long", "pred", "taxon")

#keep selected taxa
data.bio <- data.bio[which(data.bio$taxon %in% tv$a),]

write.table(data.bio,file="data_bio.tsv",sep="\t",quote=F,row.names=F)

#format data

test3 <- matrix(data.bio$pred , nrow = nrow(env.data),  ncol = nrow(data.bio)/nrow(env.data))
test3 <- data.frame(test3)
names(test3) <- unique(data.bio$taxon)

write.table(test3, file="data_to_clus.tsv", sep="\t",quote=F,row.names=F)

#Nombre de clusters max a tester
max_k <- max_k

# Initialisation des vecteurs pour stocker les indices SIH
sih_values <- rep(0, max_k)

# Calcul de l'indice SIH pour chaque nombre de clusters
for (k in 2:max_k) {
  # Exécution de clara
  clara_res <- clara(test3, k,  metric =metric,  samples = sample, sampsize = min(nrow(test3), (nrow(data.bio)/nrow(test3))+2*k))
  # Calcul de l'indice SIH
  sih_values[k] <- clara_res$silinfo$avg.width
}

# Tracé du graphique de l'indice SIH en fonction du nombre de clusters
png("Indices_SIH.png")
plot(2:max_k, sih_values[2:max_k], type = "b", xlab = "Nombre de clusters", ylab = "Indice SIH")
dev.off()
