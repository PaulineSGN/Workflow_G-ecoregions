#Ce script permet de créer un fichier nous disant pour chaque taxon si on a obtenu un modèle BRT
#Ainsi que la liste des taxon a soumettre a worms pour obtenir leur classification et trier les redondance de comptage 

#load packages
library(dplyr)
library(taxonomyCleanr)
library(stringr)
library(xlsx)
 
#load arguments
args = commandArgs(trailingOnly=TRUE) 

if (length(args)==0){
    stop("This tool needs at least one argument")
}else{
    data <- args[1]
    preds <- args[2]
    enviro <- args[3]
}

env = read.table(enviro, header=T, na.strings = "na")
occurrence_files = strsplit(data,",")
preds_files = strsplit(preds,",")

#########functions##########

`%!in%` <- Negate(`%in%`)

have_model = data.frame()
pres = 0

have.model <- function(taxon_phylum,noms_sp,comptage_sp,brt_phylum){
  for (tax in taxon_phylum) {
    if (tax %in% names(noms_sp)){
      pres = sum(comptage_sp[tax])
    }
    if (tax %in% brt_phylum$spe  ) {
      brt = c(tax,"O", pres)
      have_model = rbind(have_model,brt, make.row.names = F)}
    else {
      brt = c(tax,"N", pres)
      have_model = rbind(have_model,brt, make.row.names = F)}
  }
  colnames(have_model) = c("taxon","Model","nb_presence")
  return(have_model)}

##########Execution########
brt = NULL
for (j in 1:length(preds_files[[1]])){
    brt <- rbind(brt,read.table(preds_files[[1]][j], header = TRUE, na.strings = "na"))
}

for (i in 1:length(occurrence_files[[1]])) {
  occurrence <- NULL
  cmpt <- NULL
  taxon <- list()
  
  occurrence <- read.table(occurrence_files[[1]][i], dec = ",", sep = ";", header = TRUE, na.strings = "na")
  
  taxon_names <- names(occurrence)
  new_taxon <- taxon_names[!(taxon_names %in% names(env)) & taxon_names != "station"]
  taxon <- c(taxon, new_taxon)
  
  cmpt <- occurrence[, new_taxon]
  cmpt <- as.data.frame(cmpt)
  
  have_model <- have.model(taxon, occurrence, cmpt, brt)
}


#Taxon pour lesquels on a obtenu un modèle
#have_model2 = subset(have_model, have_model$`Model` != "N")

#Obtention d'une liste de taxon (nettoyé) ayant obtenu un modèle BRT (fichier qui sera soumis à l'outils match taxa de la base de données WoRMS pour obtenir leur classification et pouvoir trier les doublons entre les rangs taxonomique)

tax_net = NULL
for (tax in have_model$taxon){
  tax_net = rbind(tax_net,trim_taxa(tax))
}

have_model$taxon = as.character(tax_net)

have_model

#Supprimer les doublons (ici iophon et iophon sp1.)
#for(tax in have_model2$taxon){
 # if (str_ends(tax,pattern = "sp.1")==T){
  #  have_model2 = have_model2 %>% filter(have_model2$taxon != tax)
  #}
#}

#extraction de l'objet have_model2
#write.xlsx(have_model2,file = "have_model.xlsx",row.names = F)

#obtention de la liste finale à soumettre a worms
#liste_taxon = have_model2$taxon
#write.table(liste_taxon,file = "liste_taxon.txt", quote = F, row.names = F, col.names = F)


