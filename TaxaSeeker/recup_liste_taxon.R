#Ce script permet de créer un fichier nous disant pour chaque taxon si on a obtenu un modèle BRT
#Ainsi que la liste des taxon 

#load packages
library(dplyr, warn.conflicts = FALSE)
library(taxonomyCleanr, warn.conflicts = FALSE)
library(stringr, warn.conflicts = FALSE)
 
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
      brt = c(tax,"Yes", pres)
      have_model = rbind(have_model,brt, make.row.names = F)}
    else {
      brt = c(tax,"No", pres)
      have_model = rbind(have_model,brt, make.row.names = F)}
  }
  colnames(have_model) = c("Taxa","Model","Occurences")
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
have_model2 = subset(have_model, have_model$`Model` != "N")
have_model3 = subset(have_model, have_model$`Model` != "N")

#Obtention d'une liste de taxon (nettoyé) ayant obtenu un modèle BRT (fichier qui sera soumis à l'outils match taxa de la base de données WoRMS pour obtenir leur classification et pouvoir trier les doublons entre les rangs taxonomique)

have_model2$Taxa <- as.character(trim_taxa(have_model2$Taxa))

#Second nettoyage (élimination de tout les taxon qui finissent par sp1./sp2 etc qui represente un doublon)

have_model2 <- have_model2 %>% filter(!str_ends(Taxa, "sp.1|sp[0-9]"))
have_model3 <- have_model3 %>% filter(!str_ends(Taxa, "sp.1|sp[0-9]"))
have_model <- have_model %>% filter(!str_ends(Taxa, "sp.1|sp[0-9]"))

#extraction de l'objet have_model
write.csv(have_model,file = "have_model.csv", quote = F, row.names = F)

#obtention de la liste pour la suite du workflow si on n'utilise pas worms
list_taxon = have_model3$Taxa
write.table(list_taxon, file= "list_taxa.txt", quote = F, row.names = F, col.names = F)

#obtention de la liste finale à soumettre a worms
liste_taxon = have_model2$Taxa
write.table(liste_taxon,file = "list_taxa_clean.txt", quote = F, row.names = F, col.names = F)







