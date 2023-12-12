setwd('~/Desktop/Projet_AB/A4')
library('RJSONIO')
library(jsonlite)

json_file_path <- "out.json"
json_data <- jsonlite::fromJSON(json_file_path)

data_list <- list()

for (patient_name in names(json_data)) {
  patient_data <- json_data[[patient_name]]
  patient_data$Patient <- sub("Biliary Atresia Patient ", "", patient_name)
  data_list[[length(data_list) + 1]] <- patient_data
}

result_df <- do.call(rbind, data_list)
result_df <- result_df[, c("Patient", names(json_data[[1]])[1:(length(names(json_data[[1]])) - 1)])]
result_df[, 2:ncol(result_df)] <- lapply(result_df[, 2:ncol(result_df)], as.numeric)

patients_a_garder <- c(13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 39, 40, 41, 44, 47)
data_filtered <- result_df[patients_a_garder, ]

death <- c(7.25,7.51,6.13,10.2,11.97,6.85,7.61,5.74,16.59,15,15,7.28,9.28,5.38,8.85,10.5,5.02,6.5,5.51)

library(dplyr)

data_filtered1 <- as.data.frame(data_filtered)
data_filtered2 <- mutate(data_filtered1, death = death)
#supression de la colonne numéro de patient : 
data_filtered3 <- select(data_filtered2, -Patient)
#Inversion du tableau : 
data <- t(data_filtered3)
#Rajouter Score adrénergique et cholinergique et neuronale
library("GSVA")
library("clusterProfiler")

#boucle pour calculer le SA et le SC par patient et rajouter le résultat a la suite 
# Initialiser un vecteur pour stocker les scores ssGSEA
scores_ssgsea <- numeric(ncol(data))

# Boucle for pour calculer le score ssGSEA pour les 10 premières lignes par patient
for (i in seq_along(colnames(data))) {
  patient_data <- data[, i]
  gene_sets <- data[1:21]  # Exclure les gènes du patient en cours
  buffer <- ssgseaParam(exprData = patient_data,
                        geneSets = gene_sets)
  gsva_result <- gsva(buffer, method = "ssgsea")
}

# Ajouter la ligne "Score Cholinergique" au tableau
data <- rbind(data, scores_ssgsea)
rownames(data)[nrow(data)] <- "Score Cholinergique"


















test <- ssgseaParam(exprData = data_filtered,
            geneSets = data_filtered3)

testF <- gsva(test@geneSets[["ENSG00000138435"]])

