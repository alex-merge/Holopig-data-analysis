## Scripts pour diviser les fichiers TSV en plus petits fichiers pour faciliter
## le traitement.

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Paramètres généraux 
filepath = "raw_data/Quantifications_and_functional_annotations.tsv"
output_basename = "Annot_"
nbr_fichier = 20 #Nombre de fichier voulus


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#Importation of the .tsv with the header as a line
data = read.csv(filepath, sep="\t", header = FALSE)

#Parting the raw data by the number of lines
nb_rows = nrow(data)
nb_rows_per_file = ceiling(nb_rows/nbr_fichier) #round up the number
header_line = data[1,]

nb_row_inf = 1 #row number of the header line

for (id in 1:nbr_fichier) {
  nb_row_inf = nb_row_inf + 1
  nb_row_sup = id*nb_rows_per_file
  out_df = rbind(header_line, data[nb_row_inf:nb_row_sup, ])
  write.table(out_df, 
              file = paste("refined_data/",output_basename,id,".csv", sep=""), 
              row.names = FALSE, col.names = FALSE, sep=";")
  nb_row_inf = nb_row_sup
}

