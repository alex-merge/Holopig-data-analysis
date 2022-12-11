## Scripts pour diviser les fichiers TSV en plus petits fichiers pour faciliter
## le traitement.

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Paramètres généraux 
filepath = "raw_data/Quantifications_and_functional_annotations.tsv"
output_basename = "Annot_"
nbr_fichier = 20 #Nombre de fichier voulus


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
data = read.csv(filepath, sep="\t")

quant = c()
for (q in 1:nbr_fichier){
  quant = c(quant, q/nbr_fichier)
}

limits = c(0,as.integer(quantile(1:dim(data)[1], probs = quant)))

for (id in 2:length(limits)){
  out_df = data[limits[id-1]+1:limits[id],]
  write.table(out_df, 
              file = paste("data/",output_basename,id,".csv", sep=""), 
              sep=";")
}
