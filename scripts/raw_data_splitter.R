## Script for dividing TSV file into smaller CSV in order to make the post treatment
## easier

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# General parameters (can be modified by the user)

filepath        = "raw_data/Quantifications_and_functional_annotations.tsv"
output_basename = "Annot_"
nb_file         = 20       # Number of file wanted at the end


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#Importation of the .tsv with the header as a line
data = read.csv(filepath, sep="\t", header = FALSE)

#Parting the raw data by the number of lines
nb_rows = nrow(data)
nb_rows_per_file = ceiling(nb_rows/nb_file) #round up the number
header_line = data[1,]

nb_row_inf = 1 #row number of the header line

for (id in 1:nb_file) {
  nb_row_inf = nb_row_inf + 1
  nb_row_sup = id*nb_rows_per_file
  out_df = rbind(header_line, data[nb_row_inf:nb_row_sup, ])
  write.table(out_df, 
              file = paste("refined_data/",output_basename,id,".csv", sep=""), 
              row.names = FALSE, col.names = FALSE, sep=";")
  nb_row_inf = nb_row_sup
}

