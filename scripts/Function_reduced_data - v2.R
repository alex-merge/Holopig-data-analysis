# Manipulation sur le fichier réduit
filepath = "raw_data/doc_function_reduced.csv"
data = read.csv(filepath, sep="\t")
data

# Reduction of the input file
filepath = "raw_data/doc_function_reduced.csv"


write.table(log_data, file = paste(substr(tmp_dir_path, 1, nchar(tmp_dir_path)), 'P311-infos_cnvfile_input.txt', sep = '/'), quote = FALSE, row.names = TRUE, col.names = FALSE)


#functions = data$Description
#dict_function = table(functions)

#barplot(dict_function)#this is a first test

#functions_table = data.frame(data[2:36], data$Description)
function_table = data.frame(data[2:35], data$Description)
metabolic_function = c()
animal = c()
relative_abundance = c()

for (i in 1:(length(function_table)-1)){

  for (j in 1:length(function_table[,i])){
    function_table[j,i] = function_table[j,i]/data$sum[j]
    metabolic_function=c(metabolic_function,data$Description[j])
    animal=c(animal,paste0("pig_",i))
    relative_abundance=c(relative_abundance,function_table[j,i])
    
    
  }

}
test = data.frame(metabolic_function,animal,relative_abundance)

# Library
library(ggplot2)

# Heatmap 
ggplot(test, aes(metabolic_function, animal, fill= relative_abundance)) + 
  geom_tile()

