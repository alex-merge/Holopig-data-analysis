## MAKING ASSOCIATION BETWEEN THE SAMPLE AND ITS CONDITION ##
name_sample = c()
for (nb_pig in 1:34) {
  if (nb_pig < 10){
    name_sample = c(name_sample, paste("Holopig_colon0", nb_pig, ".featureCounts.tsv", sep = ""))
  } else {
    name_sample = c(name_sample, paste("Holopig_colon", nb_pig, ".featureCounts.tsv", sep = ""))
  }
}
category = c("Control", "Colistine", "Control", "Control", "Control", "Control",
             "Colistine", "Colistine", "Control", "Colistine", "Colistine", "Control",
             "Control", "Control", "Control", "Control", "Control", "Colistine",
             "Colistine", "Colistine", "Colistine", "Colistine", "Control", "Control",
             "Control", "Control", "Control", "Colistine", "Control", "Colistine",
             "Colistine", "Colistine", "Colistine", "Colistine")
data_cat_sample = data.frame(name_sample, category)

## COLLECTING THE FILE NAMES OF SPLIT TABLE ##
directory=list.files(path = "refined_data", pattern = "Annot") #Retrieve file name in path witch contains the pattern

## WORK ON EVERY SPLIT FILE ##
# Initializing table for the for-loop
columns = c("Description", name_sample, "Sum")
functions_general_table = data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(functions_general_table) = columns

for (file in directory) {
  
  # Manipulation of each split file
  filepath = paste("refined_data/", file, sep = "")
  data = read.csv(filepath, sep = ";")
  
  # Keeping columns of interest and creating a new dataframe
  functions_table_file = data.frame(data$Description, data[, grepl(".featureCounts.tsv$", names(data))], data$sum)
  colnames(functions_table_file)[colnames(functions_table_file) == 'data.Description'] <- 'Description'
  colnames(functions_table_file)[colnames(functions_table_file) == 'data.sum'] <- 'Sum'
  
  # Summing columns by the description of the function in a dataframe
  functions_light_table = aggregate(functions_table_file[, colnames(functions_table_file)[colnames(functions_table_file) != 'Description']], list(functions_table_file$Description), FUN=sum)
  colnames(functions_light_table)[colnames(functions_light_table) == 'Group.1'] <- 'Description'
  
  # Merging dataframe from the others files and the actual file
  functions_general_table = merge(functions_general_table, functions_light_table, all = TRUE)
  functions_general_table = aggregate(functions_general_table[, colnames(functions_general_table)[colnames(functions_general_table) != 'Description']], list(functions_general_table$Description), FUN=sum)
  colnames(functions_general_table)[colnames(functions_general_table) == 'Group.1'] <- 'Description'
  
}

## TREATMENT ON ALL THE FUNCTIONS
# Calculating the relative abundance
functions_general_abdrev = functions_general_table[, colnames(functions_general_table)[colnames(functions_general_table) != 'Sum']] #Exclusion of the sum column (useless)

for (name in name_sample){
  functions_general_abdrev[,name] = functions_general_table[,name] / functions_general_table$Sum #Dividing each cell by the sum of the row
}


#metabolic_function = c()
#animal = c()
#relative_abundance = c()
#categorie_concerne = c()
#somme = c()



# for (i in 1:(length(function_table)-2)){  # = nb colonne - 2
#   
#   verif_func = c()
#   verif_pig =c()
#   verif_sum = c()
#   verif_animal = c()
#   verif_categorie = c()
#   
# 
#   for (j in 1:length(function_table[,i])){ # = nb ligne dans une colonne donnée
#     #function_table[j,i] = function_table[j,i]/data$sum[j]
#     
#     #somme = c(somme, function_table$data.sum[j]) 
#     #metabolic_function=c(metabolic_function, function_table$data.Description[j])
#     #animal=c(animal,paste0("p",i))
#     #relative_abundance=c(relative_abundance,function_table[j,i])
#     #categorie_concerne = c(categorie_concerne,categorie[i])
#     
#     if (function_table$data.Description[j] %in% verif_func){
#       posiO = match(function_table$data.Description[j],verif_func)
#       verif_pig[posiO] =c(verif_pig[posiO] + function_table[j,i])
#       verif_sum[posiO] =c(verif_sum[posiO] + function_table$data.sum[j])
#       
#     }
#     else {
#       verif_func = c(verif_func,function_table$data.Description[j])
#       verif_pig =c(verif_pig,function_table[j,i])
#       verif_sum = c(verif_sum,function_table$data.sum[j])
#       verif_animal = c(verif_animal,paste0("p",i))
#       verif_categorie = c(verif_c xdf ,categorie[i])
#     }
#   }
#   verif_pig_relative = verif_pig / verif_sum# pas sûr que ce soit okk sur le jeu entier, à tester
#   #somme = c(somme, verif_sum) 
#   metabolic_function=c(metabolic_function, verif_func)
#   animal=c(animal,verif_animal)
#   relative_abundance=c(relative_abundance,verif_pig_relative)
#   categorie_concerne = c(categorie_concerne,verif_categorie)
#   
# }
# 
# 
# 
# test = data.frame(animal,metabolic_function,relative_abundance,categorie_concerne)


################################################################################

# Library
library(ggplot2)




# # Heatmap 
# ggplot(test, aes(animal,metabolic_function, fill= relative_abundance)) + 
#   #geom_tile()
#   
#   
#   geom_tile(colour = "white", size = 0.5)+
#   scale_fill_distiller(palette = "Spectral", 
#                        limits = c(min(test$relative_abundance), max(test$relative_abundance)),
#                        direction = -1) +
#   scale_y_discrete(expand=c(0, 0))+
#   scale_x_discrete(expand=c(0, 0))+
#   labs(fill="abondance relative", title = "Abondance relative des fonctions par échantillon")+
#   theme_grey(base_size=12)+
#   theme(line = element_blank(),
#         #axis.ticks.x = element_blank(),
#         #axis.text.x = element_blank(),
#         axis.text.y = element_text(face="bold.italic"),
#         #axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         panel.border=element_blank(),
#         legend.text=element_text(face="bold", size = 15),
#         plot.title = element_text(face="bold", size = 25),
#         legend.title = element_text(face="bold", size = 15),
#         legend.key.size = unit(2, "cm"),
#         strip.text.y = element_text(angle = 0, face = "bold", size = 15),
#         strip.text.x = element_text(face = "bold", size = 15))+
#   facet_grid(#rows = vars(data$regne), 
#              cols = vars(factor(test$categorie_concerne, levels = c("Control", "Colistine"))),
#              scales = "free",
#              space = "free",
#              labeller = )
# scale = 200

