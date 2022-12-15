# Manipulation sur le fichier réduit
filepath = "raw_data/doc_function_reduced.csv"
data = read.csv(filepath, sep="\t")
data


#functions = data$Description
#dict_function = table(functions)

#barplot(dict_function)#this is a first test

#functions_table = data.frame(data[2:36], data$Description)
function_table = data.frame(data[2:35], data$Description, data$sum)





categorie = c("Control", "Colistine", "Control", "Control", "Control", "Control",
              "Colistine", "Colistine", "Control", "Colistine", "Colistine", "Control",
              "Control", "Control", "Control", "Control", "Control", "Colistine",
              "Colistine", "Colistine", "Colistine", "Colistine", "Control", "Control",
              "Control", "Control", "Control", "Colistine", "Control", "Colistine",
              "Colistine", "Colistine", "Colistine", "Colistine")


metabolic_function = c()
animal = c()
relative_abundance = c()
categorie_concerne = c()
#somme = c()



for (i in 1:(length(function_table)-2)){
  
  verif_func = c()
  verif_pig =c()
  verif_sum = c()
  verif_animal = c()
  verif_categorie = c()
  

  for (j in 1:length(function_table[,i])){
    #function_table[j,i] = function_table[j,i]/data$sum[j]
    
    #somme = c(somme, function_table$data.sum[j]) 
    #metabolic_function=c(metabolic_function, function_table$data.Description[j])
    #animal=c(animal,paste0("p",i))
    #relative_abundance=c(relative_abundance,function_table[j,i])
    #categorie_concerne = c(categorie_concerne,categorie[i])
    
    if (function_table$data.Description[j] %in% verif_func){
      posiO = match(function_table$data.Description[j],verif_func)
      verif_pig[posiO] =c(verif_pig[posiO] + function_table[j,i])
      verif_sum[posiO] =c(verif_sum[posiO] + function_table$data.sum[j])
      
    }
    else {
      verif_func = c(verif_func,function_table$data.Description[j])
      verif_pig =c(verif_pig,function_table[j,i])
      verif_sum = c(verif_sum,function_table$data.sum[j])
      verif_animal = c(verif_animal,paste0("p",i))
      verif_categorie = c(verif_categorie,categorie[i])
    }
  }
  verif_pig_relative = verif_pig / verif_sum# pas sûr que ce soit okk sur le jeu entier, à tester
  #somme = c(somme, verif_sum) 
  metabolic_function=c(metabolic_function, verif_func)
  animal=c(animal,verif_animal)
  relative_abundance=c(relative_abundance,verif_pig_relative)
  categorie_concerne = c(categorie_concerne,verif_categorie)
  
}



test = data.frame(animal,metabolic_function,relative_abundance,categorie_concerne)


################################################################################

# Library
library(ggplot2)




# Heatmap 
ggplot(test, aes(animal,metabolic_function, fill= relative_abundance)) + 
  #geom_tile()
  
  
  geom_tile(colour = "white", size = 0.5)+
  scale_fill_distiller(palette = "Spectral", 
                       limits = c(min(test$relative_abundance), max(test$relative_abundance)),
                       direction = -1) +
  scale_y_discrete(expand=c(0, 0))+
  scale_x_discrete(expand=c(0, 0))+
  labs(fill="abondance relative", title = "Abondance relative des fonctions par échantillon")+
  theme_grey(base_size=12)+
  theme(line = element_blank(),
        #axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(face="bold.italic"),
        #axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border=element_blank(),
        legend.text=element_text(face="bold", size = 15),
        plot.title = element_text(face="bold", size = 25),
        legend.title = element_text(face="bold", size = 15),
        legend.key.size = unit(2, "cm"),
        strip.text.y = element_text(angle = 0, face = "bold", size = 15),
        strip.text.x = element_text(face = "bold", size = 15))+
  facet_grid(#rows = vars(data$regne), 
             cols = vars(factor(test$categorie_concerne, levels = c("Control", "Colistine"))),
             scales = "free",
             space = "free",
             labeller = )
scale = 200

