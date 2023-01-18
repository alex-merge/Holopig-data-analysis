# Library
library(ggplot2)
library(hrbrthemes)
library(tidyverse)
library(tidyr)
library(RColorBrewer)
library(reshape)


## MAKING ASSOCIATION BETWEEN THE SAMPLE AND ITS CONDITION ##
name_sample = c()
for (nb_pig in 1:34) {
  if (nb_pig < 10){
    name_sample = c(name_sample, paste("Holopig_colon0", nb_pig, ".featureCounts.tsv", sep = ""))
  } else {
    name_sample = c(name_sample, paste("Holopig_colon", nb_pig, ".featureCounts.tsv", sep = ""))
  }
}
category = c("Control", "Colistin", "Control", "Control", "Control", "Control",
             "Colistin", "Colistin", "Control", "Colistin", "Colistin", "Control",
             "Control", "Control", "Control", "Control", "Control", "Colistin",
             "Colistin", "Colistin", "Colistin", "Colistin", "Control", "Control",
             "Control", "Control", "Control", "Colistin", "Control", "Colistin",
             "Colistin", "Colistin", "Colistin", "Colistin")
data_cat_sample = data.frame(name_sample, category)


## COLLECTING THE FILE NAMES OF SPLIT TABLE ##
directory=list.files(path = "refined_data", pattern = "Annot") #Retrieve file name in path witch contains the pattern

## WORK ON PREVIOULY SPLIT FILE ##
# Initializing table for the for-loop
columns_reads = c(name_sample)
rows_reads = c(directory)
nb_reads_pig = data.frame(matrix(nrow = length(directory), ncol = length(columns_reads)))
colnames(nb_reads_pig) = columns_reads
rownames(nb_reads_pig) = rows_reads

columns = c("EC_name", name_sample, "Sum")
enzymes_general_table = data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(enzymes_general_table) = columns


for (file in directory) {
  
  # Manipulation of each split file
  filepath = paste("refined_data/", file, sep = "")
  data = read.csv(filepath, sep = ";")
  
  # Counter of total reads for each pig in each file (no filter)
  for (name_pig in name_sample) {
    sum_reads_pig = sum(na.omit(data[, colnames(data)[colnames(data) == name_pig]]))
    nb_reads_pig[file, name_pig] = sum_reads_pig
  }
  
  # If non-empty, keeping columns of interest and creating a new dataframe
  data_mod = subset(data, EC != "" & sum != 0)
  enzymes_table_file = data.frame(data_mod$EC, data_mod[, grepl(".featureCounts.tsv$", names(data_mod))], data_mod$sum)
  colnames(enzymes_table_file)[colnames(enzymes_table_file) == 'data_mod.EC'] <- 'EC_name'
  colnames(enzymes_table_file)[colnames(enzymes_table_file) == 'data_mod.sum'] <- 'Sum'
  
  # Summing columns by the name of the function in a dataframe, dropping code column
  enzymes_light_table = aggregate(x = enzymes_table_file[, colnames(enzymes_table_file)[colnames(enzymes_table_file) != 'EC_name' ]],list(enzymes_table_file$EC_name), FUN=sum)
  colnames(enzymes_light_table)[colnames(enzymes_light_table) == 'Group.1'] <- 'EC_name'
  
  # Merging dataframes with read numbers one after the other
  enzymes_general_table = merge(enzymes_general_table, enzymes_light_table, all = TRUE)
  enzymes_general_table = aggregate(enzymes_general_table[, colnames(enzymes_general_table)[colnames(enzymes_general_table) != 'EC_name']], list(enzymes_general_table$EC_name), FUN=sum)
  colnames(enzymes_general_table)[colnames(enzymes_general_table) == 'Group.1'] <- 'EC_name'
  
}


## TREATMENT ON ALL THE REMAINING FUNCTIONS
# Setting the for loop
enzymes_rel_abd = enzymes_general_table %>% select(-c(Sum)) #Exclusion of the sum column (useless in the new dataframe)


# Calculating the relative abundance
for (name in name_sample){
  total_sum_pig = sum(nb_reads_pig[[name]])
  enzymes_rel_abd[,name] = enzymes_general_table[,name] / total_sum_pig #Dividing each cell by the sum of the individual and add it to the new dataframe
}

enzymes_rel_abd_sum = separate_rows(enzymes_rel_abd, EC_name, sep = ",", convert = TRUE)
enzymes_rel_abd_sum = enzymes_rel_abd_sum[(grepl(".{+}", enzymes_rel_abd_sum$EC_name)),]
enzymes_rel_abd_sum = aggregate(enzymes_rel_abd_sum[, colnames(enzymes_rel_abd_sum)[colnames(enzymes_rel_abd_sum) != 'EC_name']], list(enzymes_rel_abd_sum$EC_name), FUN=sum)
colnames(enzymes_rel_abd_sum)[colnames(enzymes_rel_abd_sum) == 'Group.1'] <- 'EC_name'

enzymes_rel_abd_sum[(grepl("-", enzymes_rel_abd_sum$EC_name)),'EC_name'] <- 'Unclassified'

enzymes_rel_abd_sum_testColisitin = enzymes_rel_abd_sum %>% select(c(EC_name, data_cat_sample$name_sample[data_cat_sample$category == "Colistin"])) %>% mutate(nb_empty_sample_Colisitin = rowSums(. == 0)) %>% select(c(EC_name,nb_empty_sample_Colisitin))
enzymes_rel_abd_sum_testControl = enzymes_rel_abd_sum %>% select(c(EC_name, data_cat_sample$name_sample[data_cat_sample$category == "Control"])) %>% mutate(nb_empty_sample_Control = rowSums(. == 0)) %>% select(c(EC_name,nb_empty_sample_Control))
enzymes_rel_abd_sum = subset(enzymes_rel_abd_sum, (enzymes_rel_abd_sum_testColisitin$nb_empty_sample_Colisitin <= 8 & enzymes_rel_abd_sum_testControl$nb_empty_sample_Control <= 8) == TRUE)

# Add the category of the pig (colistin or control)
transposed_cat_pig = t(data_cat_sample)
dataframe_pig_cat = data.frame(transposed_cat_pig)
names(dataframe_pig_cat) = as.matrix(dataframe_pig_cat[1,]) #make the first row as the header
dataframe_pig_cat = dataframe_pig_cat[-1,] #remove the first row now



##########################################HEATMAP##################################################

name_enzyme_code = read.csv("ressources/EC_equivalences.csv", sep = ";")

# FIRST HEATMAP OF PURE FUNCTIONS

enzymes_rel_abd_sum_ord = order(rowMeans(enzymes_rel_abd_sum[,-1]), decreasing = TRUE)
enzymes_rel_abd_final = enzymes_rel_abd_sum[enzymes_rel_abd_sum_ord[1:34],]

#Reshaping the data for the heatmap
reshaped_enzymes = data.frame(matrix(nrow = 0, ncol = 4))
colnames(reshaped_enzymes) = c("EC_description","Relative_abundance","Pig_name","Category")

for (name_pig in name_sample) {
  reshaped_enzymes_pig = data.frame(enzymes_rel_abd_final$EC_name, enzymes_rel_abd_final[colnames(enzymes_rel_abd_final) == name_pig])
  reshaped_enzymes_pig['pig'] = rep(name_pig, length(enzymes_rel_abd_final$EC_name))
  reshaped_enzymes_pig['category'] = rep(dataframe_pig_cat[1,colnames(dataframe_pig_cat) == name_pig], length(enzymes_rel_abd_final$EC_name))
  colnames(reshaped_enzymes_pig)[colnames(reshaped_enzymes_pig) == name_pig] = 'relative_abundance'
  reshaped_enzymes = merge(reshaped_enzymes, reshaped_enzymes_pig, all = TRUE, by.y = c("enzymes_rel_abd_final.EC_name","relative_abundance","pig","category"), by.x = c("EC_description","Relative_abundance","Pig_name", "Category"))
}
reshaped_enzymes = reshaped_enzymes %>% mutate(Pig_name = gsub("(\\D)*", "", Pig_name))
reshaped_enzymes$EC_description = name_enzyme_code$EC_description[match(reshaped_enzymes$EC_description, name_enzyme_code$EC_name)]
reshaped_enzymes$EC_description = replace_na(reshaped_enzymes$EC_description,"Unclassified")

ggplot(reshaped_enzymes, aes(Pig_name, reorder(EC_description, Relative_abundance), fill= Relative_abundance)) +
  geom_tile(colour = "white", size = 0.5) +
  scale_fill_distiller(palette = "Spectral",
                       limits = c(min(reshaped_enzymes$Relative_abundance), max(reshaped_enzymes$Relative_abundance)),
                       direction = -1) +
  scale_y_discrete(expand=c(0, 0))+
  scale_x_discrete(expand=c(0, 0))+
  labs(fill="Sum of relative abundance", title = "Sum of relative abundance of enzymes by pig")+
  theme_grey(base_size=12)+
  theme(line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x.bottom = element_text(face="bold.italic"),
        axis.text.y = element_text(face="bold.italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border=element_blank(),
        legend.text=element_text(face="bold", size = 15),
        plot.title = element_text(face="bold", size = 20),
        legend.title = element_text(face="bold", size = 15),
        legend.key.size = unit(2, "cm"),
        strip.text.y = element_text(angle = 0, face = "bold", size = 15),
        strip.text.x = element_text(face = "bold", size = 15))+
  facet_grid(cols = vars(factor(reshaped_enzymes$Category, levels = c("Control", "Colistin"))),
             scales = "free",
             space = "free",
             labeller = )
scale = 200
ggsave(filename = "export/HM_enzymes_functions.png",
       width = 16*scale,
       height = 18*scale,
       units = "px",
       dpi=200)

reshaped_pure_functions_no_unknown = reshaped_pure_functions[-grep(pattern = 'Function unknown|Unclassified', reshaped_pure_functions$COG_cat),]

ggplot(reshaped_pure_functions_no_unknown[], aes(Pig_name, reorder(COG_cat, Relative_abundance), fill= Relative_abundance)) +
  geom_tile(colour = "white", size = 0.5) +
  scale_fill_distiller(palette = "Spectral",
                       limits = c(min(reshaped_pure_functions_no_unknown$Relative_abundance), max(reshaped_pure_functions_no_unknown$Relative_abundance)),
                       direction = -1) +
  scale_y_discrete(expand=c(0, 0))+
  scale_x_discrete(expand=c(0, 0))+
  labs(fill="Sum of relative abundance", title = "Sum of relative abundance of functions by pig without unknown functions")+
  theme_grey(base_size=12)+
  theme(line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x.bottom = element_text(face="bold.italic"),
        axis.text.y = element_text(face="bold.italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border=element_blank(),
        legend.text=element_text(face="bold", size = 15),
        plot.title = element_text(face="bold", size = 20),
        legend.title = element_text(face="bold", size = 15),
        legend.key.size = unit(2, "cm"),
        strip.text.y = element_text(angle = 0, face = "bold", size = 15),
        strip.text.x = element_text(face = "bold", size = 15))+
  facet_grid(cols = vars(factor(reshaped_pure_functions_no_unknown$Category, levels = c("Control", "Colistin"))),
             scales = "free",
             space = "free",
             labeller = )
scale = 200
ggsave(filename = "export/HM_pure_functions_no_unkown.png",
       width = 16*scale,
       height = 18*scale,
       units = "px",
       dpi=200)
