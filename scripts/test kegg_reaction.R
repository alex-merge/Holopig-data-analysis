# Library
library(ggplot2)
library(hrbrthemes)
library(tidyverse)
library(tidyr)
library(RColorBrewer)

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

columns = c("COG_cat", name_sample, "Sum")
functions_general_table = data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(functions_general_table) = columns

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
  data_mod = subset(data, KEGG_Reaction != "" & sum != 0)
  functions_table_file = data.frame(data_mod$KEGG_Reaction, data_mod[, grepl(".featureCounts.tsv$", names(data_mod))], data_mod$sum)
  colnames(functions_table_file)[colnames(functions_table_file) == 'data_mod.KEGG_Reaction'] <- 'COG_cat'
  colnames(functions_table_file)[colnames(functions_table_file) == 'data_mod.sum'] <- 'Sum'
  
  # Summing columns by the name of the function in a dataframe, dropping code column
  functions_light_table = aggregate(x = functions_table_file[, colnames(functions_table_file)[colnames(functions_table_file) != 'COG_cat' ]],list(functions_table_file$COG_cat), FUN=sum)
  colnames(functions_light_table)[colnames(functions_light_table) == 'Group.1'] <- 'COG_cat'
  
  # Merging dataframes with read numbers one after the other
  functions_general_table = merge(functions_general_table, functions_light_table, all = TRUE)
  functions_general_table = aggregate(functions_general_table[, colnames(functions_general_table)[colnames(functions_general_table) != 'COG_cat']], list(functions_general_table$COG_cat), FUN=sum)
  colnames(functions_general_table)[colnames(functions_general_table) == 'Group.1'] <- 'COG_cat'
  
}


## TREATMENT ON ALL THE REMAINING FUNCTIONS
# Setting the for loop
functions_rel_abd = functions_general_table #%>% select(-c(Sum)) #Exclusion of the sum column (useless in the new dataframe)


# Calculating the relative abundance
for (name in name_sample){
  total_sum_pig = sum(nb_reads_pig[[name]])
  functions_rel_abd[,name] = functions_rel_abd[,name] / total_sum_pig #Dividing each cell by the sum of the individual and add it to the new dataframe
}

functions_rel_abd[1,'COG_cat'] <- 'Unclassified'

# Add the category of the pig (colistin or control)
transposed_cat_pig = t(data_cat_sample)
dataframe_pig_cat = data.frame(transposed_cat_pig)
names(dataframe_pig_cat) = as.matrix(dataframe_pig_cat[1,]) #make the first row as the header
dataframe_pig_cat = dataframe_pig_cat[-1,] #remove the first row now

#functions_rel_abd = tidyr::separate_rows(functions_rel_abd,COG_cat,sep =',')
#functions_rel_abd = functions_rel_abd[order(functions_rel_abd$Sum),]
#functions_rel_abd = top_n(functions_rel_abd,60)
#functions_rel_abd = functions_rel_abd %>% select(-c(Sum)) #Exclusion of the sum column (useless in the new dataframe)


functions_rel_abd = tidyr::separate_rows(functions_rel_abd,COG_cat,sep =',')
#functions_rel_abd = functions_rel_abd[order(functions_rel_abd$Sum),]

######
functions_rel_abd_sum = separate_rows(functions_rel_abd, COG_cat, sep = ",", convert = TRUE)
functions_rel_abd_sum = functions_rel_abd_sum[(grepl(".{+}", functions_rel_abd_sum$COG_cat)),]
functions_rel_abd_sum = aggregate(functions_rel_abd_sum[, colnames(functions_rel_abd_sum)[colnames(functions_rel_abd_sum) != 'COG_cat']], list(functions_rel_abd_sum$COG_cat), FUN=sum)
colnames(functions_rel_abd_sum)[colnames(functions_rel_abd_sum) == 'Group.1'] <- 'COG_cat'
functions_rel_abd = functions_rel_abd_sum
######


functions_rel_abd_ordered = order(rowMeans(functions_rel_abd[,-1]), decreasing = TRUE)


functions_rel_abd = functions_rel_abd[functions_rel_abd_ordered[1:100],]

#functions_rel_abd = top_n(functions_rel_abd,100)
##########################################HEATMAP##################################################

# SECOND HEATMAP OF PURE FUNCTIONS + MELTED FUNCTIONS

#Reshaping the data for the heatmap
reshaped_melt_functions = data.frame(matrix(nrow = 0, ncol = 4))
colnames(reshaped_melt_functions) = c("COG_cat","Relative_abundance","Pig_name","Category")

#tidyr::separate_rows(functions_rel_abd,COG_cat,sep =',')          ##########

for (name_pig in name_sample) {
  
  cog_melt_functions = functions_rel_abd
  reshaped_melt_functions_pig = data.frame(cog_melt_functions$COG_cat, cog_melt_functions[colnames(cog_melt_functions) == name_pig])
  reshaped_melt_functions_pig['pig'] = rep(name_pig, length(cog_melt_functions$COG_cat))
  reshaped_melt_functions_pig['category'] = rep(dataframe_pig_cat[1,colnames(dataframe_pig_cat) == name_pig], length(cog_melt_functions$COG_cat))
  colnames(reshaped_melt_functions_pig)[colnames(reshaped_melt_functions_pig) == name_pig] = 'relative_abundance'
  reshaped_melt_functions = merge(reshaped_melt_functions, reshaped_melt_functions_pig, all = TRUE, by.y = c("cog_melt_functions.COG_cat","relative_abundance","pig","category"), by.x = c("COG_cat","Relative_abundance","Pig_name", "Category"))
}
reshaped_all_functions = reshaped_melt_functions %>% mutate(Pig_name = gsub("(\\D)*", "", Pig_name))
#reshaped_all_functions = merge(reshaped_pure_functions, reshaped_melt_functions, all=TRUE, by=c("COG_cat","Relative_abundance","Pig_name", "Category"))
#reshaped_all_functions = tidyr::separate_rows(reshaped_melt_functions,COG_cat,sep =',')

reshaped_all_functions = subset(reshaped_all_functions, COG_cat!= 'Unclassified')



ggplot(reshaped_all_functions, aes(Pig_name, reorder(COG_cat, Relative_abundance), fill= Relative_abundance)) +
  geom_tile(colour = "white", size = 0.5) +
  scale_fill_distiller(palette = "Spectral",
                       limits = c(min(reshaped_all_functions$Relative_abundance), max(reshaped_all_functions$Relative_abundance)),
                       direction = -1) +
  scale_y_discrete(expand=c(0, 0))+
  scale_x_discrete(expand=c(0, 0))+
  labs(fill="Relative Abundance", title = "Relative Abundance of all the functions by pig")+
  theme_grey(base_size=12)+
  theme(line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x.bottom = element_text(face="bold.italic"),
        axis.text.y = element_text(face="bold.italic", size = 5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border=element_blank(),
        legend.text=element_text(face="bold", size = 15),
        plot.title = element_text(face="bold", size = 25),
        legend.title = element_text(face="bold", size = 15),
        legend.key.size = unit(2, "cm"),
        strip.text.y = element_text(angle = 0, face = "bold", size = 15),
        strip.text.x = element_text(face = "bold", size = 15))+
  facet_grid(cols = vars(factor(reshaped_all_functions$Category, levels = c("Control", "Colistin"))),
             scales = "free",
             space = "free",
             labeller = )
scale = 200

ggsave(filename = "export/Metabolics_functions_kegg_reaction.png",
       width = 16*scale,
       height = 18*scale,
       units = "px",
       dpi=200)