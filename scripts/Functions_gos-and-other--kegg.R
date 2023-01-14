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
  data_mod = subset(data, KEGG_ko != "" & sum != 0)                    ########
  functions_table_file = data.frame(data_mod$KEGG_ko, data_mod[, grepl(".featureCounts.tsv$", names(data_mod))], data_mod$sum)    #########
  colnames(functions_table_file)[colnames(functions_table_file) == 'data_mod.KEGG_ko'] <- 'COG_cat'
  colnames(functions_table_file)[colnames(functions_table_file) == 'data_mod.sum'] <- 'Sum'
  
  # # If empty or with no functional correspondence, putting in another file to keep trace of the non-analysed data
  # data_useless <- subset(data, COG_category == "-" | COG_category == "" | sum == 0)
  # write.table(data_useless,
  #             append = TRUE,
  #             file = "refined_data/non_use_metabolics_functions.csv",
  #             row.names = FALSE, col.names = TRUE, sep = ";", dec = ".")
  
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
functions_rel_abd = functions_general_table %>% select(-c(Sum)) #Exclusion of the sum column (useless in the new dataframe)


# Calculating the relative abundance
for (name in name_sample){
  total_sum_pig = sum(nb_reads_pig[[name]])
  functions_rel_abd[,name] = functions_rel_abd[,name] / total_sum_pig #Dividing each cell by the sum of the individual and add it to the new dataframe
}

# rownames(functions_rel_abd) = as.matrix(functions_rel_abd[,1]) #make the first column as the header
# functions_rel_abd = functions_rel_abd[,-1] #remove the first column now
functions_rel_abd[1,'COG_cat'] <- 'Unknown'


# Add the category of the pig (colistin or control)
transposed_cat_pig = t(data_cat_sample)
dataframe_pig_cat = data.frame(transposed_cat_pig)
names(dataframe_pig_cat) = as.matrix(dataframe_pig_cat[1,]) #make the first row as the header
dataframe_pig_cat = dataframe_pig_cat[-1,] #remove the first row now



##########################################HEATMAP##################################################

# FIRST HEATMAP OF PURE FUNCTIONS

#Reshaping the data for the heatmap
reshaped_pure_functions = data.frame(matrix(nrow = 0, ncol = 4))
colnames(reshaped_pure_functions) = c("COG_cat","Relative_abundance","Pig_name","Category")

for (name_pig in name_sample) {
  #cog_pure_functions = functions_rel_abd[nchar(as.character(functions_rel_abd$COG_cat)) == '1' | functions_rel_abd$COG_cat == 'Unknown', ]
  cog_pure_functions = functions_rel_abd[nchar(as.character(functions_rel_abd$COG_cat)) == '9' | functions_rel_abd$COG_cat == 'Unknown', ]
  reshaped_pure_functions_pig = data.frame(cog_pure_functions$COG_cat, cog_pure_functions[colnames(cog_pure_functions) == name_pig])
  reshaped_pure_functions_pig['pig'] = rep(name_pig, length(cog_pure_functions$COG_cat))
  reshaped_pure_functions_pig['category'] = rep(dataframe_pig_cat[1,colnames(dataframe_pig_cat) == name_pig], length(cog_pure_functions$COG_cat))
  colnames(reshaped_pure_functions_pig)[colnames(reshaped_pure_functions_pig) == name_pig] = 'relative_abundance'
  reshaped_pure_functions = merge(reshaped_pure_functions, reshaped_pure_functions_pig, all = TRUE, by.y = c("cog_pure_functions.COG_cat","relative_abundance","pig","category"), by.x = c("COG_cat","Relative_abundance","Pig_name", "Category"))
}

reshaped_pure_functions = reshaped_pure_functions %>% mutate(Pig_name = gsub("(\\D)*", "", Pig_name))

ggplot(reshaped_pure_functions, aes(Pig_name, COG_cat, fill= Relative_abundance)) +
  geom_tile(colour = "white", size = 0.5) +
  scale_fill_distiller(palette = "Spectral",
                       limits = c(min(reshaped_pure_functions$Relative_abundance), max(reshaped_pure_functions$Relative_abundance )),
                       direction = -1) +
  scale_y_discrete(expand=c(0, 0))+
  scale_x_discrete(expand=c(0, 0))+
  labs(fill="Relative Abundance", title = "Relative Abundance of functions by pig")+
  theme_grey(base_size=12)+
  theme(line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x.bottom = element_text(),
        axis.text.y = element_text(face="bold.italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border=element_blank(),
        legend.text=element_text(face="bold", size = 15),
        plot.title = element_text(face="bold", size = 25),
        legend.title = element_text(face="bold", size = 15),
        legend.key.size = unit(2, "cm"),
        strip.text.y = element_text(angle = 0, face = "bold", size = 15),
        strip.text.x = element_text(face = "bold", size = 15))+
  facet_grid(cols = vars(factor(reshaped_pure_functions$Category, levels = c("Control", "Colistin"))),
             scales = "free",
             space = "free",
             labeller = )
scale = 200








################################################
################################################
##############################################"





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
  data_mod = subset(data, GOs != "" & sum != 0)                    ########
  functions_table_file = data.frame(data_mod$GOs, data_mod[, grepl(".featureCounts.tsv$", names(data_mod))], data_mod$sum)    #########
  colnames(functions_table_file)[colnames(functions_table_file) == 'data_mod.GOs'] <- 'COG_cat'
  colnames(functions_table_file)[colnames(functions_table_file) == 'data_mod.sum'] <- 'Sum'
  
  # # If empty or with no functional correspondence, putting in another file to keep trace of the non-analysed data
  # data_useless <- subset(data, COG_category == "-" | COG_category == "" | sum == 0)
  # write.table(data_useless,
  #             append = TRUE,
  #             file = "refined_data/non_use_metabolics_functions.csv",
  #             row.names = FALSE, col.names = TRUE, sep = ";", dec = ".")
  
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
functions_rel_abd = functions_general_table %>% select(-c(Sum)) #Exclusion of the sum column (useless in the new dataframe)


# Calculating the relative abundance
for (name in name_sample){
  total_sum_pig = sum(nb_reads_pig[[name]])
  functions_rel_abd[,name] = functions_rel_abd[,name] / total_sum_pig #Dividing each cell by the sum of the individual and add it to the new dataframe
}

# rownames(functions_rel_abd) = as.matrix(functions_rel_abd[,1]) #make the first column as the header
# functions_rel_abd = functions_rel_abd[,-1] #remove the first column now
functions_rel_abd[1,'COG_cat'] <- 'Unknown'


# Add the category of the pig (colistin or control)
transposed_cat_pig = t(data_cat_sample)
dataframe_pig_cat = data.frame(transposed_cat_pig)
names(dataframe_pig_cat) = as.matrix(dataframe_pig_cat[1,]) #make the first row as the header
dataframe_pig_cat = dataframe_pig_cat[-1,] #remove the first row now



##########################################HEATMAP##################################################

# FIRST HEATMAP OF PURE FUNCTIONS

#Reshaping the data for the heatmap
reshaped_pure_functions = data.frame(matrix(nrow = 0, ncol = 4))
colnames(reshaped_pure_functions) = c("COG_cat","Relative_abundance","Pig_name","Category")

for (name_pig in name_sample) {
  #cog_pure_functions = functions_rel_abd[nchar(as.character(functions_rel_abd$COG_cat)) == '1' | functions_rel_abd$COG_cat == 'Unknown', ]
  cog_pure_functions = functions_rel_abd[nchar(as.character(functions_rel_abd$COG_cat)) == '11' | functions_rel_abd$COG_cat == 'Unknown', ]
  reshaped_pure_functions_pig = data.frame(cog_pure_functions$COG_cat, cog_pure_functions[colnames(cog_pure_functions) == name_pig])
  reshaped_pure_functions_pig['pig'] = rep(name_pig, length(cog_pure_functions$COG_cat))
  reshaped_pure_functions_pig['category'] = rep(dataframe_pig_cat[1,colnames(dataframe_pig_cat) == name_pig], length(cog_pure_functions$COG_cat))
  colnames(reshaped_pure_functions_pig)[colnames(reshaped_pure_functions_pig) == name_pig] = 'relative_abundance'
  reshaped_pure_functions = merge(reshaped_pure_functions, reshaped_pure_functions_pig, all = TRUE, by.y = c("cog_pure_functions.COG_cat","relative_abundance","pig","category"), by.x = c("COG_cat","Relative_abundance","Pig_name", "Category"))
}

reshaped_pure_functions = reshaped_pure_functions %>% mutate(Pig_name = gsub("(\\D)*", "", Pig_name))

ggplot(reshaped_pure_functions, aes(Pig_name, COG_cat, fill= Relative_abundance)) +
  geom_tile(colour = "white", size = 0.5) +
  scale_fill_distiller(palette = "Spectral",
                       limits = c(min(reshaped_pure_functions$Relative_abundance), max(reshaped_pure_functions$Relative_abundance )),
                       direction = -1) +
  scale_y_discrete(expand=c(0, 0))+
  scale_x_discrete(expand=c(0, 0))+
  labs(fill="Relative Abundance", title = "Relative Abundance of functions by pig")+
  theme_grey(base_size=12)+
  theme(line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x.bottom = element_text(),
        axis.text.y = element_text(face="bold.italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border=element_blank(),
        legend.text=element_text(face="bold", size = 15),
        plot.title = element_text(face="bold", size = 25),
        legend.title = element_text(face="bold", size = 15),
        legend.key.size = unit(2, "cm"),
        strip.text.y = element_text(angle = 0, face = "bold", size = 15),
        strip.text.x = element_text(face = "bold", size = 15))+
  facet_grid(cols = vars(factor(reshaped_pure_functions$Category, levels = c("Control", "Colistin"))),
             scales = "free",
             space = "free",
             labeller = )
scale = 200








############################################################################
#############################################################################
################################################################################

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
  data_mod = subset(data, KEGG_Pathway != "" & sum != 0)                    ########
  functions_table_file = data.frame(data_mod$KEGG_Pathway, data_mod[, grepl(".featureCounts.tsv$", names(data_mod))], data_mod$sum)    #########
  colnames(functions_table_file)[colnames(functions_table_file) == 'data_mod.KEGG_Pathway'] <- 'COG_cat'
  colnames(functions_table_file)[colnames(functions_table_file) == 'data_mod.sum'] <- 'Sum'
  
  # # If empty or with no functional correspondence, putting in another file to keep trace of the non-analysed data
  # data_useless <- subset(data, COG_category == "-" | COG_category == "" | sum == 0)
  # write.table(data_useless,
  #             append = TRUE,
  #             file = "refined_data/non_use_metabolics_functions.csv",
  #             row.names = FALSE, col.names = TRUE, sep = ";", dec = ".")
  
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
functions_rel_abd = functions_general_table %>% select(-c(Sum)) #Exclusion of the sum column (useless in the new dataframe)


# Calculating the relative abundance
for (name in name_sample){
  total_sum_pig = sum(nb_reads_pig[[name]])
  functions_rel_abd[,name] = functions_rel_abd[,name] / total_sum_pig #Dividing each cell by the sum of the individual and add it to the new dataframe
}

# rownames(functions_rel_abd) = as.matrix(functions_rel_abd[,1]) #make the first column as the header
# functions_rel_abd = functions_rel_abd[,-1] #remove the first column now
functions_rel_abd[1,'COG_cat'] <- 'Unknown'


# Add the category of the pig (colistin or control)
transposed_cat_pig = t(data_cat_sample)
dataframe_pig_cat = data.frame(transposed_cat_pig)
names(dataframe_pig_cat) = as.matrix(dataframe_pig_cat[1,]) #make the first row as the header
dataframe_pig_cat = dataframe_pig_cat[-1,] #remove the first row now



##########################################HEATMAP##################################################

# FIRST HEATMAP OF PURE FUNCTIONS

#Reshaping the data for the heatmap
reshaped_pure_functions = data.frame(matrix(nrow = 0, ncol = 4))
colnames(reshaped_pure_functions) = c("COG_cat","Relative_abundance","Pig_name","Category")

for (name_pig in name_sample) {
  #cog_pure_functions = functions_rel_abd[nchar(as.character(functions_rel_abd$COG_cat)) == '1' | functions_rel_abd$COG_cat == 'Unknown', ]
  cog_pure_functions = functions_rel_abd[nchar(as.character(functions_rel_abd$COG_cat)) == '8' | functions_rel_abd$COG_cat == 'Unknown', ]
  reshaped_pure_functions_pig = data.frame(cog_pure_functions$COG_cat, cog_pure_functions[colnames(cog_pure_functions) == name_pig])
  reshaped_pure_functions_pig['pig'] = rep(name_pig, length(cog_pure_functions$COG_cat))
  reshaped_pure_functions_pig['category'] = rep(dataframe_pig_cat[1,colnames(dataframe_pig_cat) == name_pig], length(cog_pure_functions$COG_cat))
  colnames(reshaped_pure_functions_pig)[colnames(reshaped_pure_functions_pig) == name_pig] = 'relative_abundance'
  reshaped_pure_functions = merge(reshaped_pure_functions, reshaped_pure_functions_pig, all = TRUE, by.y = c("cog_pure_functions.COG_cat","relative_abundance","pig","category"), by.x = c("COG_cat","Relative_abundance","Pig_name", "Category"))
}

reshaped_pure_functions = reshaped_pure_functions %>% mutate(Pig_name = gsub("(\\D)*", "", Pig_name))

ggplot(reshaped_pure_functions, aes(Pig_name, COG_cat, fill= Relative_abundance)) +
  geom_tile(colour = "white", size = 0.5) +
  scale_fill_distiller(palette = "Spectral",
                       limits = c(min(reshaped_pure_functions$Relative_abundance), max(reshaped_pure_functions$Relative_abundance )),
                       direction = -1) +
  scale_y_discrete(expand=c(0, 0))+
  scale_x_discrete(expand=c(0, 0))+
  labs(fill="Relative Abundance", title = "Relative Abundance of functions by pig")+
  theme_grey(base_size=12)+
  theme(line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x.bottom = element_text(),
        axis.text.y = element_text(face="bold.italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border=element_blank(),
        legend.text=element_text(face="bold", size = 15),
        plot.title = element_text(face="bold", size = 25),
        legend.title = element_text(face="bold", size = 15),
        legend.key.size = unit(2, "cm"),
        strip.text.y = element_text(angle = 0, face = "bold", size = 15),
        strip.text.x = element_text(face = "bold", size = 15))+
  facet_grid(cols = vars(factor(reshaped_pure_functions$Category, levels = c("Control", "Colistin"))),
             scales = "free",
             space = "free",
             labeller = )
scale = 200











##############################################################
################################################################
###################################################################


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
  data_mod = subset(data, KEGG_Module != "" & sum != 0)                    ########
  functions_table_file = data.frame(data_mod$KEGG_Module, data_mod[, grepl(".featureCounts.tsv$", names(data_mod))], data_mod$sum)    #########
  colnames(functions_table_file)[colnames(functions_table_file) == 'data_mod.KEGG_Module'] <- 'COG_cat'
  colnames(functions_table_file)[colnames(functions_table_file) == 'data_mod.sum'] <- 'Sum'
  
  # # If empty or with no functional correspondence, putting in another file to keep trace of the non-analysed data
  # data_useless <- subset(data, COG_category == "-" | COG_category == "" | sum == 0)
  # write.table(data_useless,
  #             append = TRUE,
  #             file = "refined_data/non_use_metabolics_functions.csv",
  #             row.names = FALSE, col.names = TRUE, sep = ";", dec = ".")
  
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
functions_rel_abd = functions_general_table %>% select(-c(Sum)) #Exclusion of the sum column (useless in the new dataframe)


# Calculating the relative abundance
for (name in name_sample){
  total_sum_pig = sum(nb_reads_pig[[name]])
  functions_rel_abd[,name] = functions_rel_abd[,name] / total_sum_pig #Dividing each cell by the sum of the individual and add it to the new dataframe
}

# rownames(functions_rel_abd) = as.matrix(functions_rel_abd[,1]) #make the first column as the header
# functions_rel_abd = functions_rel_abd[,-1] #remove the first column now
functions_rel_abd[1,'COG_cat'] <- 'Unknown'


# Add the category of the pig (colistin or control)
transposed_cat_pig = t(data_cat_sample)
dataframe_pig_cat = data.frame(transposed_cat_pig)
names(dataframe_pig_cat) = as.matrix(dataframe_pig_cat[1,]) #make the first row as the header
dataframe_pig_cat = dataframe_pig_cat[-1,] #remove the first row now



##########################################HEATMAP##################################################

# FIRST HEATMAP OF PURE FUNCTIONS

#Reshaping the data for the heatmap
reshaped_pure_functions = data.frame(matrix(nrow = 0, ncol = 4))
colnames(reshaped_pure_functions) = c("COG_cat","Relative_abundance","Pig_name","Category")

for (name_pig in name_sample) {
  #cog_pure_functions = functions_rel_abd[nchar(as.character(functions_rel_abd$COG_cat)) == '1' | functions_rel_abd$COG_cat == 'Unknown', ]
  cog_pure_functions = functions_rel_abd[nchar(as.character(functions_rel_abd$COG_cat)) == '7' | functions_rel_abd$COG_cat == 'Unknown', ]
  reshaped_pure_functions_pig = data.frame(cog_pure_functions$COG_cat, cog_pure_functions[colnames(cog_pure_functions) == name_pig])
  reshaped_pure_functions_pig['pig'] = rep(name_pig, length(cog_pure_functions$COG_cat))
  reshaped_pure_functions_pig['category'] = rep(dataframe_pig_cat[1,colnames(dataframe_pig_cat) == name_pig], length(cog_pure_functions$COG_cat))
  colnames(reshaped_pure_functions_pig)[colnames(reshaped_pure_functions_pig) == name_pig] = 'relative_abundance'
  reshaped_pure_functions = merge(reshaped_pure_functions, reshaped_pure_functions_pig, all = TRUE, by.y = c("cog_pure_functions.COG_cat","relative_abundance","pig","category"), by.x = c("COG_cat","Relative_abundance","Pig_name", "Category"))
}

reshaped_pure_functions = reshaped_pure_functions %>% mutate(Pig_name = gsub("(\\D)*", "", Pig_name))

ggplot(reshaped_pure_functions, aes(Pig_name, COG_cat, fill= Relative_abundance)) +
  geom_tile(colour = "white", size = 0.5) +
  scale_fill_distiller(palette = "Spectral",
                       limits = c(min(reshaped_pure_functions$Relative_abundance), max(reshaped_pure_functions$Relative_abundance )),
                       direction = -1) +
  scale_y_discrete(expand=c(0, 0))+
  scale_x_discrete(expand=c(0, 0))+
  labs(fill="Relative Abundance", title = "Relative Abundance of functions by pig")+
  theme_grey(base_size=12)+
  theme(line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x.bottom = element_text(),
        axis.text.y = element_text(face="bold.italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border=element_blank(),
        legend.text=element_text(face="bold", size = 15),
        plot.title = element_text(face="bold", size = 25),
        legend.title = element_text(face="bold", size = 15),
        legend.key.size = unit(2, "cm"),
        strip.text.y = element_text(angle = 0, face = "bold", size = 15),
        strip.text.x = element_text(face = "bold", size = 15))+
  facet_grid(cols = vars(factor(reshaped_pure_functions$Category, levels = c("Control", "Colistin"))),
             scales = "free",
             space = "free",
             labeller = )
scale = 200










##################################################################
###################################################################
####################################################################



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
  data_mod = subset(data, KEGG_Reaction != "" & sum != 0)                    ########
  functions_table_file = data.frame(data_mod$KEGG_Reaction, data_mod[, grepl(".featureCounts.tsv$", names(data_mod))], data_mod$sum)    #########
  colnames(functions_table_file)[colnames(functions_table_file) == 'data_mod.KEGG_Reaction'] <- 'COG_cat'
  colnames(functions_table_file)[colnames(functions_table_file) == 'data_mod.sum'] <- 'Sum'
  
  # # If empty or with no functional correspondence, putting in another file to keep trace of the non-analysed data
  # data_useless <- subset(data, COG_category == "-" | COG_category == "" | sum == 0)
  # write.table(data_useless,
  #             append = TRUE,
  #             file = "refined_data/non_use_metabolics_functions.csv",
  #             row.names = FALSE, col.names = TRUE, sep = ";", dec = ".")
  
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
functions_rel_abd = functions_general_table %>% select(-c(Sum)) #Exclusion of the sum column (useless in the new dataframe)


# Calculating the relative abundance
for (name in name_sample){
  total_sum_pig = sum(nb_reads_pig[[name]])
  functions_rel_abd[,name] = functions_rel_abd[,name] / total_sum_pig #Dividing each cell by the sum of the individual and add it to the new dataframe
}

# rownames(functions_rel_abd) = as.matrix(functions_rel_abd[,1]) #make the first column as the header
# functions_rel_abd = functions_rel_abd[,-1] #remove the first column now
functions_rel_abd[1,'COG_cat'] <- 'Unknown'


# Add the category of the pig (colistin or control)
transposed_cat_pig = t(data_cat_sample)
dataframe_pig_cat = data.frame(transposed_cat_pig)
names(dataframe_pig_cat) = as.matrix(dataframe_pig_cat[1,]) #make the first row as the header
dataframe_pig_cat = dataframe_pig_cat[-1,] #remove the first row now



##########################################HEATMAP##################################################

# FIRST HEATMAP OF PURE FUNCTIONS

#Reshaping the data for the heatmap
reshaped_pure_functions = data.frame(matrix(nrow = 0, ncol = 4))
colnames(reshaped_pure_functions) = c("COG_cat","Relative_abundance","Pig_name","Category")

for (name_pig in name_sample) {
  #cog_pure_functions = functions_rel_abd[nchar(as.character(functions_rel_abd$COG_cat)) == '1' | functions_rel_abd$COG_cat == 'Unknown', ]
  cog_pure_functions = functions_rel_abd[nchar(as.character(functions_rel_abd$COG_cat)) == '7' | functions_rel_abd$COG_cat == 'Unknown', ]
  reshaped_pure_functions_pig = data.frame(cog_pure_functions$COG_cat, cog_pure_functions[colnames(cog_pure_functions) == name_pig])
  reshaped_pure_functions_pig['pig'] = rep(name_pig, length(cog_pure_functions$COG_cat))
  reshaped_pure_functions_pig['category'] = rep(dataframe_pig_cat[1,colnames(dataframe_pig_cat) == name_pig], length(cog_pure_functions$COG_cat))
  colnames(reshaped_pure_functions_pig)[colnames(reshaped_pure_functions_pig) == name_pig] = 'relative_abundance'
  reshaped_pure_functions = merge(reshaped_pure_functions, reshaped_pure_functions_pig, all = TRUE, by.y = c("cog_pure_functions.COG_cat","relative_abundance","pig","category"), by.x = c("COG_cat","Relative_abundance","Pig_name", "Category"))
}

reshaped_pure_functions = reshaped_pure_functions %>% mutate(Pig_name = gsub("(\\D)*", "", Pig_name))

ggplot(reshaped_pure_functions, aes(Pig_name, COG_cat, fill= Relative_abundance)) +
  geom_tile(colour = "white", size = 0.5) +
  scale_fill_distiller(palette = "Spectral",
                       limits = c(min(reshaped_pure_functions$Relative_abundance), max(reshaped_pure_functions$Relative_abundance )),
                       direction = -1) +
  scale_y_discrete(expand=c(0, 0))+
  scale_x_discrete(expand=c(0, 0))+
  labs(fill="Relative Abundance", title = "Relative Abundance of functions by pig")+
  theme_grey(base_size=12)+
  theme(line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x.bottom = element_text(),
        axis.text.y = element_text(face="bold.italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border=element_blank(),
        legend.text=element_text(face="bold", size = 15),
        plot.title = element_text(face="bold", size = 25),
        legend.title = element_text(face="bold", size = 15),
        legend.key.size = unit(2, "cm"),
        strip.text.y = element_text(angle = 0, face = "bold", size = 15),
        strip.text.x = element_text(face = "bold", size = 15))+
  facet_grid(cols = vars(factor(reshaped_pure_functions$Category, levels = c("Control", "Colistin"))),
             scales = "free",
             space = "free",
             labeller = )
scale = 200