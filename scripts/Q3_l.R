# Libraries
library(ellipse)
library(RColorBrewer)

# Use of the mtcars data proposed by R
data <- cor(mtcars)

# Build a Pannel of 100 colors with Rcolor Brewer
my_colors <- brewer.pal(5, "Spectral")
my_colors <- colorRampPalette(my_colors)(100)

# Order the correlation matrix
ord <- order(data[1, ])
data_ord <- data[ord, ord]
plotcorr(data_ord , col=my_colors[data_ord*50+50] , mar=c(1,1,1,1)  )

### ----------------------------------------
rm(list=ls())

library(tidyr)


df = read.csv(file="refined_data/Annot_1.csv",
              sep=";",
              stringsAsFactors = TRUE,
              nrows = 1000,
              na.strings = c("NA", "-"))
dim_raw = dim(df)[1]
df = df[,c(2:35,44)]
df = drop_na(df)
print(paste0((dim(df)[1]/dim_raw)*100, "% of the rows have been conserved (",dim(df)[1],")"))
head(df)
colnames(df) = c(paste0("Piglet_", seq(34)), "Gene")
df = df[c(dim(df)[2], 1:(dim(df)[2])-1)] #pour mettre la colonne gène à gauche. 

# Computing relative abundance
for (i in 2:35){
  sum_col = sum(df[,i])
  if (sum_col != 0){
    df[,i] = df[,i]/sum_col
  }
}

# Aggregating 
df = aggregate(df[-1], by = list(df$Gene), FUN = sum)
colnames(df)[1] = "Gene"

gene_names = levels(df$Gene) #on stocke la lite des genes

df2 = as.data.frame(t(df)) #on transpose df (indiv en ligne et genes en colonne)
colnames(df2) = df2[1,]
df2 = df2[-1,]
df2$Piglet = rownames(df2)

metabo = read.csv(file="raw_data/HOLOPIG METABOLITES clean.csv",
                  sep=";",
                  stringsAsFactors = TRUE,
                  na.strings = "NA")
metabo = subset(metabo, select = -c(Group)) #supression de la colonne group

metabo$Piglet = paste0("Piglet_",as.character(metabo$Piglet))
metabo_names = colnames(metabo[-1]) #on stocke la liste des métabolites

data = merge(df2, metabo, by = "Piglet", all = T)
data[,-1] = as.numeric(unlist(data[,-1]))
data = Filter(function(x)!all(is.na(x)), data) # Removing empty columns

#library(ellipse)
library(RColorBrewer)
cor_data = cor(data[,-1])
#final_data = melt(cor_data)

cor_data = subset(cor_data, rownames(cor_data) %in% metabo_names, select = as.factor(gene_names))
#thres = 0.95
#cor_data[cor_data[] < thres] = NA

col_ord = order( colMeans(cor_data) )
row_ord = order( rowMeans(cor_data) )
#ord_data <- cor_data[ord, ord]


heatmap(cor_data[row_ord,col_ord], Colv = NA, Rowv = NA)

#plotcorr(data_ord , col=my_colors[data_ord*50+50] , mar=c(1,1,1,1)  )

library(reshape2)
cor_data = as.data.frame(cor_data)
cor_data$metabolites <- rownames(cor_data)
cor_data = cor_data[c(dim(cor_data)[2], 1:(dim(cor_data)[2])-1)] #pour mettre la colonne gène à gauche. 
#cor_data$mean = rowMeans(cor_data[,-1])


final_data = melt(cor_data)
final_data$metabolites = as.factor(final_data$metabolites)

#calculating the average correlation coefficient of each metabolite
mean_per_metab = data.frame()
for (i in levels(final_data$metabolites)){
mean_per_metab = rbind(mean_per_metab,c(i, colMeans(subset(final_data, metabolites == i, select = value))))
}
colnames(mean_per_metab) = c("metabolite", "mean")

# keeping only the metabolites most correlated to the genes 
mean_per_metab = subset(mean_per_metab, mean >= 0.1)
final_data2= subset(final_data, metabolites %in% mean_per_metab$metabolite)


#CREATION HEATMAP
library(ggplot2)
library(RColorBrewer)
library(tidyr)

ggplot(final_data2, aes(variable, metabolites, fill= value)) +
  geom_tile(colour = "white", size = 0.5)+
  scale_fill_gradientn(name = "Correlation coefficient",# definition of the color gradient.
                       #trans = "log10", #chosing a logarythmic color scale
                       breaks = c(-1, -0.5, 0, 0.5 ,1), 
                       labels = c(-1, -0.5, 0, 0.5 ,1),
                       #na.value = "#43A3B9",# plotting the 0 in the lowest color of the gradient
                       colours = c("#43A3B9",  "#D0EC97", "#FFFFDF", "#FD945D", "#DC4E51", "#370303"),
  )+
  
  scale_x_discrete(expand=c(0, 0))+ # Removing blank space between data and X axis
  labs(fill="correlation coefficient", # Setting legend and title name
       title = "Correlation coefficients between genes and metabolites")+
  theme_grey(base_size=12)+ # Setting the base font size
  theme(line = element_blank(), # Removing the lines in between the boxes
        axis.ticks.x = element_blank(), # Removing ticks on X axis
        axis.text.x = element_text(angle = 90, face="bold.italic", size = 30), # Removing X labels
        axis.text.y = element_text(face="bold.italic", size = 30), # Setting Y labels to bold
        axis.title.x = element_blank(), # Removing X axis title
        axis.title.y = element_blank(), # Removing Y axis title
        panel.border=element_blank(), # Removing space between figure and panel
        legend.text=element_text(face="bold", size = 15), # Customizing legend font
        plot.title = element_text(face="bold", size = 50), # Same with main title
        legend.title = element_text(face="bold", size = 15), # Same with legend
        legend.key.size = unit(2, "cm"), # Setting the legend scale
        strip.text.y = element_text(angle = 0, face = "bold", size = 15),
        strip.text.x = element_text( face = "bold", size = 15))
        
  #facet_grid(rows = vars(filtered_data$kingdom), 
             #cols = vars(factor(filtered_data$cat, 
                                #levels = c("Control", "Colistin"))),
             #scales = "free",
             #space = "free")

#SAVE HEATMAP  
scale = 200 # Set the scale for the exported image

ggsave(filename = "export/heatmap_correlation.png",
       plot = last_plot(),
       width = 40*scale,
       height = 30*scale,
       units = "px",
       dpi=200)

