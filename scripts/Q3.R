### ----------------------------------------
rm(list=ls())

library(tidyr)


df = read.csv(file="raw_data/Quantifications_and_functional_annotations.tsv",
              sep="\t",
              stringsAsFactors = TRUE,
              #nrows = 1000,
              na.strings = c("NA", "-"))
dim_raw = dim(df)[1]
df = df[,c(2:35,44)]
df = drop_na(df)
print(paste0((dim(df)[1]/dim_raw)*100, "% of the rows have been conserved (",dim(df)[1],")"))
head(df)
colnames(df) = c(paste0("Piglet_", seq(34)), "Gene")
df = df[c(dim(df)[2], 1:(dim(df)[2])-1)]

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

gene_names = levels(df$Gene)

df2 = as.data.frame(t(df))
colnames(df2) = df2[1,]
df2 = df2[-1,]
df2$Piglet = rownames(df2)

metabo = read.csv(file="raw_data/HOLOPIG METABOLITES clean.csv",
                  sep=";",
                  stringsAsFactors = TRUE,
                  na.strings = "NA")
metabo = subset(metabo, select = -c(Group))

metabo$Piglet = paste0("Piglet_",as.character(metabo$Piglet))
metabo_names = colnames(metabo[-1])

data = merge(df2, metabo, by = "Piglet", all = T)
data[,-1] = as.numeric(unlist(data[,-1]))
data = Filter(function(x)!all(is.na(x)), data) # Removing empty columns

#library(ellipse)
library(RColorBrewer)
cor_data = cor(data[,-1])
final_data = melt(cor_data)

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
final_data = melt(cor_data)

