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
library(tidyr)


df = read.csv(file="raw_data/Quantifications_and_functional_annotations.tsv",
              sep="\t",
              stringsAsFactors = TRUE,
              na.strings = "NA",
              nrows = 1000)
df = df[,1:35]
head(df)
colnames(df) = c("Gene", paste0("Piglet_", seq(34)))
for (i in 2:35){
  df[,i] = df[,i]/sum(df[,i])
}
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
data = merge(df2, metabo, by = "Piglet", all = T)
data[,-1] = as.numeric(unlist(data[,-1]))

library(ellipse)
library(RColorBrewer)
cor_data = cor(data[,-1])

# Build a Pannel of 100 colors with Rcolor Brewer
my_colors <- brewer.pal(5, "Spectral")
my_colors <- colorRampPalette(my_colors)(100)

ord <- order(cor_data[1, ])
data_ord <- cor_data[ord, ord]
plotcorr(data_ord , col=my_colors[data_ord*50+50] , mar=c(1,1,1,1)  )
