library(dplyr)
library(tidyr)

df = read.csv(file="raw_data/quantification_by_contig_lineage_all_reduced.csv",
              sep="\t",
              na.strings = "None")

df_abond = df[,1:7]
for (n in 1:34){
  if (n < 10){
    col_name=paste("nb_reads_Holopig_colon0",n,"_quantif_percontig", sep="")
  }
  else {
    col_name=paste("nb_reads_Holopig_colon",n,"_quantif_percontig", sep="")
  }
  df_abond[paste("abond",n,sep="")] = df[col_name]/sum(df[col_name])
}

write.table(df_abond, 
            file = "refined_data/taxonomic_abondance.csv",
            col.names = TRUE,
            row.names = FALSE,
            sep=";")

rm(list=ls())

df = read.csv(file="refined_data/taxonomic_abondance.csv",
              sep=";",
              stringsAsFactors = TRUE,
              na.strings = "NA")

genre = unique(df[["genre"]])
genre = genre[-1]

df_bygenre = data.frame(genre, row.names=genre)
for (name in genre){
  for (col in colnames(df)[8:41]){
    tmp = sum(subset(df, genre == name)[[col]])
    # print(tmp)
    df_bygenre[name,col] = tmp
  }
}

df_bp = expand.grid(genre=items, ind=paste0("ind",seq(1:34)))
values = c()
for (col in colnames(df_bygenre)[2:35]){
  values = c(values, df_bygenre[[col]])
}
df_bp$values = values

genreXbacterie = levels( droplevels( subset(df, regne == "Bacteria")$genre ) )[-1]
genreXarchees = levels( droplevels( subset(df, regne == "Archaea")$genre ) )[-1]
genreXeuk = levels( droplevels( subset(df, regne == "Eukaryota")$genre ) )[-1]
genreXvirus = levels( droplevels( subset(df, regne == "Viruses")$genre ) )[-1]

library(ggplot2)
library("RColorBrewer")
library(scales)
thres = -1
df_bp$values = -log(df_bp$values)
ggplot(subset(df_bp, values >= thres), aes(genre, ind, fill= values)) + 
  geom_tile()+
  theme_grey()

create_graph = function(genres, df=df_bp, thres=-1, log = F, title="Title"){
  data = subset(df_bp, genre %in% genres & values >= thres)
  if (log == T){
    print(data$values[1:20])
    data$values = -log(data$values)
  }
  ggplot(data, aes(ind, genre, fill= values)) + 
    geom_tile()+
    scale_fill_distiller(palette = "RdYlBu", 
                         limits = c(min(data$values), max(data$values))) +
    theme_grey()+
    xlab(label = "Individus")+
    ylab(label = "Genres")+
    scale_x_discrete(guide = guide_axis(check.overlap = TRUE))+
    labs(fill="Abondance relative", title = title)
    #theme(aspect.ratio = 9/21)
}
ggsave(filename = "export/HM_taxo.png",
       width = 4200,
       height = 1800,
       units = "px")
