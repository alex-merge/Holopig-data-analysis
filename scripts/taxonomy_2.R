library(dplyr)
library(tidyr)

#===============================================================================
# Refining abundance data : phyloXsample
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

#===============================================================================
# Refining abundance data : abundance by genre and sample.
rm(list=ls())

df = read.csv(file="refined_data/taxonomic_abondance.csv",
              sep=";",
              stringsAsFactors = TRUE,
              na.strings = "NA")

# Suppression des Règnes non identifiés
df = droplevels(subset(df, regne %in% c("Bacteria", "Archaea", "Eukaryota", "Viruses")))

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



df_bp = expand.grid(genre=genre, ind=as.character(seq(1:34)))
values = c()
for (col in colnames(df_bygenre)[2:35]){
  values = c(values, df_bygenre[[col]])
}
df_bp$values = values

bacteries = levels( droplevels( subset(df, regne == "Bacteria")$genre ) )
archees = levels( droplevels( subset(df, regne == "Archaea")$genre ) )
euka = levels( droplevels( subset(df, regne == "Eukaryota")$genre ) )
virus = levels( droplevels( subset(df, regne == "Viruses")$genre ) )
categorie = c("Control", "Colistine", "Control", "Control", "Control", "Control",
              "Colistine", "Colistine", "Control", "Colistine", "Colistine", "Control",
              "Control", "Control", "Control", "Control", "Control", "Colistine",
              "Colistine", "Colistine", "Colistine", "Colistine", "Control", "Control",
              "Control", "Control", "Control", "Colistine", "Control", "Colistine",
              "Colistine", "Colistine", "Colistine", "Colistine")

for (row in 1:dim(df_bp)[1]){
  if (df_bp[row, "genre"] %in% bacteries){df_bp[row, "regne"] = "Bacteria"}
  if (df_bp[row, "genre"] %in% archees){df_bp[row, "regne"] = "Archaea"}
  if (df_bp[row, "genre"] %in% euka){df_bp[row, "regne"] = "Eukaryota"}
  if (df_bp[row, "genre"] %in% virus){df_bp[row, "regne"] = "Viruses"}
  df_bp[row, "type"] = categorie[as.double(df_bp[row, "ind"])]
}

df_bp$type = as.factor(df_bp$type)
df_bp$regne = as.factor(df_bp$regne)
df_bp = df_bp[, c(4, 5, 1, 2, 3)]

write.table(df_bp, 
            file = "refined_data/taxonomic_abondance_by_genre.csv",
            col.names = TRUE,
            row.names = FALSE,
            sep=";")

#===============================================================================
# Taxo complète
rm(list = ls())
df = read.csv(file="refined_data/taxonomic_abondance.csv",
              sep=";",
              stringsAsFactors = TRUE,
              na.strings = "NA")

df$racine = as.factor(rep("racine",dim(df)[1]))
df = df[, c(dim(df)[2], 1:dim(df)[2]-1)]

base_taxo = df[, 1:8]
categorie = c("Controle", "Colistine", "Controle", "Controle", "Controle", "Controle",
              "Colistine", "Colistine", "Controle", "Colistine", "Colistine", "Controle",
              "Controle", "Controle", "Controle", "Controle", "Controle", "Colistine",
              "Colistine", "Colistine", "Colistine", "Colistine", "Controle", "Controle",
              "Controle", "Controle", "Controle", "Colistine", "Controle", "Colistine",
              "Colistine", "Colistine", "Colistine", "Colistine")

n_df = data.frame()
for (sample in 1:34){
  tmp = base_taxo
  tmp$ind = as.factor(rep(sample, dim(base_taxo)[1]))
  tmp$cat = as.factor(rep(categorie[sample], dim(base_taxo)[1]))
  tmp$abond = df[[8+sample]]
  n_df = rbind(n_df, tmp)
}

write.table(n_df, 
            file = "refined_data/taxonomic_abondance_complete.csv",
            col.names = TRUE,
            row.names = FALSE,
            sep=";")

#===============================================================================
# Creating the heatmap with abundance by genre with the complete csv file.
rm(list=ls())

library(ggplot2)
library(RColorBrewer)

df = read.csv(file="refined_data/taxonomic_abondance_complete.csv",
              sep=";",
              stringsAsFactors = TRUE,
              na.strings = "NA")
df$ind = as.factor(df$ind)

thres = 0.001
selected_classes = c("Archaea", "Bacteria", "Eukaryota", "Viruses")

# Filtering data and getting the log of the values.
data = subset(df, abond >= thres & regne %in% selected_classes)
data = drop_na(data[c("regne","genre","ind","cat","abond")])
data$abond = log(data$abond)

ggplot(data, aes(ind, genre, fill= abond)) + 
  geom_tile(colour = "white", size = 0.5)+
  scale_fill_distiller(palette = "Spectral", 
                       limits = c(min(data$abond), max(data$abond)),
                       direction = -1) +
  scale_y_discrete(expand=c(0, 0))+
  scale_x_discrete(expand=c(0, 0))+
  labs(fill="Logarithme de\nl'abondance relative", title = "Abondance relative des genres par échantillon")+
  theme_grey(base_size=12)+
  theme(line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
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
  facet_grid(rows = vars(data$regne), 
             cols = vars(factor(data$cat, levels = c("Controle", "Colistine"))),
             scales = "free",
             space = "free",
             labeller = )
scale = 200
ggsave(filename = "export/HM_taxo_complete(complete_csv).png",
       width = 16*scale,
       height = 18*scale,
       units = "px",
       dpi=200)


#===============================================================================
# Creating the heatmap with abundance by genre.

library(ggplot2)
library(RColorBrewer)
#library(scales)
thres = 0.001

# Filtering data and getting the log of the values.
data = subset(df_bp, values >= thres)
data$values = log(data$values)

ggplot(data, aes(ind, genre, fill= values)) + 
  geom_tile(colour = "white", size = 0.5)+
  scale_fill_distiller(palette = "Spectral", 
                       limits = c(min(data$values), max(data$values)),
                       direction = -1) +
  scale_y_discrete(expand=c(0, 0))+
  scale_x_discrete(expand=c(0, 0))+
  labs(fill="Logarithme de\nl'abondance relative", title = "Abondance relative des genres par échantillon")+
  theme_grey(base_size=12)+
  theme(line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
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
  facet_grid(rows = vars(data$regne), 
             cols = vars(factor(data$type, levels = c("Control", "Colistine"))),
             scales = "free",
             space = "free",
             labeller = )
scale = 200
ggsave(filename = "export/HM_taxo_complete.png",
       width = 16*scale,
       height = 18*scale,
       units = "px",
       dpi=200)

#===============================================================================
## Dendogram graph
library(tidyverse)

rm(list=ls())
thres = 0.001

df = read.csv(file="refined_data/taxonomic_abondance_complete.csv",
              sep=";",
              stringsAsFactors = T,
              na.strings = "NA")
df$ind = as.factor(df$ind)

data = subset(df, ind == "1")[,1:8]

values_t = rep(0, dim(data)[1])
values_c = rep(0, dim(data)[1])
for (i in levels(df$ind)){
  if (subset(df, ind == i)$cat[1] == "Controle"){
    values_t = values_t + subset(df, ind == i)$abond
  }
  else{
    values_c = values_c + subset(df, ind == i)$abond
  }
}
data$mean_temoin = values_t
data$mean_coli = values_c
data = subset(data, mean_temoin >= thres & mean_coli >= thres)

data_tmp = data[, 1:7]
#data = drop_na(data)
colnames(data_tmp)=paste0("level",seq(7))

edges_level1_2 <- data_tmp %>% select(level1, level2) %>% unique %>% rename(from=level1, to=level2)
edges_level2_3 <- data_tmp %>% select(level2, level3) %>% unique %>% rename(from=level2, to=level3)
edges_level3_4 <- data_tmp %>% select(level3, level4) %>% unique %>% rename(from=level3, to=level4)
edges_level4_5 <- data_tmp %>% select(level4, level5) %>% unique %>% rename(from=level4, to=level5)
edges_level5_6 <- data_tmp %>% select(level5, level6) %>% unique %>% rename(from=level5, to=level6)
edges_level6_7 <- data_tmp %>% select(level6, level7) %>% unique %>% rename(from=level6, to=level7)
edges=drop_na(rbind(edges_level1_2, edges_level2_3, edges_level3_4, edges_level4_5,
                    edges_level5_6, edges_level6_7)
              )

unwanted_cat = c("", "Unable to found taxonomy consensus", "Unknown")
edges = subset(edges, !(from %in% unwanted_cat) & !(to %in% unwanted_cat))
edges = droplevels(edges)


values_t = c()
values_c = c()
for (name in unique(c(as.character(edges$from), as.character(edges$to)))){
  if (name %in% levels(data$genre)){
    values_t = c(values_t, mean(subset(data, genre == name)$mean_t))
    values_c = c(values_c, mean(subset(data, genre == name)$mean_c))
  }else {
    values_t = c(values_t, NA)
    values_c = c(values_c, NA)
  }
}

df_t = data.frame(name = unique(c(as.character(edges$from), as.character(edges$to))), 
                  abond = values_t)
df_t$group = edges$from[ match( df_t$name, edges$to ) ]
df_c = data.frame(name = unique(c(as.character(edges$from), as.character(edges$to))), 
                  abond = values_c)
df_c$group = edges$from[ match( df_c$name, edges$to ) ]

library(ggraph)
library(igraph)
library(RColorBrewer) 

t_graph <- graph_from_data_frame( edges, vertices=df_t )

# Make the plot
ggraph(t_graph, layout = 'dendrogram', circular = TRUE) + 
  geom_edge_diagonal(colour="grey") +
  scale_edge_colour_distiller(palette = "RdPu") +
  geom_node_text(aes(x = x*1.15, y=y*1.15, filter = leaf, label=name, colour=group), size=2.7, alpha=1) +
  geom_node_point(aes(filter = leaf, x = x*1.07, y=y*1.07, colour=group, size=abond, alpha=0.2)) +
  scale_colour_manual(values= rep( brewer.pal(9,"Paired") , 30)) +
  scale_size_continuous( range = c(0.1,10) ) +
  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(c(0,0,0,0),"cm"),
  ) +
  expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3))

#===============================================================================
## Barplot graph
library(ggplot2)

rm(list=ls())

df = read.csv(file="refined_data/taxonomic_abondance_complete.csv",
              sep=";",
              stringsAsFactors = T,
              na.strings = "NA")
df$ind = as.factor(df$ind)

thres = 0.001
selected_classes = c("Archaea", "Bacteria", "Eukaryota", "Viruses")

# Filtering data and getting the log of the values.
data = subset(df, abond >= thres & regne %in% selected_classes)
data = drop_na(data[c("regne","genre","ind","cat","abond")])
data$abond = log(data$abond)
data = droplevels(data[order(data$abond, decreasing = T),])
data$genre = factor(data$genre, levels = unique(data$genre))

ggplot(data, aes(abond, genre)) + 
  geom_boxplot(fill="slateblue", alpha=0.2)+
  scale_y_discrete(expand=c(0, 0))+
  #scale_x_discrete(expand=c(0, 0))+
  labs(fill="Logarithme de\nl'abondance relative", title = "Abondance relative des genres par échantillon")+
  xlab("Logarithme de l'abondance relative")+
  theme_grey(base_size=12)+
  theme(line = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(face="bold"),
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
  facet_grid(rows = vars(data$regne), 
             cols = vars(factor(data$cat, levels = c("Controle", "Colistine"))),
             scales = "free",
             space = "free",
             labeller = )
scale = 200
ggsave(filename = "export/BP_taxo_complete(complete_csv).png",
       width = 16*scale,
       height = 18*scale,
       units = "px",
       dpi=200)
