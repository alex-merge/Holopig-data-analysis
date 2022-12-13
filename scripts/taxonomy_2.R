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

bacteries = levels( droplevels( subset(df, regne == "Bacteria")$genre ) )[-1]
archees = levels( droplevels( subset(df, regne == "Archaea")$genre ) )[-1]
euka = levels( droplevels( subset(df, regne == "Eukaryota")$genre ) )[-1]
virus = levels( droplevels( subset(df, regne == "Viruses")$genre ) )[-1]
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


library(ggplot2)
library(RColorBrewer)
library(scales)
thres = 0.001
# df_bp$values = -log(df_bp$values)



# Heatmap des archées
data = subset(df_bp, genre %in% archees & values >= thres)
ggplot(data, aes(ind, genre, fill= values)) + 
  geom_tile(colour = "white", size = 0.7)+
  scale_fill_distiller(palette = "Spectral", 
                       limits = c(min(data$values), max(data$values)),
                       direction = -1) +
  scale_y_discrete(expand=c(0, 0))+
  scale_x_discrete(expand=c(0, 0))+
  labs(fill="Abondance relative", title = "Abondance relative des archées par échantillon")+
  theme_grey(base_size=10)+
  theme(line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border=element_blank(),
        legend.text=element_text(face="bold"),
        plot.title = element_text(face="bold"),
        legend.title = element_text(face="bold"))
scale = 150
ggsave(filename = "export/HM_taxo_archees.png",
       width = 21*scale,
       height = 9*scale,
       units = "px",
       dpi=200)

# Heatmap des bactéries
data = subset(df_bp, genre %in% bacteries & values >= thres)
ggplot(data, aes(ind, genre, fill= values)) + 
  geom_tile(colour = "white", size = 0.5)+
  scale_fill_distiller(palette = "Spectral", 
                       limits = c(min(data$values), max(data$values)),
                       direction = -1) +
  scale_y_discrete(expand=c(0, 0))+
  scale_x_discrete(expand=c(0, 0))+
  labs(fill="Abondance relative", title = "Abondance relative des bactéries par échantillon")+
  theme_grey(base_size=10)+
  theme(line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border=element_blank(),
        legend.text=element_text(face="bold"),
        plot.title = element_text(face="bold"),
        legend.title = element_text(face="bold"))
scale = 200
ggsave(filename = "export/HM_taxo_bacteries.png",
       width = 20*scale,
       height = 11*scale,
       units = "px",
       dpi=200)

# Heatmap totale
data = subset(df_bp, values >= thres)
ggplot(data, aes(ind, genre, fill= values)) + 
  geom_tile(colour = "white", size = 0.5)+
  scale_fill_distiller(palette = "Spectral", 
                       limits = c(min(data$values), max(data$values)),
                       direction = -1) +
  scale_y_discrete(expand=c(0, 0))+
  scale_x_discrete(expand=c(0, 0))+
  labs(fill="Abondance relative", title = "Abondance relative des genres par échantillon")+
  theme_grey(base_size=12)+
  theme(line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
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
