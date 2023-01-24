
#lire tableau sorties Q3: les gènes les plus corrélés avec les métabolites
gene_cor_df = read.csv(file="refined_data/filtered_metabo_X_gene_cor.csv",
                   sep=";",
                   stringsAsFactors = TRUE,
                   nrow = 100000,
                   na.strings = c("NA", "-"))


#lire tableau des gènes
gene_df = read.csv(file="raw_data/Quantifications_and_functional_annotations.tsv",
                   sep="\t",
                   stringsAsFactors = TRUE,
                   #col.names = c("seed_cluster", "Preferred_name"),
                   nrows = 100000,
                   na.strings = c("NA", "-"))

#conserver les colonne d'interêt: contig et noms de gene
gene_df = subset(gene_df, select = c(seed_cluster , Preferred_name ) )
colnames(gene_df) = c("Contig", "Gene")

# conserver les lignes pour lesquelles le nom de gène se trouve dans gene_cor
gene_contig_filtered = subset(gene_df, Gene %in% gene_cor_df$Gene)

#suprimer les CDS pour ne garder que le nom du contig
splitted_CDS = strsplit(as.character(gene_contig_filtered$Contig), ".C")
gene_contig_filtered =  data.frame(t(data.frame(splitted_CDS)),gene_contig_filtered)
gene_contig_filtered =  data.frame(gene_contig_filtered[,1], gene_contig_filtered[,4])
colnames(gene_contig_filtered) = c("Contig", "Gene")

library(tidyr)
library(dplyr) 

# ajouter la colonne métabolite (duplique les ligne pour chaque metabolite corrélé)
gene_contig_meta = merge(gene_contig_filtered, gene_cor_df, group_by=  gene_contig_filtered$gene, all.x=TRUE)


# trier les lignes par ordre décroissant de facteur de correlation
ordered_gene_contig_meta <- gene_contig_meta[order(gene_contig_meta$Correlation, decreasing = T), ]




### RESEAU

ordered_gene_cor <- gene_cor_df[order(gene_cor_df$Correlation, decreasing = T), ]



ordered_gene_cor_filtered = ordered_gene_cor[1:100,]
nodes <- data.frame(
  name=c(unique(ordered_gene_cor_filtered[,1]), unique(ordered_gene_cor_filtered[,2])),
  carac=c( rep("metabolite",length(unique(ordered_gene_cor_filtered[,1]))), rep("gene",length(unique(ordered_gene_cor_filtered[,2]))))
)


library(igraph)

# create the network object
network <- graph_from_data_frame(d=ordered_gene_cor_filtered[,1:2],vertices=nodes, directed=F) 

# Make a palette of 3 colors
library(RColorBrewer)
coul  <- brewer.pal(3, "Set1") 

# Create a vector of color
my_color <- coul[as.numeric(as.factor(V(network)$carac))]


# plot it
png(file = "aaaa.png",
    width = 900 , height = 900,
  )
plot(network,
     edge.width= scale(ordered_gene_cor[,3], center = TRUE, scale = TRUE)+ 2,
     #layout=layout.circle,
     vertex.color = my_color )
dev.off()

scale = 200 # Set the scale for the exported image

ggsave(filename = "ksjfb.png")

savePlot(filename = "Network_100",
         type = "png"
       #width = 16*scale,
       #height = 18*scale,
       #units = "px",
       #dpi=200
)


bb = png(file = "Network_100.Rdata")
save(bb, file = "Network_100.png")
