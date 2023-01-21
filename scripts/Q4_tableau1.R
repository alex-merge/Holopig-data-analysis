
#lire tableau sorties Q3: les gènes les plus corrélés avec les métabolites
gene_cor_df = read.csv(file="refined_data/filtered_metabo_X_gene_cor.csv",
                   sep=";",
                   stringsAsFactors = TRUE,
                   nrow = 100000,
                   na.strings = c("NA", "-"))
colnames(gene_cor_df) =c("Metabolite",  "Gene", "Correlation")


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





library(tidyr)
library(dplyr) 

# ajouter la colonne métabolite (duplique les ligne pour chaque metabolite corrélé)
gene_contig_meta = merge(gene_contig_filtered, gene_cor_df, group_by=  gene_contig_filtered$gene, all.x=TRUE)




# trier les lignes par ordre décroissant de facteur de correlation
ordered_gene_contig_meta <- gene_contig_meta[order(gene_contig_meta$Correlation, decreasing = T), ]















