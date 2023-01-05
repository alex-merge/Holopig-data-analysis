
# Loading raw data as a data.frame called raw_df
raw_df = read.csv(file="raw_data/quantification_by_contig_lineage_all_reduced.csv",
              sep="\t",
              na.strings = "None")

# Selecting 7th first columns of the raw_df and calling it abundance_df
abundance_df = df[,1:7]

# 
for (n in 1:34){
  if (n < 10){
    col_name=paste("nb_reads_Holopig_colon0",n,"_quantif_percontig", sep="")
  }
  else {
    col_name=paste("nb_reads_Holopig_colon",n,"_quantif_percontig", sep="")
  }
  abundance_df[paste("abond",n,sep="")] = df[col_name]/sum(df[col_name])
}

write.table(df_abond, 
            file = "refined_data/taxonomic_abondance.csv",
            col.names = TRUE,
            row.names = FALSE,
            sep=";")

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