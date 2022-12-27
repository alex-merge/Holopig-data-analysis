## A - Chargement du jeu de données

df = read.csv(file="raw_data/quantification_by_contig_lineage_all_reduced.csv",
              sep="\t")


## sommer les abondances relatives de toutes les ligne du même genre


rm(dat_genre, var, var2, test, abondance)
dat_genre= data.frame(0)

test=0
for (var in c(1:dim.data.frame(df)[1])){
  for (var2 in c(1:dim.data.frame(dat_genre)[2])){
    if (df[var,]$genre == dat_genre[1,var2]){
      dat_genre[2,var2]=as.integer(dat_genre[2,var2]) + as.integer(df[var,]$nb_reads_Holopig_colon01_quantif_percontig)  
      test=0
      break
    }
    else {test=1}
  }   
  if (test==1){
    dat_genre= cbind(dat_genre, c(df[var,]$genre , df[var,]$nb_reads_Holopig_colon01_quantif_percontig))
    colnames(dat_genre[var2])= df[var,]$genre
    test = 0
  }
}
colnames(dat_genre)= dat_genre[1,]
abondances= c(t(dat_genre[2,]))


barplot(t(abondances))


#barplot des règnes

rm(dat_regne, var, var2, test, abondances)
dat_regne= data.frame(0)

test=0
for (var in c(1:dim.data.frame(df)[1])){
  for (var2 in c(1:dim.data.frame(dat_regne)[2])){
    if (df[var,]$regne == dat_regne[1,var2]){
      dat_regne[2,var2]=as.integer(dat_regne[2,var2]) + as.integer(df[var,]$nb_reads_Holopig_colon01_quantif_percontig)  
      test=0
      break
    }
    else {test=1}
  }   
  if (test==1){
    dat_regne= cbind(dat_regne, c(df[var,]$regne , df[var,]$nb_reads_Holopig_colon01_quantif_percontig))
    colnames(dat_regne[var2])= df[var,]$regne
    test = 0
  }
}
colnames(dat_regne)= dat_regne[1,]
abondances= c(t(dat_regne[2,]))


barplot(t(abondances), names.arg = colnames(dat_regne))
abondances_hors_bacteries= abondances[-3]
barplot(t(abondances_hors_bacteries), names.arg = colnames(dat_regne[-3]))
