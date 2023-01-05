@@ -0,0 +1,94 @@
  # Load ggplot2
  library(ggplot2)


## A - Chargement du jeu de données

df = read.csv(file="refined_data/taxonomic_abondance.csv",
              sep=";")

moyenne = data.frame(mean = mean(c(t(df[1,-(1:7)]))))
for (var in c(2:dim.data.frame(df)[1])){
  moyenne = rbind (moyenne, mean(c(t(df[var,-(1:7)]))))
}
df = cbind(df, moyenne)
## sommer les abondances relatives de toutes les ligne du même genre

rm(dat_genre, var, var2, test, abondance)
dat_genre= data.frame("NA" = 0)

test=0
for (var in c(1:dim.data.frame(df)[1])){
  if (is.na(df[var,]$genre)){
    dat_genre[2,1]=as.integer(dat_genre[2,1]) + df[var,]$mean
    test=0
    next
  }
  for (var2 in c(1:dim.data.frame(dat_genre)[2])){
    if (df[var,]$genre == dat_genre[1,var2]){
      dat_genre[2,var2]=as.integer(dat_genre[2,var2]) + as.integer(df[var,]$mean)  
      test=0
      break
    }
    else {test=1}
  }   
  if (test==1){
    dat_genre= cbind(dat_genre, c(df[var,]$genre , df[var,]$mean))
    colnames(dat_genre[var2])= df[var,]$genre
    test = 0
  }
}

#FAIRE LES GRAPHES
# Create data
data <- data.frame(
  genre=c(t(dat_genre[1,])) ,  
  value=c(t(dat_genre[2,]))
)

# Barplot
ggplot(data, aes(x=genre, y=value)) + 
  geom_bar(stat = "identity") +
  coord_flip()




#BARPLOT DES REGNES

rm(dat_regne, var, var2, test, abondances)
dat_regne= data.frame("NA" = 0)

test=0
for (var in c(1:dim.data.frame(df)[1])){
  if (is.na(df[var,]$regne)){
    dat_regne[2,1]=as.integer(dat_regne[2,1]) + df[var,]$mean
    test=0
    next
  }
  for (var2 in c(1:dim.data.frame(dat_regne)[2])){
    if (df[var,]$regne == dat_regne[1,var2]){
      dat_regne[2,var2]=as.integer(dat_regne[2,var2]) + df[var,]$mean  
      test=0
      break
    }
    else {test=1}
  }   
  if (test==1){
    dat_regne= cbind(dat_regne, c(df[var,]$regne , df[var,]$mean))
    colnames(dat_regne[var2])= df[var,]$regne
    test = 0
  }
}

#FAIRE LES GRAPHES
# Create data
data <- data.frame(
  regne=c(t(dat_regne[1,])) ,  
  value=c(t(dat_regne[2,]))
)

# Barplot
ggplot(data, aes(x=regne, y=value)) + 
  geom_bar(stat = "identity") +
  coord_flip()