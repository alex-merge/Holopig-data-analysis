df = read.csv(file="raw_data/quantification_by_contig_lineage_all_reduced.csv",
              sep="\t")

## Source : https://stackoverflow.com/questions/15343338/how-to-convert-a-data-frame-to-tree-structure-object-such-as-dendrogram
## recursion function
traverse <- function(a,i,innerl){
  if(i < (ncol(df))){
    alevelinner <- as.character(unique(df[which(as.character(df[,i])==a),i+1]))
    desc <- NULL
    for(b in alevelinner) desc <- c(desc,traverse(b,i+1,innerl))
    il <- NULL; if(innerl==TRUE) il <- a
    (newickout <- paste("(",paste(desc,collapse=","),")",il,sep=""))
  }
  else { (newickout <- a) }
}

## data.frame to newick function
df2newick <- function(df, innerlabel=FALSE){
  alevel <- as.character(unique(df[,1]))
  newick <- NULL
  for(x in alevel) newick <- c(newick,traverse(x,1,innerlabel))
  (newick <- paste("(",paste(newick,collapse=","),");",sep=""))
}

tree = df2newick(df = df[1:2,c("regne", "phylum", "classe")])

## see: https://phylot.biobyte.de/