# Reading data of mouse embryo
setwd("D:/Work/Project/dataset/sc-RNA/mouse embryo/")
files = list.files(".", pattern="*.txt$")

data = lapply(files, function(i){
  data = read.table(i, head=TRUE, comment.char="")
  list(reads=data$reads, RPKM=data$RPKM, gene=data$X.Gene_symbol)
})
