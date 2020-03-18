library(purrr)


# Reading data of mouse embryo
setwd("D:/Work/Project/dataset/sc-RNA/mouse embryo/")
files = list.files(".", pattern="*.txt$")
cell_ID = sapply(files, function(i){substr(i, 1, 10)})

data = lapply(files, function(i){
  data = read.table(i, head=TRUE, comment.char="")
  list(reads=data$reads, RPKM=data$RPKM, gene=data$X.Gene_symbol)
})


# list transpose
result = transpose(data)
result$gene = as.character(result$gene[[1]])
result$cell_ID = as.character(cell_ID)
result$reads = data.frame(Reduce(rbind, result$reads))
result$RPKM = data.frame(Reduce(rbind, result$RPKM))
names(result$reads) = result$gene
names(result$RPKM) = result$gene
row.names(result$reads) = cell_ID
row.names(result$RPKM) = cell_ID

if(FALSE){
  reads = data.frame(matrix(unlist(lapply(data, function(i){
    i$reads
  })), nrow=length(data), byrow=TRUE))
  names(reads) = data[[1]]$gene
  
  RPKM = data.frame(matrix(unlist(lapply(data, function(i){
    i$RPKM
  })), nrow=length(data), byrow=TRUE))
  names(RPKM) = data[[1]]$gene
}


# loading labels
series = read.table("GSE45719_series_matrix.txt", head=F)
fea_name = series[, 1]
series = data.frame(t(series[, 2:ncol(series)]))
names(series) = fea_name
series = data.frame(ID=series$ID_REF, 
                    source_name=series$`!Sample_source_name_ch1`, 
                    type=series$`!Sample_description`)

strvec_split = function(strvec, n){
  if(n == 1){
    sapply(as.character(strvec), function(i){
      strsplit(i, " ")[[1]][n]})
  }else if(n > 1){
    sapply(as.character(strvec), function(i){
      do.call(paste, as.list(c(strsplit(i, " ")[[1]][1:n], sep=" ")))
      })
  }else if(!is.integer(n) || n < 1){
    stop("n should be inter and should greater than 1!")
  }
}

# labels edit
result$label[1:100] = strvec_split(series$source_name[1:100], 1)
result$label[101:105] = as.character(series$source_name[101:105])
result$label[106:113] = strvec_split(series$source_name[106:113], 1)
result$label[114:276] = strvec_split(series$source_name[114:276], 2)
result$label[277:280] = as.character(series$source_name[277:280])
result$label[281:317] = strvec_split(series$source_name[281:317], 1)
result$label = unlist(result$label)

# single cell mark
result$singlecell[1:92] = TRUE
result$singlecell[101:280] = TRUE
result$singlecell[289:288] = TRUE

# save data
saveRDS(result, file="./mouse_embryo.rds")