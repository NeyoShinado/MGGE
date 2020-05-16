
# input data matrix of N*P size
informative_gene_selection <- function(data, delta=0.8){

  # normalize expr's scale
  counts = data
  data = sweep(data, 2, colMeans(data), '-')
  
  # construce P*P kernel matrix of genes
  P = dim(data)[2]
  kernel = t(data) %*% data
  res = eigen(kernel, symmetric=TRUE)
  lambda = res$values
  vec = res$vectors
  
  # building genes' weight
  weight = vec %*% matrix(lambda, P, 1)
  
  rank_id = order(weight, decreasing=TRUE)
  weight = sort(weight, decreasing=TRUE)
  weight = cumsum(weight) / sum(weight)
  
  id = which(weight >= delta)[1]
  select_gene = rank_id[1:id]
  message("## selecting ", length(select_gene), " informative genes for downstream analysis...\n")
  
  select_data = counts[, select_gene]
  return(select_data)
  
}