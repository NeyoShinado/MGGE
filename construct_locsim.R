
construct_locsim <- function(X, droprate, NN){
  X = as.matrix(X)
  N = dim(X)[1]
  P = dim(X)[2]
  
  
  gene_weight = colSums(droprate)
  gene_weight = exp(-gene_weight)
  gene_weight = gene_weight / sum(gene_weight)
  norm_cell = rowSums(droprate) + 1e-10
  #* not use yet
  
  local_sim = t(sapply(1:N, function(i){
    valid_X = X * t(matrix(droprate[i, ], P, N))
    local_sim = (valid_X[i, ] * gene_weight) %*% t(valid_X)
    
    id = order(local_sim, decreasing = TRUE)[1:NN]
    local_sim[id] = 1
    local_sim[setdiff(c(1:N), id)] = 0
    
    #local_sim = local_sim / sqrt(sum(local_sim^2)) #* F-norm. sum-norm
    return(local_sim)
  }))

  local_sim = (local_sim + t(local_sim)) / 2
  
  return(local_sim)
}


