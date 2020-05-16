# construct neighobr vec & dist_mat for data


constructW <- function(X, k){
  
  N = dim(X)[1]
  
  # normalize
  norm_x = sqrt(colSums(t(X * X)))
  X = X / norm_x
  
  W = X %*% t(X)
  dist = -W
  
  neighbors = list()
  neighbors$neighbors = t(sapply(1:N, function(i){
    nei = rep(0, N)
    nei[order(dist[i, ])[1:(k+1)]] = 1
    return(nei)
  }))
  
  neighbors$dist_list = lapply(1:N, function(i){
    dist[i, ]
  })
  
  sparse_W = t(sapply(1:N, function(i){
    res = rep(0, N)
    id = order(W[i, ], decreasing = TRUE)[1:(k+1)]
    #res[id] = W[i, id]
    res[id] = 1
    
    return(res)
  }))
  
  res = list(neighbors = neighbors, W = sparse_W)
  
  return(res)
}