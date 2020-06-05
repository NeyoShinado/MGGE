
constructW <- function(X, NN, K){
  X = as.matrix(X)
  N = dim(X)[1]
  P = dim(X)[2]
  
  ## LPI algorithm chose npc
  if(min(N, P) <= 300){
    npc = min(N, P)
  }else{
    npc = 300
  }
  
  # normalize data
  cell_norms = matrix(sqrt(rowSums(X^2)), N, P)
  X = X / cell_norms
  # normalize genes not necessary
  #norm_gene = sqrt(colSums(X^2))
  #norm_gene = t(matrix(norm_gene, length(norm_gene), N))
  #X = X / norm_gene
  
  S = X %*% t(X)
  
  neighbors = list()
  neighbors$neighbors = t(sapply(1:N, function(i){
    nei = rep(0, N)
    nei[order(S[i, ], decreasing=TRUE)[1:(NN+1)]] = 1
    return(nei)
  }))
  
  neighbors$dist_list = lapply(1:N, function(i){
    S[i, ]
  })
  
  S = t(sapply(1:N, function(i){
    res = rep(0, N)
    id = order(S[i, ], decreasing=TRUE)[1:(NN+1)]
    res[id] = S[i, id]
    return(res)
  }))
  
  for(i in c(1:N)){
    for(j in c(i:N)){
      S[i, j] = max(S[i, j], S[j, i]) 
      S[j, i] = S[i, j]
    }
  }
  
  if(FALSE){
    S = X %*% t(X)
    S = sapply(1:N, function(i){
      temp = S[, i]
      id = order(temp, decreasing = TRUE)[1:NN]
      temp[setdiff(c(1:N), id)] = 0
      return(temp)
    })
    S = (S + t(S)) / 2
    D = Matrix(diag(rowSums(S)), sparse=TRUE)
    L = D - S
    
    X_mean = colSums(D %*% X) / sum(D)
    X = X - t(matrix(X_mean, P, N))
    
    # svd representive fea
    svd_res = svd(X %*% t(X))
    pro_X = svd_res$u[, c(1:npc)] %*% diag(sqrt(svd_res$d))
    
    # discriminated fea
    W_lpi = eigen(t(pro_X) %*% L %*% pro_X)$vectors[, (N-K):(N-1)]
    W = pro_X %*% W_lpi
    S = W %*% t(W)
    local_sim = S
  }
  
  
  if(FALSE){
    gene_weight = colSums(droprate)
    gene_weight = exp(-gene_weight)
    gene_weight = gene_weight / sum(gene_weight)
    norm_cell = rowSums(droprate) + 1e-10
    #* not use yet
    
    local_sim = t(sapply(1:N, function(i){
      valid_X = X * t(matrix(droprate[i, ], P, N))
      #local_sim = (valid_X[i, ] * gene_weight) %*% t(valid_X)
      local_sim = (valid_X[i, ] %*% t(valid_X))
      
      id = order(local_sim, decreasing = TRUE)[1:NN]
      #local_sim[id] = 1
      local_sim[setdiff(c(1:N), id)] = 0
      
      #local_sim = local_sim / sqrt(sum(local_sim^2)) #* F-norm. sum-norm
      return(local_sim)
    }))
    
    local_sim = (local_sim + t(local_sim)) / 2
    
  }
  
  res = list(S=S, neighbors=neighbors)
  return(res)
}


