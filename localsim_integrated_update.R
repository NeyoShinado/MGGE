# clustering vars updating
library(Matrix)
library(bignmf)


var_update <- function(lg_X, K, npc, lambda1=0.1, lambda2=0.1, 
                       thre_J=1e-4, drop_thre=0.5, iteration=1, 
                       clust_iteration=100, imp_iteration=100, 
                       output="result/", res_save=TRUE){
  # lg_X without pseudo count
  cat("## Initializing for program...\n")
  if(!dir.exists(paste0(output, "clust_res"))){
    dir.create(paste0(output, "clust_res"))
  }
  if(!dir.exists(paste0(output, "imp_res"))){
    dir.create(paste0(output, "imp_res"))
  }
  
  N = dim(lg_X)[1]
  P = dim(lg_X)[2]
  # generate cluster set & neighbors
  res = constructW(lg_X, 20)    # 20 nearest neighbors by default
  W = res$W
  neighbors = res$neighbors

  L = diag(colSums(t(W))) - W
  
  H = eigen(L)$vectors[, (N-K):(N-1)]
  H = Re(H)
  clust = kmeans(H, K)$cluster
  #clust = sapply(1:N, function(i){    #* not positive
  #  return(order(H[i, ], decreasing = TRUE)[1])
  #})
  
if(FALSE){
  pars = sapply(1:dim(lg_X)[2], function(i) {
    count_nzero = setdiff(lg_X[, i], 0)
    mu = mean(count_nzero)
    mu[is.na(mu)] = 0
    sd = sd(count_nzero)
    sd[is.na(sd)] = 0
    cv = sd/mu
    cv[is.na(cv)] = 0
    return(c(mu, cv))
  })
  
  mu = pars[1, ]
  cv = pars[2, ]
  high_var_genes = which(mu >= quantile(mu, 0.25) & cv >= quantile(cv, 0.25))
  if(length(high_var_genes) < 500){
    high_var_genes = 1:dim(lg_X)[2]}
  count_hv = lg_X[, high_var_genes]
  pca = prcomp(count_hv)
  eigs = (pca$sdev)^2
  var_cum = cumsum(eigs)/sum(eigs)
  pr_num = max(which.max(var_cum) > 0.4, K)
  mat_pcs = pca$x[, 1:pr_num]
  dist = sapply(1:N, function(i){
    d = sapply(1:i, function(j){
      sse = sum((mat_pcs[i, ] - mat_pcs[j, ])^2)
      sqrt(sse)
    })
    return(c(d, rep(0, N - i)))
  })
  dist = dist + t(dist)
  neighbors = list()
  neighbors$neighbors = sapply(1:N, function(i){
    nei = rep(0, N)
    nei[order(dist[i, ])[1:round(N/4)]] = 1
    return(nei)
  })
  neighbors$dist_list = lapply(1:N, function(i){
    dist[i, ]
  })    #* 
}
  
  #H = matrix(runif(dim(lg_X)[1]*K), N, K)
  
  imp_iter = 1
  while(imp_iter <= iteration){
    imp_res = imputing_update(log(exp(lg_X) + 0.01), clust, neighbors, imp_iteration=imp_iteration, out_file=paste0(output, "imp_res/", imp_iter, "th_impres.rds"), res_save=res_save)
    imp_X = imp_res$imp_X
    imp_X = log(exp(imp_X) - 0.01)   # log-tran 1.01 for imputing, 1.0 for clustering
    
    imp_X = Matrix(imp_X, sparse=TRUE)    
    
    # accumulate update with unitized local_sim
    if(imp_iter == 1){
      local_sim = imp_res$local_sim
    }else{
      local_sim = imp_res$local_sim + local_sim
      local_sim = t(apply(local_sim, 1, function(i){
        if(max(i) != min(i)){
          i = (i - min(i)) / (max(i) - min(i))
        }
        return(i)
      }))
    }
    local_sim = 0.5*(local_sim + t(local_sim))
    
    # dropout weight based on drop_genes' number
    if(imp_iter == 1){
      droprate = imp_res$droprate
    }else{
      droprate = imp_res$droprate + droprate
    }  
    D = rowSums(droprate)
    D = diag((P - D) / P)   #* dropout weight
    D = Matrix(D, sparse=TRUE)
    
    imp_iter = imp_iter + 1
  }
  
  D = Matrix(diag(1, N), sparse=TRUE)   #**
  ### Clustering update
  cat(paste0("## running ", imp_iter - 1, "th vars update via clustering and imputation...\n"))

  
  #imp_X = informative_gene_selection(imp_X, delta=0.7)   #* select informative genes for clustering
  

  H[H <= 0] = 1e-10
  W = matrix(runif(N * npc), N, npc) 
  V = matrix(runif(P * npc), npc, P) 
  clust_res = clustering_update(imp_X, K, npc, lambda1, lambda2, W=W, V=V, H=H, D=D, S=local_sim, iteration=clust_iteration,
                                out_file=paste0(output, "clust_res/localsim_integrated_clustres.rds"), res_save=res_save)
  #clust_res = clustering_update(imp_X, K, npc, lambda1, lambda2, W=W, V=V, H=H, D=D, L=L, iteration=clust_iteration,
   #                             out_file=paste0(output, "clust_res/localsim_integrated_clustres.rds"), res_save=res_save)
  
  H = clust_res$H
  V = clust_res$V
  W = clust_res$W
  cluster = clust_res$cluster
  neighbors = clust_res$neighbors
  res_J = clust_res$J
  
  
  res = list(cluster=cluster, neighbors=neighbors, 
             imp_X=imp_X, droprate=droprate, local_sim = local_sim, 
             J = res_J, W = W, V = V, H = H)
  cat("# MGGE iteration complete!\n")
  saveRDS(res, paste0(output, "MGGE_res.rds"))
  return(res)
}


clustering_update <- function(lg_X, K, npc, lambda1, lambda2, W=NULL, 
                              V=NULL, H=NULL, D=NULL, S=NULL, iteration=50, 
                              thre=1e-4, out_file="result/clust_res/res.rds", 
                              res_save=TRUE){
  cat("## clustering vars updating...\n")
  
  # initialization
  J_set = NULL


  iter = 1
  N = dim(W)[1]
  
  while(iter <= iteration){
    cat(paste0("### Updating ", iter, "th round of ", iteration, "\n"))
    cat(paste0(strrep("#", round(30*iter/iteration)), " ", 
               100*signif(iter/iteration, digits=2), "%\n"))
    cat("### updating W...\n")
    #*
    LH = H %*% t(H)
    sec_ord_g = lambda2 * LH %*% W + D^2 %*% W %*% V %*% t(V)
    g_W = -D^2 %*% lg_X %*% t(V) + sec_ord_g#* 0.36s
    
    ## the N*k * N*k Hessian matrix is corresponse to Wij which expand by col
    # Hessian = kronecker(Matrix(diag(1, npc), sparse=TRUE), Matrix(lambda2 * LH, sparse=TRUE)) + 
    #   kronecker(V %*% t(V), Matrix(D^2, sparse=TRUE))  #* sparse 0.66s
    #      Hessian = Matrix(Hessian, sparse=TRUE)   #* 10s
    # sec_ord_g = Hessian %*% matrix(W, ncol=1)   
    
    W = W * ((D^2 %*% lg_X %*% t(V)) / sec_ord_g)
    
    
    cat("### updating V...\n")
    V = V * ((t(W) %*% D^2 %*% lg_X) / (t(W) %*% D^2 %*% W %*% V))
    
    
    cat("### updating H...\n")
    #*
    d = colSums(W^2)
    d = matrix(d, nrow=N, ncol=N)
    d = d + t(d) - 2*W %*% t(W)
    #H = eigen(0.5*lambda1 * d + lambda2 * S)$vectors[, (N-K+1):N]
    D_s = Matrix(diag(colSums(S)), sparse = TRUE)
    H = H * ((lambda2 * S %*% H) / ((lambda1 / 2 * d + lambda2 * D_s) %*% H))
    #H[H <= 0] = 1e-10
    
    
    # cost calculation
    # if verbose
    J = norm(D %*% (lg_X - W %*% V), 'F')^2 + lambda1 * sum(diag(t(W) %*% LH %*% W)) + 
      lambda2 * sum(diag(t(H) %*% S %*% H))
    cat("### Current cost:", J, "\n")
    
    
    if(iter == 1){
      J_set = J
    }else{
      var_J = abs(J - J_set[length(J_set)])
      J_set = c(J_set, J)
      
      
      # convergence check
      if(var_J <= thre){
        cluster = sapply(1:N, function(i){
          which.max(H[i, ])
        })   #* apply(H, 1, which.max)
        
        dist = colSums(W^2)
        dist = matrix(d, nrow=N, ncol=N)
        dist = d + t(d) - 2*W %*% t(W)
        dist_list = lapply(1:N, function(i){
          dist[i, ]
        })
        
        neighbors = list()
        neighbors$dist_list = dist_list
        neighbors$neighbors = t(sapply(1:N, function(i){
          nei = rep(0, N)
          id = union(order(dist_list[[i]][1:round(N/4)]), 
                     which(cluster==cluster[i]))
          nei[id] = 1
          return(nei)
        }))
        res = list(W = W, V = V, H = H, cluster = cluster, 
                   neighbors = neighbors, J_set = J_set)
        if(res_save){
          saveRDS(res, out_file)
        }
        return(res)
      }
    }
    iter = iter + 1
  }
  
  #cluster = kmeans(H, K)$cluster
  cluster = sapply(1:N, function(i){
    which.max(H[i, ])
  })
  dist = sapply(1:N, function(i){
    d = sapply(1:i, function(j){
      dist = sum((W[i, ] - W[j, ])^2)   #* F-norm
      return(dist)
    })
    return(c(d, rep(0, N-i)))
  })
  
  dist = dist + t(dist)
  dist_list = lapply(1:N, function(i){
    dist[i, ]
  })
  neighbors = list()
  neighbors$dist_list = dist_list
  neighbors$neighbors = t(sapply(1:N, function(i){
    nei = rep(0, N)
    id = union(order(dist_list[[i]][1:round(N/10)]), 
               which(cluster==cluster[i]))
    nei[id] = 1
    return(nei)
  }))
  res = list(W = W, V = V, H = H, cluster = cluster, 
             neighbors = neighbors, J_set = J_set)
  if(res_save){
    saveRDS(res, out_file)
  }
  return(res)
  
}


imputing_update <- function(lg_X, cluster, neighbors, drop_thre=0.5, out_file="result/imp_res/res.rds", ncores=1, imp_iteration=100, res_save=TRUE){
  # return log-tran X
  #* lg 1.01
  # not parallel version on imputation
  #* output file
  #* mix-model fit iteration
  imp_res = imputation_wlabel_model(lg_X, cluster, neighbors=neighbors, point=log(1.01),
                                    drop_thre=drop_thre, ncores=ncores, imp_iteration=imp_iteration)
  imp_X = imp_res$count_imp
  droprate = imp_res$droprate
  local_sim = imp_res$local_sim
  if(res_save){
    saveRDS(imp_res, out_file)
  }
  rm(imp_res)
  
  return(list(imp_X = imp_X, local_sim = local_sim, droprate = droprate))
}

