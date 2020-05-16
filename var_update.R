# clustering vars updating
library(Matrix)
library(bignmf)


var_update <- function(lg_X, K, npc, lambda1=0.1, lambda2=0.1, 
                       thre_J=1e-4, drop_thre=0.5, iteration=5, clust_iteration=50, imp_iteration=100, output="result/", res_save=TRUE){
  message("## Initializing for program...\n")
  # create output_dir
    if(!dir.exists(paste0(output, "clust_res"))){
    dir.create(paste0(output, "clust_res"))
  }
  if(!dir.exists(paste0(output, "imp_res"))){
    dir.create(paste0(output, "imp_res"))
  }
  
  H = matrix(runif(dim(lg_X)[1]*K), N, K)
  imp_iter = 0
  
  
  ## Clustering update
  time_cum = system.time(
  while(imp_iter <= iteration){
    message(paste0("## running ", imp_iter, "th vars update via clustering and imputation...\n"))
    if(imp_iter == 0){
      clust_res = clustering_update(log2(10^lg_X - 1.01 + 1.0), K, npc, lambda1, lambda2, iteration=clust_iteration,
                                    out_file=paste0(output, "clust_res/", imp_iter, "th_clustres.rds"), res_save=res_save)
    }else{#* should not use imp_X
      clust_res = clustering_update(imp_X, K, npc, lambda1, lambda2, W=W, V=V, H=H, D=D, L=L, iteration=clust_iteration,
                                    out_file=paste0(output, "clust_res/", imp_iter, "th_clustres.rds"), res_save=res_save)
      H = clust_res$H
    }
    V = clust_res$V
    W = clust_res$W
    cluster = clust_res$cluster
    neighbors = clust_res$neighbors
    res_J = clust_res$J
  
      
    ## convergence check
    if(imp_iter == 0){
      J_set = res_J
    }else{
      var_J = res_J - J_set[length(J_set)]
      J_set = c(J_set, res_J)
    }
    
    if(exists("var_J")){
      if(abs(var_J) <= thre_J){
        res = list(cluster=cluster, neighbors=neighbors, cost_series=J_set,
                   imp_X=imp_X, droprate=droprate, local_sim=local_sim)
        saveRDS(res, paste0(output, "MGGE_res.rds"))
        return(res)
      }
    } 
    
    imp_iter = imp_iter + 1
    
    
    ## imputing update
    imp_res = imputing_update(lg_X, cluster, neighbors, imp_iteration=imp_iteration, out_file=paste0(output, "imp_res/", imp_iter, "th_impres.rds"), res_save=res_save)
    imp_X = imp_res$imp_X
    imp_X = log2(10^imp_X - 1.01 + 1.0)   #* log-tran 1.01 for imputing, 1.0 for clustering
    
    
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
    L = diag(rowSums(local_sim)) - local_sim
    
    
    # dropout weight based on drop_genes' number
    droprate = rowSums(imp_res$droprate)
    D = diag((P - droprate) / P)   #* dropout weight
    D = Matrix(D, sparse=TRUE)
  })
  res = list(cluster=cluster, neighbors=neighbors, cost_series=J_set,
             imp_X=imp_X, droprate=droprate, local_sim = local_sim, J_set=J_set, time_consume = time_cum)
  message("# MGGE iteration complete!\n")
  saveRDS(res, paste0(output, "MGGE_res_informative_gene.rds"))
  return(res)
}


clustering_update <- function(lg_X, K, npc, lambda1, lambda2, W=NULL, V=NULL, H=NULL, 
                              D=NULL, L=NULL, iteration=50, thre=1e-4, 
                              out_file="result/clust_res/res.rds", res_save=TRUE){
  message("## clustering vars updating...\n")
  
  
  # initialization
  J_set = NULL
  
  if(is.null(L)){
    # initialize H & W, V  
    message("## Pre-iteration Performing NMF on data...")
    res = bignmf(lg_X, npc)
    W = res$W
    V = res$H
    cluster = kmeans(W, K)$cluster
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
      id = union(order(dist_list[[i]])[1:round(N/10)], 
                 which(cluster==cluster[i]))
      nei[id] = 1
      return(nei)
    }))
    J = norm(lg_X - (W %*% V), type="F")^2
    
    res = list(W = W, V = V, cluster = cluster,
               neighbors = neighbors, J = J)
    if(res_save){
      saveRDS(res, out_file)
    }
    return(res)
  }else{
    iter = 1
    N = dim(W)[1]
    
    while(iter <= iteration){
      message(paste0("### Updating ", iter, "th round of ", iteration, "\n"))
      message(paste0(strrep("#", round(30*iter/iteration)), " ", 
                 100*signif(iter/iteration, digits=2), "%\n"))
      message("### updating W...\n")
      
      #*
      LH = H %*% t(H)
      sec_ord_g = lambda2 * LH %*% W + D^2 %*% W %*% V %*% t(V)
      g_W = -D^2 %*% lg_X %*% t(V) + sec_ord_g#* 0.36s
      ## the N*k * N*k Hessian matrix is corresponse to Wij which expand by col
      #* too large vector to expand if N & npc is big
      # Hessian = kronecker(Matrix(diag(1, npc), sparse=TRUE), Matrix(lambda2 * LH, sparse=TRUE)) + 
        # kronecker(V %*% t(V), Matrix(D^2, sparse=TRUE))  #* sparse 0.66s
      # sec_ord_g = Hessian %*% matrix(W, ncol=1)   
      
      W = W * ((D^2 %*% lg_X %*% t(V)) / sec_ord_g)   
      
      message("### updating V...\n")
      V = V * ((t(W) %*% D^2 %*% lg_X) / (t(W) %*% D^2 %*% W %*% V))
      
      message("### updating H...\n")
      #*
      d = sapply(1:N, function(i){
        norm(W[i, ], "2")^2
      })
      d = matrix(d, nrow=N, ncol=N)
      d = d + t(d) - 2*W %*% t(W)
      H = eigen(0.5*lambda1 * d + lambda2 * L)$vectors[, (N-K+1):N]
      
      
      # cost calculation
      # if verbose
      J = norm(D %*% (lg_X - W %*% V), 'F')^2 + lambda1 * sum(diag(t(W) %*% LH %*% W)) + 
        lambda2 * sum(diag(t(H) %*% L %*% H))
      message("### Current cost:", J, "\n")
      
      
      if(iter == 1){
        J_set = J
      }else{
        var_J = abs(J - J_set[length(J_set)])
        J_set = c(J_set, J)
        
        # convergence check
        if(var_J <= thre){
          cluster = kmeans(H, K)$cluster   #* apply(H, 1, which.max)
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
    
    cluster = kmeans(H, K)$cluster   #* apply(H, 1, which.max)
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
}


imputing_update <- function(lg_X, cluster, neighbors, drop_thre=0.5, out_file="result/imp_res/res.rds", ncores=1, imp_iteration=100, res_save=TRUE){
  # return log-tran X
  #* lg 1.01
  # not parallel version on imputation
  #* output file
  #* mix-model fit iteration
  imp_res = imputation_wlabel_model(lg_X, cluster, neighbors=neighbors, point=log2(1.01),
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

