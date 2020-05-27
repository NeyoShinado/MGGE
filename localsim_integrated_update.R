# clustering vars updating
library(Matrix)
library(bignmf)


var_update <- function(lg_X, K, npc, gene_id, M, lambda1=0.1, lambda2=0.1, 
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
    # generate cluster set & neighbors
    NN = as.integer(log10(N))
    NN = min(50, max(10, NN))
    res = constructW(lg_X, NN)
    W = res$W
    neighbors = res$neighbors
  
    L = diag(colSums(t(W))) - W
    
    H = eigen(L)$vectors[, (N-K):(N-1)]
    H = Re(H)
    #* 4\guide clust for imp
    clust = kmeans(H, K)$cluster
    #clust = rep(1, N)
    
    #clust = sapply(1:N, function(i){    
    #  return(order(H[i, ], decreasing = TRUE)[1])
    #})
    
    
    #H = matrix(runif(dim(lg_X)[1]*K), N, K)
    
    imp_iter = 1
    local_sim = list()
    for(i in c(1:M)){
    #- par-iteration is unnessary for MGGE any more
    cat(paste0("## carrying ", i, "th local_sim constructing...\n"))
    id = (1+(i-1)*1000):min(length(gene_id), i*1000)
    imp_res = imputing_update(log(exp(lg_X[, id]) + 0.01), clust, neighbors, 
                              imp_iteration=imp_iteration, 
                              out_file=paste0(output, "imp_res/", imp_iter, "th_impres.rds"), 
                              res_save=res_save)
    imp_X = imp_res$imp_X
    #imp_X = log(exp(imp_X) - 0.01)   # log-tran 1.01 for imputing, 1.0 for clustering
    #imp_X = Matrix(imp_X, sparse=TRUE)    
    
    # accumulate update with normalized local_sim
    #local_sim = imp_res$local_sim
    
    
    # dropout weight based on drop_genes' number
    droprate = imp_res$droprate


  #** 3\test of best local_sim
  local_sim[[i]] = construct_locsim(imp_X, droprate, NN)
  }
  
  
  D = Matrix(diag(1, N), sparse=TRUE)   #** dropout weight
  ### Clustering update
  cat(paste0("## running ", imp_iter - 1, "th vars update via clustering and imputation...\n"))

  
  #imp_X = informative_gene_selection(imp_X, delta=0.7)   #- select informative genes for clustering
  
  P = length(id)
  H[H <= 0] = 1e-10
  W = matrix(runif(N * npc), N, npc) 
  V = matrix(runif(P * npc), npc, P) 
  #* 7\ using imp_X
  clust_res = clustering_update(lg_X[, id], K, npc, lambda1, lambda2, W=W, V=V, 
                                H=H, D=D, local_sim=local_sim, iteration=clust_iteration,
                                out_file=paste0(output, "clust_res/localsim_integrated_clustres.rds"), 
                                res_save=res_save)
  #clust_res = clustering_update(imp_X, K, npc, lambda1, lambda2, W=W, V=V, H=H, D=D, L=L, iteration=clust_iteration,
   #                             out_file=paste0(output, "clust_res/localsim_integrated_clustres.rds"), res_save=res_save)
  
  H = clust_res$H
  V = clust_res$V
  W = clust_res$W
  S = clust_res$S
  alpha = clust_res$alpha
  cluster = clust_res$cluster
  #neighbors = clust_res$neighbors
  res_J = clust_res$J
  lambda1 = clust_res$lambda1
  lambda2 = clust_res$lambda2
  J_DR = clust_res$J_DR
  J_HE = clust_res$J_HE
  J_LE = clust_res$J_LE
  

  res = list(cluster=cluster, guide_cluster = clust, lambda1 = lambda1, lambda2 = lambda2,
             droprate=droprate, J_DR = J_DR, J_HE = J_HE, J_LE = J_LE, J = res_J, 
             W = W, V = V, H = H, S_set = local_sim, S = S, weight = alpha)  

#  res = list(cluster=cluster, neighbors=neighbors, 
#             imp_X=imp_X, droprate=droprate, local_sim = local_sim, 
#             J = res_J, W = W, V = V, H = H)
  cat("# MGGE iteration complete!\n")
  #saveRDS(res, paste0(output, "MGGE_res.rds"))
  return(res)
}


clustering_update <- function(lg_X, K, npc, lambda1, lambda2, W=NULL, 
                              V=NULL, H=NULL, D=NULL, local_sim=NULL, iteration=50, 
                              thre=1e-4, out_file="result/clust_res/res.rds", 
                              res_save=TRUE){
  cat("## clustering vars updating...\n")
  
  # initialization
  J_set = NULL
  J_DR = NULL
  J_HE = NULL
  J_LE = NULL
  par1 = NULL
  par2 = NULL
  iter = 1
  N = dim(W)[1]
  L = diag(colSums(S)) - S
  Xe = sqrt(eigen(lg_X %*% t(lg_X))$values)[1]
  Le = eigen(L)$values[1]
  S = Reduce('+', local_sim) / length(local_sim)
  alpha = rep(1/M, M)
  rm(L)
  
  
  for(iter in 1:iteration){
    cat(paste0("### Updating ", iter, "th round of ", iteration, "\n"))
    cat(paste0(strrep("#", round(30*iter/iteration)), " ", 
               100*signif(iter/iteration, digits=2), "%\n"))
    cat("### updating W...\n")
    #
    LH = H %*% t(H)
    sec_ord_g = lambda1 * LH %*% W + D^2 %*% W %*% V %*% t(V)
    g_W = -D^2 %*% lg_X %*% t(V) + sec_ord_g#- 0.36s
    
    ## the N*k * N*k Hessian matrix is corresponse to Wij which expand by col
    # Hessian = kronecker(Matrix(diag(1, npc), sparse=TRUE), Matrix(lambda2 * LH, sparse=TRUE)) + 
    #   kronecker(V %*% t(V), Matrix(D^2, sparse=TRUE))  #- sparse 0.66s
    #      Hessian = Matrix(Hessian, sparse=TRUE)   #- 10s
    # sec_ord_g = Hessian %*% matrix(W, ncol=1)   
    
    W = W * ((D^2 %*% lg_X %*% t(V)) / sec_ord_g)
    
    
    cat("### updating V...\n")
    V = V * ((t(W) %*% D^2 %*% lg_X) / (t(W) %*% D^2 %*% W %*% V))
    
    #* 5\ can not normalize W & V
    #norms = NormalizeUV(W, V)
    #W = norms$U
    #V = norms$V
    #rm(norms)
    
    
    cat("### updating H...\n")
    #d = -2 * W %*% t(W)
    d = colSums(W^2)
    d = matrix(d, nrow=N, ncol=N)
    d = d + t(d) - 2*W %*% t(W)
    
    #- H = eigen(0.5*lambda1 * d + lambda2 * S)$vectors[, (N-K+1):N]
    D_s = Matrix(diag(colSums(S)), sparse = TRUE)
    H = H * ((lambda2 * S %*% H) / ((lambda1 / 2 * d + lambda2 * D_s) %*% H))
    
    norms = rowSums(H)    #* 5\row based H normalize
    norms[norms == 0] = 1e-10
    H = H / matrix(norms, N, K)
    
    
    cat('### updating S...\n')
    d_H = matrix(rowSums(H^2), N, N) + t(matrix(rowSums(H^2), N, N)) - 2* H %*% t(H)
    A = lapply(1:length(local_sim), function(i){
      local_sim[[i]] * alpha[i]
      })
    A = Reduce('+', A)
    S = t(sapply(1:N, function(i){
      id = which(A[i, ] > 0)
      ad = (A[i, id] - 0.25*lambda2*d_H[i, id]) / sum(alpha)
      res = proj_sim(ad)
      return(res)
    }))
    

    cat('### updating alpha...\n')
    alpha = sapply(1:M, function(i){
      alpha[i] = 0.5 / norm(S - local_sim[[i]], 'F')
    })
    
    
    #** 6\update parameters
    if(iter %% 10 == 0){
      #browser()
      He = sqrt(eigen(H %*% t(H))$values)[1]
      Ve = sqrt(eigen(V %*% t(V))$values)[1]
      lambda1 = Xe / (He + Ve)
      
      de = eigen(d)$values[1]
      lambda2 = (lambda1 * eigen(d)$values[1]) / (2 * Le)
      par1[(iter-9):iter] = lambda1
      par2[(iter-9):iter] = lambda2
    }
    
    
    # cost calculation
    # if verbose
    J_DR[iter] = norm(D %*% (lg_X - W %*% V), 'F')^2
    J_HE[iter] = sum(diag(t(W) %*% LH %*% W))
    J_LE[iter] = sum(diag(t(H) %*% S %*% H))
    J_set[iter] = J_DR[iter] + lambda1 * J_HE[iter] + lambda2 * J_LE[iter]
    cat("### Current cost:", J_set[iter], "\n")
    
    
    if(iter > 1){
      var_J = abs(J_set[iter] - J_set[iter-1])
      
      # convergence check
      if(var_J <= thre){
        cluster = sapply(1:N, function(i){
          which.max(H[i, ])
        })
        
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
  }
  
  #cluster = kmeans(H, K)$cluster
  cluster = sapply(1:N, function(i){
    which.max(H[i, ])
  })

  
  if(FALSE){
    dist = colSums(W^2)
    dist = matrix(dist, nrow=N, ncol=N)
    dist = dist + t(dist) - 2*W %*% t(W)
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
  }
  

  
  res = list(W = W, V = V, H = H, cluster = cluster,
             J = J_set, J_DR = J_DR, J_HE = J_HE, J_LE = J_LE,
             lambda1 = par1, lambda2 = par2, alpha = alpha, S = S)   #* lambda save
  if(res_save){
    saveRDS(res, out_file)
  }
  return(res)
  
}


imputing_update <- function(lg_X, cluster, neighbors, drop_thre=0.5, out_file="result/imp_res/res.rds", ncores=1, imp_iteration=100, res_save=TRUE){
  # return log-tran X
  # lg 1.01
  # not parallel version on imputation
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

