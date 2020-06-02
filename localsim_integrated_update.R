# clustering vars updating
library(Matrix)
library(bignmf)


var_update <- function(lg_X, K, npc, clust_set, guide_cluster, sigmac=3, sigmag=2, lambda1=0.1, lambda2=0.1, 
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
  NN = 5
  #NN = as.integer(log10(N))         #* 10 at most
  #NN = min(50, max(10, 10 * NN))
  res = constructW(lg_X, NN)
  W = res$W
  neighbors = res$neighbors
  

  #clust = rep(1, N)
  
  #clust = sapply(1:N, function(i){    
  #  return(order(H[i, ], decreasing = TRUE)[1])
  #})
  
  
  #H = matrix(runif(dim(lg_X)[1]*K), N, K)
  
  imp_iter = 1
  #* 13\ ensemble imp_X
  local_sim = list()
  M = length(clust_set)
  droprate = 0
  imp_X = 0


  for(i in c(1:M)){
    #- par-iteration is unnessary for MGGE any more
    cat(paste0("## carrying ", i, "th cluster imputation...\n\n"))
    ###clust = clust_set[[i]]
    
    imp_res = imputing_update(log10(10^(lg_X) + 0.01), guide_cluster, neighbors, 
                              imp_iteration=imp_iteration, 
                              out_file=paste0(output, "imp_res/", imp_iter, "th_impres.rds"), 
                              res_save=res_save)
    
    #imp_X = log(exp(imp_X) - 0.01)   # log-tran 1.01 for imputing, 1.0 for clustering
    #imp_X = Matrix(imp_X, sparse=TRUE)    
    
    # accumulate update with normalized local_sim
    #local_sim = imp_res$local_sim
    

    imp_X = imp_X + imp_res$imp_X
    droprate = droprate + imp_res$droprate
    
    #** 3\test of best local_sim / imp_local_sim
    if(FALSE){
      local_sim[[i]] = imp_res$local_sim
      local_sim[[i]] = t(sapply(1:N, function(j){
        res = rep(0, N)
        id = order(local_sim[[i]][j, ], decreasing=TRUE)[1:NN]
        res[id] = local_sim[[i]][j, ][id]
        return(res)
      }))
      local_sim[[i]]  = (local_sim[[i]] + t(local_sim[[i]])) / 2
    }
    ###local_sim[[i]] = construct_locsim(imp_res$imp_X, droprate, NN, K)
    local_sim = clust_set
    
  }
  imp_X = imp_X / M
  
  #* 12\weighted cell
  droprate = droprate / M
  Dc = rowSums(droprate)
  Dg = colSums(droprate)
  #* scale factor 2 for small dataset and 3 for the big one
  Dc = Matrix(diag(exp(-sigmac * Dc / P)), sparse=TRUE)
  Dg = Matrix(diag(exp(-sigmag * Dg / N)), sparse=TRUE)
  
  #* no dropout weight
  #Dc = Matrix(diag(1, N), sparse=TRUE)
  #Dg = Matrix(diag(1, P), sparse=TRUE)
  
  
  ### Clustering update
  cat(paste0("## running ", imp_iter - 1, "th vars update via clustering and imputation...\n"))
  
  
  #imp_X = informative_gene_selection(imp_X, delta=0.7)   #- select informative genes for clustering
  
  map = unique(cluster)
  H = sapply(1:K, function(i){
    res = rep(0, N)
    id = which(cluster == map[i])
    res[id] = 1
    return(res)
  })
  W = matrix(runif(N * npc), N, npc) 
  V = matrix(runif(P * npc), npc, P) 
  

  #* 7\ using imp_X
  clust_res = clustering_update(imp_X, K, npc, lambda1, lambda2, W=W, V=V, 
                                H=H, Dc=Dc, Dg = Dg, local_sim=local_sim, iteration=clust_iteration,
                                out_file=paste0(output, "clust_res/localsim_integrated_clustres.rds"), 
                                res_save=res_save)
  #clust_res = clustering_update(imp_X, K, npc, lambda1, lambda2, W=W, V=V, H=H, Dc=Dc, L=L, iteration=clust_iteration,
  #                             out_file=paste0(output, "clust_res/localsim_integrated_clustres.rds"), res_save=res_save)
  
  H = clust_res$H
  V = clust_res$V
  W = clust_res$W
  S = clust_res$S
  Dc = clust_res$Dc
  alpha = clust_res$alpha
  cluster = clust_res$cluster
  #neighbors = clust_res$neighbors
  res_J = clust_res$J
  lambda1 = clust_res$lambda1
  lambda2 = clust_res$lambda2
  J_DR = clust_res$J_DR
  J_HE = clust_res$J_HE
  J_LE = clust_res$J_LE
  
  
  res = list(cluster=cluster, guide_cluster = guide_cluster, lambda1 = lambda1, lambda2 = lambda2,
             droprate=droprate, J_DR = J_DR, J_HE = J_HE, J_LE = J_LE, J = res_J, 
             W = W, V = V, H = H, Dc = Dc, S_set = local_sim, S = S, weight = alpha, imp_X = imp_X)  
  
  #  res = list(cluster=cluster, neighbors=neighbors, 
  #             imp_X=imp_X, droprate=droprate, local_sim = local_sim, 
  #             J = res_J, W = W, V = V, H = H)
  cat("# MGGE iteration complete!\n")
  #saveRDS(res, paste0(output, "MGGE_res.rds"))
  return(res)
}


clustering_update <- function(lg_X, K, npc, lambda1, lambda2, W=NULL, 
                              V=NULL, H=NULL, Dc=NULL, Dg=NULL, local_sim=NULL, iteration=50, 
                              thre=1e-4, out_file="result/clust_res/res.rds", 
                              res_save=TRUE){
  cat("## clustering vars updating...\n")
  
  # initialization
  M = length(local_sim)
  N = dim(W)[1]
  J_set = NULL
  J_DR = NULL
  J_HE = NULL
  J_LE = matrix(0, iteration, M)
  par1 = NULL
  par2 = NULL
  iter = 1
  Xe = sqrt(eigen(lg_X %*% t(lg_X))$values)[1]
  Ds = diag(colSums(local_sim[[1]]))
  alpha = matrix(1/M, iteration, M)
  
  
  for(iter in 1:iteration){
    cat(paste0("### Updating ", iter, "th round of ", iteration, "\n"))
    cat(paste0(strrep("#", round(30*iter/iteration)), " ", 
               100*signif(iter/iteration, digits=2), "%\n"))
    cat("### updating W...\n")
    #
    LH = H %*% t(H)
    sec_ord_g = lambda1 * LH %*% W + Dc^2 %*% W %*% V %*% Dg^2 %*% t(V)
    #g_W = -Dc^2 %*% lg_X %*% t(V) + sec_ord_g#- 0.36s
    ## the N*k * N*k Hessian matrix is corresponse to Wij which expand by col
    # Hessian = kronecker(Matrix(diag(1, npc), sparse=TRUE), Matrix(lambda2 * LH, sparse=TRUE)) + 
    #   kronecker(V %*% t(V), Matrix(Dc^2, sparse=TRUE))  #- sparse 0.66s
    #      Hessian = Matrix(Hessian, sparse=TRUE)   #- 10s
    # sec_ord_g = Hessian %*% matrix(W, ncol=1)   

    W = W * ((Dc^2 %*% lg_X %*% Dg^2 %*% t(V)) / sec_ord_g)
    
    
    cat("### updating V...\n")
    V = V * ((t(W) %*% Dc^2 %*% lg_X %*% Dg^2) / (t(W) %*% Dc^2 %*% W %*% V %*% Dg^2))
    
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
    S_sum = 0
    for(i in c(1:M)){
      S_sum = S_sum + alpha[i] * local_sim[[i]]
    }
    
  
    #D_s = Matrix(diag(colSums(S)), sparse = TRUE)
    H = H * ((lambda2 * S_sum %*% H) / ((lambda2 * Ds + lambda1 / 2 * d) %*% H))
    
    norms = rowSums(H)    #* 5\row based H normalize
    norms[norms == 0] = 1e-10
    H = H / matrix(norms, N, K)
    
    
    cat('### updating alpha...\n')
    J_LE[iter, ] = t(sapply(1:M, function(i){
      J_LEv = sum(diag(t(H) %*% local_sim[[i]] %*% H))
      return(J_LEv)
    }))
    
    J_HE[iter] = sum(diag(t(W) %*% LH %*% W))
    
    # integration update
    if(FALSE){
      cat('### updating S...\n')
      d_H = matrix(rowSums(H^2), N, N) + t(matrix(rowSums(H^2), N, N)) - 2* H %*% t(H)
      A = lapply(1:length(local_sim), function(i){
        local_sim[[i]] * alpha[i]
      })
      A = Reduce('+', A)
      S = t(sapply(1:N, function(i){
        res = rep(0, N)
        id = which(A[i, ] > 0)
        #* 10\update L with H
        ad = (A[i, id]) / sum(alpha)   #* not H constraint
        #ad = (A[i, id] - 0.25*lambda2*d_H[i, id]) / sum(alpha)
        res[id] = proj_sim(ad)
        return(res)
      }))
      S = (S + t(S)) / 2
      
      alpha[iter, ] = t(sapply(1:M, function(i){
        res = 0.5 / sqrt(J_LE[iter, i])
        return(res)
      }))
      alpha[iter, ] = alpha[iter, ] / sum(alpha[iter, ])
    }
    
    
    #** 6\update parameters
    if(FALSE){  #iter %% 10 == 0){
      He = sqrt(eigen(H %*% t(H))$values)[1]
      Ve = sqrt(eigen(V %*% t(V))$values)[1]
      lambda1 = Xe / (He + Ve)
      
      de = eigen(d)$values[1]
      Se = eigen(S_sum)$values[1]
      lambda2 = (lambda1 * de) / (2 * Se)
      par1[iter] = lambda1
      par2[iter] = lambda2
      #par1[(iter-9):iter] = lambda1
      #par2[(iter-9):iter] = lambda2
    }
    
    lambda2 = 0.5 / sqrt(J_LE[iter, i])
    lambda1 = 0.5 / sqrt(J_HE[iter])
    par1[iter] = lambda1
    par2[iter] = lambda2
    
    # cost calculation
    # if verbose
    J_DR[iter] = norm(Dc^2 %*% (lg_X - W %*% V) %*% Dg^2, 'F')^2
    J_set[iter] = J_DR[iter] + lambda1 * J_HE[iter] + lambda2 * alpha[iter, ] %*% J_LE[iter, ]
    cat("### Current cost:", J_set[iter], "\n")
    
    
    if(iter > 1){
      var_J = abs(J_set[iter] - J_set[iter-1])
      
      # convergence check
      if(var_J <= thre){
        cluster = sapply(1:N, function(i){
          which.max(H[i, ])
        })
        
        res = list(W = W, V = V, H = H, cluster = cluster,
                   J = J_set, J_DR = J_DR, J_HE = J_HE, J_LE = J_LE, 
                   lambda1 = par1, lambda2 = par2, alpha = alpha, S = S)
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
  S = d
  
  res = list(W = W, V = V, H = H, cluster = cluster, Dc = diag(Dc),
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
  imp_res = imputation_wlabel_model(lg_X, cluster, neighbors=neighbors, point=log10(1.01),
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

