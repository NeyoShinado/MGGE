# clustering vars updating
library(Matrix)
library(bignmf)


var_update <- function(lg_X, K, npc, S, neighbors, guide_cluster, gt, sigmac=3, sigmag=2, lambda1=0.1, lambda2=0.1, 
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

  
  imp_iter = 1
  droprate = list()
  imp_X = 0
  
  
### imputing_update
  cat(paste0("## carrying ", imp_iter, "th cluster imputation...\n\n"))
  #*clust = clust_set[[i]]
  imp_res = imputing_update(log10(10^(lg_X) + 0.01), guide_cluster, neighbors, 
                            imp_iteration=imp_iteration, 
                            out_file=paste0(output, "imp_res/", imp_iter, "th_impres.rds"), 
                            res_save=res_save)
  
  #imp_X = log10(10^(imp_X) - 0.01)   # log-tran 1.01 for imputing, 1.0 for clustering
  
  ## arrange imp_output
  #local_sim = imp_res$local_sim
  imp_X = imp_res$imp_X
  droprate$data = imp_res$droprate
  mean_drop = apply(droprate, 2, mean)
  # L\M\H means high\middle\low droprate
  # 0.6\0.9\1 percent of drop
  Lid = which(meandrop < 0.4)
  Mid = which(meandrop < 0.65 & meandrop >= 0.4)
  Hid = which(meandrop >= 0.65)
  droprate$Lid = Lid
  droprate$Mid = MLid
  droprate$Hid = Hid
  
  
  # divide gene
  split_data = list()
  split_data[[1]] = imp_X[, Lid]
  split_data[[2]] = imp_X[, Mid]
  split_data[[3]] = droprate[, Hid]
  
  
  #* 12\weighted cell
  Dc = rowSums(droprate)
  Dg = colSums(droprate)[Mid]
  Dc = Matrix(diag(exp(-sigmac * Dc / P)), sparse=TRUE)
  Dg = Matrix(diag(exp(-sigmag * Dg / N)), sparse=TRUE)    # scale factor 2 for small dataset and 3 for the big one
  
  
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
  
  
  ### Clustering update
  cat(paste0("## running ", imp_iter, "th vars update via clustering and imputation...\n"))

  
  # init vars
  #H = matrix(runif(dim(lg_X)[1]*K), N, K)
  H = eigen(diag(colSums(clust_set[[1]])) - clust_set[[1]])$vector
  H = H[, (N-K):(N-1)]
  H[H<=0] = 1e-10
  W = matrix(runif(N * npc), N, npc) 
  V = lapply(1:3, function(i){
    pi = dim(split_data[[i]])[2]
    Vi = matrix(runif(npc * pi), npc, pi)
    return(Vi)
  })
  
  
  
  #* 7\ using split imp_X & dropout
  clust_res = clustering_update(split_data, K, npc, lambda1, lambda2, W=W, V=V, droprate=droprate,
                                H=H, Dc=Dc, Dg = Dg, local_sim=local_sim, iteration=clust_iteration,
                                out_file=paste0(output, "clust_res/localsim_integrated_clustres.rds"), 
                                res_save=res_save)
  #clust_res = clustering_update(imp_X, K, npc, lambda1, lambda2, W=W, V=V, H=H, Dc=Dc, L=L, iteration=clust_iteration,
  #                             out_file=paste0(output, "clust_res/localsim_integrated_clustres.rds"), res_save=res_save)
  
  H = clust_res$H
  V = clust_res$V
  W = clust_res$W
  dw = clust_res$dw
  alpha = clust_res$alpha
  cluster = clust_res$cluster
  nmi = clust_res$nmi
  res_J = clust_res$J
  lambda1 = clust_res$lambda1
  lambda2 = clust_res$lambda2
  J_DR = clust_res$J_DR
  J_HE = clust_res$J_HE
  J_LE = clust_res$J_LE
  #neighbors = clust_res$neighbors
  
  
  res = list(cluster=cluster, guide_cluster = guide_cluster, lambda1 = lambda1, lambda2 = lambda2,
             droprate=droprate, J_DR = J_DR, J_HE = J_HE, J_LE = J_LE, J = res_J, nmi = nmi,
             W = W, V = V, H = H, Dc = Dc, Dg = Dg, S_set = local_sim, dw = dw, weight = alpha, imp_X = imp_X)  
  
  #  res = list(cluster=cluster, neighbors=neighbors, 
  #             imp_X=imp_X, droprate=droprate, local_sim = local_sim, 
  #             J = res_J, W = W, V = V, H = H)
  cat("# MGGE iteration complete!\n")
  #saveRDS(res, paste0(output, "MGGE_res.rds"))
  return(res)
}


clustering_update <- function(lg_X, K, npc, lambda1, lambda2, W=NULL, V=NULL, droprate=NULL,
                              H=NULL, Dc=NULL, Dg=NULL, S=NULL, iteration=50, 
                              thre=1e-4, out_file="result/clust_res/res.rds", 
                              res_save=TRUE){
  cat("## clustering vars updating...\n")

  # initialization
  M = length(lg_X)
  N = dim(W)[1]
  J_set = NULL
  J_DR = matrix(0, iteration, M)
  J_HE = NULL
  J_LE = NULL
  par1 = NULL
  par2 = NULL
  nmi = matrix(0, 1, 3)
  colnames(nmi) = c("nmiH", "nmiW_L", "nmiW_sp")
  Xe = c(sqrt(eigen(lg_X[[1]] %*% t(lg_X[[1]]))$values)[1], sqrt(eigen(lg_X[[2]] %*% t(lg_X[[2]]))$values)[1], 
         sqrt(eigen(lg_X[[3]] %*% t(lg_X[[3]]))$values)[1])
  Se = eigen(S)$values[1]
  Ds = diag(colSums(S))
  L = Ds - S
  alpha = matrix(1/3, iteration, 3)
  Lid = droprate$Lid
  Mid = droprate$Mid
  Hid = droprate$Hid
  drop_H = droprate[, Hid]
  rm(droprate)
  
  
  for(iter in 1:iteration){
    cat(paste0("### Updating ", iter, "th round of ", iteration, "\n"))
    cat(paste0(strrep("#", round(30*iter/iteration)), " ", 
               100*signif(iter/iteration, digits=2), "%\n"))
    cat("### updating W...\n")
    #
    LH = H %*% t(H)
    pl = (alpha[iter, 1] + alpha[iter, 3]) * Matrix(diag(1, N), sparse=TRUE) + alpha[iter, 2] * Dc^2
    pr = alpha[iter, 1] * V[[1]] %*% t(V[[1]]) + alpha[iter, 2] * V[[2]] %*% Dg^2 %*% t(V[[2]]) + 
      alpha[iter, 3] * V[[3]] %*% t(V[[3]])
    sec_ord_g = lambda1 * LH %*% W + pl %*% W %*% pr
    nW = alpha[iter, 1] * lg_X[[1]] %*% t(V[[1]]) + alpha[iter, 2] * Dc^2 %*% lg_X[[2]] %*% Dg^2 %*% t(V[[2]]) + 
            alpha[iter, 3] * lg_X[[3]] %*% t(V[[3]])
    W = W * ( nW/ sec_ord_g)
    rm(pr, pl, nW)
    #g_W = -Dc^2 %*% lg_X %*% t(V) + sec_ord_g#- 0.36s
    ## the N*k * N*k Hessian matrix is corresponse to Wij which expand by col
    # Hessian = kronecker(Matrix(diag(1, npc), sparse=TRUE), Matrix(lambda2 * LH, sparse=TRUE)) + 
    #   kronecker(V %*% t(V), Matrix(Dc^2, sparse=TRUE))  #- sparse 0.66s
    #      Hessian = Matrix(Hessian, sparse=TRUE)   #- 10s
    # sec_ord_g = Hessian %*% matrix(W, ncol=1)   
    
    
    cat("### updating V...\n")
    V[[1]] = V[[1]] * ((t(W) %*% lg_X[[1]]) / (t(W) %*% W %*% V[[1]]))
    V[[2]] = V[[2]] * ((t(W) %*% Dc^2 %*% lg_X[[2]] %*% Dg^2) / (t(W) %*% Dc^2 %*% W %*% V[[2]] %*% Dg^2))
    V[[3]] = V[[3]] * ((t(W) %*% lg_X[[3]]) / (t(W) %*% W %*% V[[3]]))
    
    #* 5\ can not normalize W & V
    #norms = NormalizeUV(W, V)
    #W = norms$U
    #V = norms$V
    #rm(norms)
    
    
    cat("### updating H...\n")
    #d = -2 * W %*% t(W)
    d = rowSums(W^2)
    d = matrix(d, nrow=N, ncol=N)
    d = d + t(d) - 2*W %*% t(W)
    
    #- H = eigen(0.5*lambda1 * d + lambda2 * S)$vectors[, (N-K+1):N]
    H = H * ((lambda2 * S %*% H) / ((lambda2 * Ds + lambda1 / 2 * d) %*% H))
    
    norms = rowSums(H)    #* 5\row based H normalize
    norms[norms == 0] = 1e-10
    H = H / matrix(norms, N, K)
    
    
    cat('### updating alpha...\n')
    J_DR[iter, ] = t(sapply(1:M, function(i){
      if(i == 2){
        cost = norm(Dc^2 %*% (lg_X[[i]] - W %*% V[[i]]) %*% Dg^2, 'F')^2
      }else{
        cost = norm(lg_X[[i]] - W %*% V[[i]], 'F')^2
      }
      
      return(cost)
    }))
    
    alpha[iter, ] = t(sapply(1:M, function(i){
      res = 0.5 / sqrt(J_DR[iter, i])
      return(res)
    }))
    alpha[iter, ] = alpha[iter, ] / sum(alpha[iter, ])
    
    # integration update
    if(FALSE){
      cat('### updating S...\n')
      d_H = matrix(rowSums(H^2), N, N) + t(matrix(rowSums(H^2), N, N)) + 2* H %*% t(H)
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
  
    }
    
    
    #** 6\update parameters    option one
    if(iter %% 10 == 0){
      He = sqrt(eigen(H %*% t(H))$values)[1]
      Ve = c(sqrt((eigen(V[[1]] %*% t(V[[1]]))$values)[1]), sqrt((eigen(V[[2]] %*% t(V[[2]]))$values)[1]),
                  sqrt((eigen(V[[3]] %*% t(V[[3]]))$values)[1]))
      lambda1 = (alpha[iter, ] * Xe) / (He + alpha[iter, ] * Ve)   
      
      de = eigen(d)$values[1]
      lambda2 = (lambda1 * de) / (2 * Se)
      #par1[iter] = lambda1
      #par2[iter] = lambda2
      par1[(iter-9):iter] = lambda1
      par2[(iter-9):iter] = lambda2
    }
    
    
    # cost calculation
    # if verbose
    J_LE[iter] = sum(diag(t(H) %*% L %*% H))
    J_HE[iter] = sum(diag(t(W) %*% LH %*% W))
    J_set[iter] = alpha[iter, ] * J_DR[iter, ] + lambda1 * J_HE[iter] + lambda2 * J_LE[iter]
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
    
    
    # recording nmi
    #cluster = specc(as.matrix(d), K)@.Data
    #nmiW_sp = NMI(cluster, gt)
    nmiW_sp = 0
    
    F = eigen(d)$vectors
    F = F[, (N-K+1):(N)]
    clust = kmeans(F, K)$cluster
    nmiW_L = NMI(clust, gt)
    
    cluster = sapply(1:N, function(i){
      which.max(H[i, ])
    })
    nmiH = NMI(cluster, gt)
    
    nmi = rbind(nmi, c(nmiH, nmiW_L, nmiW_sp))
    
  }
  
  
  res = list(W = W, V = V, H = H, cluster = cluster, Dc = diag(Dc), Dg = diag(Dg),
             J = J_set, J_DR = J_DR, J_HE = J_HE, J_LE = J_LE, nmi = nmi,
             lambda1 = par1, lambda2 = par2, alpha = alpha, dw = d)   #* lambda save
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

