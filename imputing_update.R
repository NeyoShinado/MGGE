library(parallel)
library(kernlab)
library(doParallel)
library(penalized)

#* modified for rare cluster
# filter genes with low potentiality being a dropout gene
find_va_genes = function(parslist, subcount){
  point = log10(1.01)
  valid_genes = which((colSums(subcount) > point * nrow(subcount)) &
                         complete.cases(t(parslist)))
  if(length(valid_genes) == 0) return(valid_genes)
  # find out genes that violate assumption
  mu = parslist["mu", ]
  sgene1 = which(mu <= log10(1+1.01))    #% exclude less mean value
  # sgene2 = which(mu <= log2(10+1.01) & mu - parslist[,5] > log2(1.01))
  
  dcheck1 = dgamma(mu+1, shape = parslist["alpha", ], rate = parslist["beta", ])
  dcheck2 = dnorm(mu+1, mean = parslist["mu", ], sd = parslist["sigma", ])
  sgene3 = which(dcheck1 >= dcheck2 & mu <= 1)    #% exclude large dropout rate but less mean value
  sgene = union(sgene1, sgene3)
  valid_genes = setdiff(valid_genes, sgene)
  return(valid_genes)
}


impute_nnls = function(N_gene, cellid, subcount, id_nei, geneid_drop, 
                       geneid_obs, nbs, distc, imp_iteration){
  yobs = subcount[id_nei, ]
  if (length(geneid_drop) == 0 | length(geneid_drop) == N_gene) {
    yimpute = yobs
    local_sim = rep(0, length(nbs))
    return(list(yimpute=yimpute, cell_local_sim=local_sim)) 
  }  
  yimpute = rep(0, N_gene)
  local_sim = rep(0, length(nbs))
  
  xx = matrix(subcount[-id_nei, geneid_obs], nrow=length(geneid_obs))   # training cells
  yy = matrix(subcount[id_nei, geneid_obs])   # predict cell
  ximpute = matrix(subcount[-id_nei, geneid_drop], nrow=N_gene - length(geneid_obs))
  num_thre = 500   #* why select cells
  if(ncol(xx) >= min(num_thre, nrow(xx))){
    if (num_thre >= ncol(xx)){
      new_thre = round((2*ncol(xx)/3))
    }else{ new_thre = num_thre}
    filterid = order(distc[-1])[1: new_thre]
    xx = xx[filterid, , drop = FALSE]
    ximpute = ximpute[, filterid, drop = FALSE]
  }
  set.seed(cellid)
  # cat("### training nnls model...\n")
  nnls = penalized(yy, penalized = xx, unpenalized = ~0,
                   positive = TRUE, lambda1 = 0, lambda2 = 0, 
                   maxiter = imp_iteration, trace = FALSE)   #* 3000 iteration maybe too large
  if(length(geneid_drop) > 1){
    ynew = penalized::predict(nnls, penalized = ximpute, unpenalized = ~0)[1]
  }else{
    ynew = penalized::predict(nnls, penalized = ximpute, unpenalized = ~0)[,1]
  }
  yimpute[geneid_drop] = t(ynew)
  yimpute[geneid_obs] = yobs[geneid_obs]
  maxobs = apply(subcount, 2, max)
  yimpute[yimpute > maxobs] = maxobs[yimpute > maxobs]
  coeff = coefficients(nnls, "all")   
  local_sim[which(nbs == 1)] = coeff
  #local_sim[as.logical(nbs)] = coeff
  return(list(yimpute=yimpute, cell_local_sim=local_sim))
}


# imputation with cluster & low-dim fea dist
# input clustering cluster vec & low-dim fea neighbors list of logical neighbor matrix and dist matrix
# output imputed lg_data, exp-based droprate & regressive coeff(local sim)
imputation_wlabel_model = function(count, cluster, neighbors, point, drop_thre, ncores=1, imp_iteration=100){
  if(!(class(cluster) %in% c("character", "numeric", "integer"))){
    stop("cell_labels should be a character or integer vector!")
  }
  
  count = as.matrix(count)
  N = nrow(count)
  P = ncol(count)
  count_imp = count
  droprate = matrix(0, N, P)
  local_sim = matrix(0, N, N)
  
  dist_list = neighbors$dist_list  
  neighbors = neighbors$neighbors
  
  
  # fitting mixture model & dropout indentification
  nclust = sum(!is.na(unique(cluster)))
  cl = makeCluster(ncores, outfile="")
  #registerDoParallel(cl)
  
  # N * P matrix of droprate
  #clust_droprate = lapply(1:nclust, function(cc){
  for(cc in 1:nclust){
    cat(paste("### estimating dropout probability for type", cc, "...\n"))
    # fit_model by_gene
    cells = which(cluster == cc)
    pars_bycell = get_mix_parameters(count = count[cells, , drop = FALSE], 
                       point = log10(1.01), ncores = ncores)
    if(length(cells) <= 1){ next }   #* rare cluster
    cat("### searching for valid genes ...\n")
    valid_genes = find_va_genes(pars_bycell, subcount = count[cells, ])  #* diff clust may have diff va_genes
    if(length(valid_genes) <= 10){ next }   # not obvious dropout gene
    
    subcount = count[cells, valid_genes, drop = FALSE]
    N_gene = length(valid_genes)
    N_cell = nrow(subcount)
    pars_bycell = pars_bycell[, valid_genes, drop = FALSE]
    
    # dropout probability
    clust_droprate = sapply(1:N_gene, function(i) {
      wt = calculate_weight(subcount[, i], pars_bycell[, i])
      return(wt[, 1])
    })
    mucheck = sweep(subcount, MARGIN = 2, pars_bycell["mu", ], FUN = ">")
    clust_droprate[mucheck & clust_droprate > drop_thre] = 0  # N_cc * va_gene vector
    
    # finding dropout genes
    setA = lapply(1:N_cell, function(cellid){
      which(clust_droprate[cellid, ] > drop_thre)
    })
    setB = lapply(1:N_cell, function(cellid){
      which(clust_droprate[cellid, ] <= drop_thre)
    })
    gc()
    
    # imputation based on neighborss
    cat(paste("### imputing dropout values for cluster", cc, "...\n"))
    cellid = NULL
    # parallel model
    if(FALSE){
      subres = foreach(cellid = 1:N_cell, .packages = c("penalized"), 
                       .combine = rbind, .export = c("impute_nnls")) %dopar% {
                         ##sink(paste0(out_dir, "log.txt"), append=TRUE))
                         ##print(paste("imputing dropout values for type", cc, "\n")
                         if (cellid %% 10 == 0) {gc()}
                         if (cellid %% 100 == 0) {cat(paste0("### imputing ", cellid, "th cell of ", N_cell, "\n"))}
                         nbs = neighbors[cells[cellid], ]   # to be imputed cellid of nei_cell
                         cell_nei = which(which(nbs==1) == cells[cellid])
                         if (length(which(nbs==1)) == 0) {
                           cat("### there's no neighbors of", cellid, "th cell.\n")
                           return(NULL)
                         }
                         imp_subcount = count[which(nbs==1), valid_genes, drop=FALSE]   # Nnei * N_gene
                         geneid_drop = setA[[cellid]]
                         geneid_obs = setB[[cellid]]
                         imputation = try(impute_nnls(N_gene, cellid = cellid, imp_subcount, cell_nei, geneid_drop, 
                                                      geneid_obs, nbs, distc = dist_list[[cells[cellid]]], imp_iteration=imp_iteration), silent = TRUE)
                         if (class(imputation) == "try-error") {
                           # print(y)
                           y = subcount[cellid, , drop = FALSE]
                           cell_local_sim = rep(0, N)
                         }else{
                           y = imputation$yimpute
                           cell_local_sim = imputation$cell_local_sim
                         }
                         
                         res = list(clust_y=y, clust_local_sim=cell_local_sim)
                         return(res)
                       }
      
    }
    subres = lapply(1:N_cell, function(cellid) {
         ##sink(paste0(out_dir, "log.txt"), append=TRUE))
         ##print(paste("imputing dropout values for type", cc, "\n")
         if (cellid %% 10 == 0) {gc()}
         if (cellid %% 100 == 0) {cat(paste0("### imputing ", cellid, "th cell of ", N_cell, "\n"))}
         nbs = neighbors[cells[cellid], ]   # to be imputed cellid of nei_cell
         cell_nei = which(which(nbs==1) == cells[cellid])
         if (length(which(nbs==1)) == 0) {
           cat("### there's no neighbors of", cellid, "th cell.\n")
           return(NULL)
         }
         imp_subcount = count[which(nbs==1), valid_genes, drop=FALSE]   # Nnei * N_gene
         geneid_drop = setA[[cellid]]
         geneid_obs = setB[[cellid]]
         imputation = try(impute_nnls(N_gene, cellid = cellid, imp_subcount, cell_nei, geneid_drop, 
                                      geneid_obs, nbs, distc = dist_list[[cells[cellid]]]), silent = TRUE)
         if (class(imputation) == "try-error") {
           # print(y)
           y = subcount[cellid, , drop = FALSE]
           cell_local_sim = rep(0, N)
         }else{
           y = imputation$yimpute
           cell_local_sim = imputation$cell_local_sim
         }
         
         res = list(clust_y=y, clust_local_sim=cell_local_sim)
         return(res)
       })
    
    
    clust_y = t(sapply(subres, function(i){return(i[[1]])}))
    clust_local_sim = t(sapply(subres, function(i){return(i[[2]])}))
    
    # update clust_output info
    local_sim[cells, ] = clust_local_sim
    count_imp[cells, valid_genes] = clust_y
    droprate[cells, valid_genes] = clust_droprate
  }
  
  #stopCluster(cl)
  count_imp[count_imp < point] = point
  # local-sim uniation
  local_sim = t(apply(local_sim, 1, function(i){
    if(max(i) != min(i)){
      #i = (i - min(i)) / (max(i) - min(i))
      i = i / sum(i)
    }
    return(i)
  }))
  local_sim = (local_sim + t(local_sim)) / 2
  
  return(list(count_imp = count_imp, droprate=droprate, local_sim=local_sim))

}
