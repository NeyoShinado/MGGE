library(rsvd)

npc_cal <- function(X, K, var_thre=0.4){
  # input var_thre for small dataset with cell number < 1000
  # or will be set 0.8 by default for big dataset with noise
  
  if(is.data.frame(X)){
    X = matrix(X)
  }else if(!is.matrix(X)){
    stop("Type of argument X must be dataframe of matrix!")
  }
  
  if(length(nrow(X)) <= 1000){
    pca = prcomp(X)
    eigs = (pca$sdev)^2
    var_cum = cumsum(eigs)/sum(eigs)
    if(max(var_cum) <= var_thre){
      npc = length(var_cum)
    }else{
      npc = which.max(var_cum > var_thre)
      npc = max(npc, K)
    }
  }else{
    var_thre = 0.8
    pca = rpca(X, k=1000, center=TRUE, scale=FALSE)
    eigs = (pca$sdev)^2
    var_cum = cumsum(eigs)/sum(eigs)
    if(max(var_cum) <= var_thre){
      npc = length(var_cum)
    }else{
      npc = which.max(var_cum > var_thre)
      npc = max(npc, K)
    }
  }
  return(npc)  
}