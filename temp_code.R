signif(norm(D %*% (lg_X - W %*% V), 'F')^2 , digits = 2)
signif(lambda1 * sum(diag(t(W) %*% LH %*% W)), digits = 3)
signif(lambda2 * sum(diag(t(H) %*% L %*% H)), digits=3)


nmi = sapply(0:20, function(i){
  clust_res = readRDS(paste0("result/clust_res/", i, "th_clustres.rds"))
  NMI(clust_res$cluster, gt)
})

nzero = sapply(1:N, function(i){
  length(which(clust_res$H[i, ] == 0))
})
hist(nzero)


### loading dataset
setwd("E:/Project/Paper_debug/Clustering algorithm/MGGE/")
data = readRDS("./dataset/mouse_embryo.rds")
gt = data$label[which(data$singlecell)]
data = data$reads[which(data$singlecell), ]

library(SingleCellExperiment)
library(SC3)
library(scater)

clust_res = sapply(1:9, function(i){
  clust = readRDS(paste0("result/clust_res/", i, "th_clustres.rds"))$J_set
  return(J_set)
  #droprate = rowSums(res$droprate)
  #D = (P - droprate) / P   #* dropout weight
})
J_set = as.vector(J_set)
plot(J_set[2:450])


clust = list()
for(lambda1 in 10^seq(-3, 5, 2)){
  for(lambda2 in 10^seq(-4, -2)){
    out_dir = paste0("result/lambda_", lambda1, "_", lambda2)
    clust[[paste0(lambda1, "_", lambda2)]] = sapply(1:10, function(i){
      res = readRDS(paste0(out_dir, "/clust_res/", i, "th_clustres.rds"))
      cluster = res$cluster
      nmi = NMI(cluster, gt)
      return(nmi=nmi)
    })
  }
}


J = list()
for(lambda1 in 10^seq(-3, 5, 2)){
  for(lambda2 in 10^seq(-4, -2)){
    out_dir = paste0("result/lambda_", lambda1, "_", lambda2)
    J[[paste0(lambda1, "_", lambda2)]] = as.vector(sapply(1:10, function(i){
      res = readRDS(paste0(out_dir, "/clust_res/", i, "th_clustres.rds"))
      J = res$J_set
      return(J)
    }))
  }
}


nmi = sapply(1:5, function(i){
  cluster = kmeans(as.matrix(lg_X), 12, iter.max = 1000)$cluster
  nmi = NMI(cluster, gt)
  return(nmi)
})


files = list.files("result/", pattern="lambda_*")
nmi = lapply(files, function(i){
  try({cluster <- readRDS(paste0("result/", i, "/MGGE_res.rds"))$cluster
  nmi<- NMI(cluster, gt)
  return(nmi)})
})