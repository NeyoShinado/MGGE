setwd("E:/Project/Paper_debug/Clustering algorithm/MGGE/")
source("data_normalize.R")
source("npc_cal.R")
source("localsim_integrated_update.R")
source("imputing_update.R")
#source("var_update.R")
source("informative_gene_selection.R")
source("NormalizeUV.R")
source("constructW.R")

library(Matrix)
library(parallel)
library(aricode)
library(Seurat)
library(R.matlab)


# dataset info
datas = c("Biase/Biase.rds", "Deng_GSE45719/mouse_embryo.rds", "Zeisel/Zeisel.rds", 
          "mouse1/mouse1.rds", "mouse2/mouse2.rds", "human3/human3.rds", "Hrvatin/Hrvatin.rds")
data_path = "E:/Project/dataset/Bioinformatics/scRNA/Selected data/"

for(i in c(3:6)){
  try({
  s = Sys.time()
  
  X = readRDS(paste0(data_path, datas[i]))
  gt = X$gt
  X = X$expr
  N = nrow(X)
  P = ncol(X)
  message(paste0("## Loading raw data of ",N , "*", P, " and labels...\n"))
  
  
  ## filter zero gene & cell and remain high_var gene
  ## as well as normalization & log-transform without pseudo count
  message("## Data nomalization and log-transformation...\n")
  lg_X = data_normalize(X, N, P, gt, mu_probs=0.2, cv_probs=0.2)
  gt = lg_X$gt
  lg_X = as.matrix(lg_X$count_hv)
  
  #lg_X = informative_gene_selection(lg_X, delta=0.7)
  
  data_object = CreateSeuratObject(counts = Matrix(t(lg_X), sparse=TRUE),
                            project = "MGGE", min.cells = 3)
  ## Normalizing
  data_object = NormalizeData(data_object, normalization.method = "LogNormalize", scale.factor = 10000)
  
  ## Highly variable feature selection
  #** 2\select about thousand genes
  data_object = FindVariableFeatures(data_object, selection.method = "vst", nfeatures = 3000)
  id = VariableFeatures(data_object)
  if(all(id %in% colnames(lg_X))){
  }else{
    id = id[-which(!(id %in% colnames(lg_X)))]
  }
  
  message("## Finally choose ", length(id), " genes for MGGE...")

  lg_X = lg_X[, id]
  rm(X, data_object)

  
  ##* cluster number estimation
  #K = cluster_number_estimation(lg_X)
  K = length(unique(gt))
  message("## Estimated cluster number:", K, "\n")
  
  
  ## dimension reduction(npc) calculation
  npc = npc_cal(as.matrix(lg_X), K, var_thre=0.8)
  message("## Caluculated dim reduction npc:", npc, "\n")
  
  
  ## cluster update
  # clustering update, 200 iteration by default
  #* stop update when Jbefore < Jafter
  files = list.files("./scImpute/", pattern="*.R")
  sapply(paste0("./scImpute/", files), source)
  
  
  #** 3\test of gt local_sim
  if(FALSE){
    mask = unique(gt)
    F = sapply(1:K, function(i){
      F = rep(0, N)
      id = which(gt == mask[i])
      F[id] = 1
      return(F)
    })
    S = F %*% t(F)    
  }
  
  output = strsplit(datas[i], split="/")[[1]][1]
  S = readMat(paste0("./dataset/", output, "_W.mat"))$W
  
  
  res <- var_update(lg_X, K, npc, S=S, lambda1=2, lambda2=2, 
                    iteration=1, clust_iteration=300, imp_iteration=3000, 
                    res_save=FALSE)
  clust = res$cluster
  nmi = NMI(clust, gt)
  ari = ARI(clust, gt)
  res$nmi = nmi
  res$ari = ari
  message("## NMI: ", nmi, "   ARI: ", ari, "\n")
  
  
  time = Sys.time() - s
  message("## Consume", time, "seconds.\n")
  
  # save res
  output = strsplit(datas[i], split="/")[[1]][1]
  
  saveRDS(res, paste0("result/test/", output, 'normal_W.rds'))
  })
}
