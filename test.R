setwd("E:/Project/Paper_debug/Clustering algorithm/MGGE/")
source("data_normalize.R")
source("npc_cal.R")
source("localsim_integrated_update.R")
source("imputing_update.R")
#source("var_update.R")
source("informative_gene_selection.R")
source("constructW.R")


# dataset info
datas = c("Biase/Biase.rds", "Deng_GSE45719/mouse_embryo.rds", "Zeisel/Zeisel.rds", 
          "Baron/mouse1.rds", "Baron/mouse2.rds", "Baron/human3.rds", "Hrvatin/Hrvatin.rds")
data_path = "E:/Project/dataset/Bioinformatics/scRNA/Selected data/"

for(i in 1:6){
#try({
  s = Sys.time()
  
  X = readRDS(paste0(data_path, datas[i]))
  gt = X$gt
  X = X$expr
  N = nrow(X)
  P = ncol(X)
  message(paste0("## Loading raw data of ",N , "*", P, " and labels...\n"))
  
  
  ## filter zero gene & cell and remain high_var gene
  ## as well as normalization & log-transform
  message("## Data nomalization and log-transformation...\n")
  lg_X = data_normalize(X, N, P, gt, mu_probs=0.2, cv_probs=0.2)    # return log_count without pseudo count
  gt = lg_X$gt
  lg_X = as.matrix(lg_X$count_hv)
  #lg_X = informative_gene_selection(lg_X, delta=0.7)
  rm(X)
  
  ##* cluster number estimation
  #K = cluster_number_estimation(lg_X)
  K = length(unique(gt))
  message("## Estimated cluster number:", K, "\n")
  
  
  ## dimension reduction(npc) calculation
  npc = npc_cal(as.matrix(lg_X), K, var_thre=0.8)
  message("## Caluculated dim reduction npc:", npc, "\n")
  
  
  ## vars initialization
  library(Matrix)
  library(parallel)
  library(aricode)
  
  
  ## cluster update
  # clustering update, 200 iteration by default
  #* stop update when Jbefore < Jafter
  files = list.files("./scImpute/", pattern="*.R")
  sapply(paste0("./scImpute/", files), source)
  
  
  res <- var_update(lg_X, K, npc, lambda1=10, lambda2=10, 
                    iteration=1, clust_iteration=100, imp_iteration=1000, 
                    res_save=FALSE)
  clust = res$cluster
  nmi = NMI(clust, gt)
  ari = ARI(clust, gt)
  res$nmi = nmi
  res$ari = ari
  message("## nmi: ", nmi, "   art: ", ari, "\n")
    
  
  time = Sys.time() - s
  message("## Consume", time, "seconds.\n")
  
  #res = list(clust=clust, nmi=nmi, ari=ari, time=time_consume, H = H)
  output = strsplit(datas[i], split="/")[[1]][1]
  
  saveRDS(res, paste0("result/", output, '.rds'))
#})
}
