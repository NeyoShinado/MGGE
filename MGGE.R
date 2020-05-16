# demo program of MGGE(Multi-Granularity Graph Embedding)
# (row-oriented)

## Input
# x -- imputed or origin log-tranformed N*P data matrix
# D -- dropout weighted fix diagonal matrix for N cells
# L -- integrated laplacian graph of global clustered-sim & local imputed-sim
# K -- number of cluster
# npc -- number of reduced dims decided by data var cumsum
# lambda1 -- par of low dim feature embedding  
# lambda2 -- par of integrated graph embedding
# imputed_flag -- pointer of imputed data type or unimputed data type
## Output

# load functions
setwd("E:/Project/Paper_debug/Clustering algorithm/MGGE/")
source("data_normalize.R")
source("npc_cal.R")
source("localsim_integrated_update.R")
source("imputing_update.R")
#source("var_update.R")


# dataset info
files = c("Biase/Biase.rds", "Deng_GSE45719/mouse_embryo.rds", "Zeisel/Zeisel.rds", 
          "Baron/mouse1.rds", "Baron/mouse2.rds", "Baron/human3.rds", "Hrvatin/Hrvatin.rds")
data_path = "~/Data/sc-data/"


# Perprocess
## loading raw data(mouse embryo by default)
res = lapply(files, function(i){
  tryCatch(
    function(i){s = Sys.time()
    
    X = readRDS(paste0(data_path, i))
    gt = X$gt
    X = X$expr
    N = nrow(X)
    P = ncol(X)
    cat(paste0("## Loading raw data of ",N , "*", P, " and labels...\n"))
    
    
    ## filter zero gene & cell and remain 75% high_var gene
    ## as well as normalization & log-transform
    cat("## Data nomalization and log-transform...\n")
    lg_X = data_normalize(X, N, P, gt, mu_probs=0.4, cv_probs=0.25)  
    gt = lg_X$gt
    lg_X = as.matrix(lg_X$count_hv)
    lg_X = Matrix(lg_X, sparse=TRUE)
    rm(X)
    
    
    ##* cluster number estimation
    #K = cluster_number_estimation(lg_X)
    K = length(unique(gt))
    cat("## Estimated cluster number:", K, "\n")
    
    
    ## dimension reduction(npc) calculation
    npc = npc_cal(as.matrix(lg_X), K, var_thre=0.6)
    cat("## Caluculated dim reduction npc:", npc, "\n")
    
    
    ## vars initialization
    library(bignmf)
    library(Matrix)
    library(parallel)
    
    
    ## cluster update
    # clustering update, 200 iteration by default
    #* stop update when Jbefore < Jafter
    files = list.files("./scImpute/", pattern="*.R")
    sapply(paste0("./scImpute/", files), source)
    
    
    time_consume = system.time({
      res <- var_update(lg_X, K, npc, lambda1=lambda1, lambda2=lambda2, 
                        iteration=1, clust_iteration=1000, imp_iteration=3000, 
                        res_save=FALSE)
      clust = res$cluster
      nmi = NMI(clust, gt)
      ari = ARI(clust, gt)
      
    })
    cat("## Consume", time_consume[3], "seconds.\n")
    
    time = Sys.time() - s
    return(list(clust=clust, nmi=nmi, ari=ari, time=time_consume))
    }
  , error=function(e){return(e)}
  )
  
})

names(res) = c("Biase", "Deng", "Zeisel", "mouse1", "mouse2", "human3", "Hrvatin")
saveRDS(res, "result/res.rds")
