

### loading dataset
setwd("E:/Project/Paper_debug/Clustering algorithm/MGGE/")
data = readRDS("./dataset/mouse_embryo.rds")
gt = data$label[which(data$singlecell)]
# row oriented data
data$expr = data$reads[which(data$singlecell), ]
data$datatype = "count"


if(FALSE){
  ## running SAFE
  message("# Performing SAFE clustering...")
  SAFE_PATH = "E:/Project/Paper_debug/Bioinformatics/sigle cell RNA Clustering/SAFE/R"
  files = list.files(SAFE_PATH)
  for (file in files){
    source(paste(SAFE_PATH, file, sep="/"))
  } 
}


library(Seurat)
library(dplyr)
library(Matrix)
library(methods)



## redifine
cidr <- function(inputTags, datatype, nPC.cidr, nCluster, SEED){
  message("# Performing CIDR Clustering...")
  library(cidr)
  
  cidrOUTPUT <- NULL
  
  if (datatype == "count"){
    cidrOUTPUT <- scDataConstructor(as.matrix(t(inputTags)), tagType = "raw")
  }else if (datatype == "CPM" || datatype == "FPKM" || datatype == "RPKM" || datatype == "TPM"){
    cidrOUTPUT <- scDataConstructor(as.matrix(t(inputTags)), tagType = "cpm")
  }
  cidrOUTPUT <- determineDropoutCandidates(cidrOUTPUT)
  cidrOUTPUT <- wThreshold(cidrOUTPUT)
  cidrOUTPUT <- scDissim(cidrOUTPUT)
  cidrOUTPUT <- scPCA(cidrOUTPUT)
  cidrOUTPUT <- nPC(cidrOUTPUT)
  
  ### Define nPC
  if(!is.null(nPC.cidr)) {
    cidrOUTPUT@nPC <- nPC.cidr
  } else {
    nPC.cidr <- cidrOUTPUT@nPC
  }
  
  ### Clustering by CIDR.
  # The optimal clustering number is determined automatically
  cidrOUTPUT <- scCluster(cidrOUTPUT, nPC = nPC.cidr, nCluster=nCluster)
  return(cidrOUTPUT)  
}


SC3 <- function(inputTags, datatype, gene_filter, svm_num_cells, SEED){
  # inputTags P*N matrix
  library(SingleCellExperiment)
  library(SC3)
  library(e1071)
  library(parallel)
  library(doParallel)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(doRNG)
  library(foreach)
  
  exp_cell_exprs <- NULL
  sc3OUTPUT <- NULL
  
  # cell expression
  if (datatype == "count") {
    ### For count data, it would be normalized by the total cound number and then log2 transformed
    exp_cell_exprs <- SingleCellExperiment(assays = list(counts = inputTags))
    normcounts(exp_cell_exprs) <- t(t(inputTags)/colSums(inputTags))*1000000
    logcounts(exp_cell_exprs) <- log2(normcounts(exp_cell_exprs) + 1)
  } else if (datatype == "CPM" || datatype == "FPKM" || datatype == "RPKM" || datatype == "TPM") {
    ### For CPM, FPKM, RPKM or TPM data, it would be log2 transformed
    exp_cell_exprs <- SingleCellExperiment(assays = list(normcounts = inputTags))
    logcounts(exp_cell_exprs) <- log2(normcounts(exp_cell_exprs) + 1)
  }
  
  rowData(exp_cell_exprs)$feature_symbol <- rownames(exp_cell_exprs)
  exp_cell_exprs <- exp_cell_exprs[!duplicated(rowData(exp_cell_exprs)$feature_symbol), ]
  
  ### Estimating optimal number of clustering
  exp_cell_exprs <- sc3_estimate_k(exp_cell_exprs)
  optimal_K <- metadata(exp_cell_exprs)$sc3$k_estimation
  
  ### Clustering by SC3 at the optimal K
  if (ncol(inputTags) < svm_num_cells){
    #print(optimal_K)
    exp_cell_exprs <- sc3(exp_cell_exprs, ks = optimal_K, biology = FALSE, gene_filter = gene_filter, n_cores = 1, rand_seed = SEED)
  } else if (ncol(inputTags) >= svm_num_cells){
    ### Runing SVM
    exp_cell_exprs <- sc3(exp_cell_exprs, ks = optimal_K, biology = FALSE, gene_filter = gene_filter,
                          svm_max = svm_num_cells, svm_num_cells = svm_num_cells, n_cores = 1, rand_seed = SEED)
    exp_cell_exprs <- sc3_run_svm(exp_cell_exprs, ks = optimal_K)
  }
  
  ### Exporting SC3 results
  p_Data <- colData(exp_cell_exprs)
  col_name <- paste("sc3_", optimal_K, "_clusters", sep = '')
  sc3OUTPUT <- p_Data[, grep(col_name, colnames(p_Data))]
  return(sc3OUTPUT)
}


tSNE_kmeans <- function(inputTags, datatype, saver, dimensions, perplexity=30, ncluster, var_genes, SEED){
  library(Rtsne)
  library(ADPclust)
  library(SAVER)
  library(doParallel)
  library(S4Vectors)
  library(stats)
  input_lcpm <- NULL
  tsne_input <- NULL
  tsne_output <- NULL
  tsne_kmeansOUTPUT <- NULL
  adpOUTPUT <- NULL
  
  ### Data tranformation
  if (datatype == "count") {
    if(saver == TRUE){
      cl <- makeCluster(12, outfile = "")
      registerDoParallel(cl)
      inputTags_corrected <- saver(inputTags)
      stopCluster(cl)
      tsne_input <- log2(inputTags_corrected$estimate + 1)
    } else{
      ### If the input data is original count data or CPM, it would be tranformed to CPM
      input_lcpm <- log2(t(t(inputTags)/colSums(inputTags))*1000000+1)
      tsne_input <- input_lcpm
    }
  } else if (datatype == "CPM" || datatype == "FPKM" || datatype == "RPKM" || datatype == "TPM") {
    ### If the input data is FPKM or RPKM, we use the transformed TPM data generated before as the input
    tsne_input <- log2(inputTags + 1)
  }
  
  if (is.null(var_genes)){
    set.seed(SEED)
    
    tsne_output <- Rtsne(t(tsne_input), dims = dimensions, perplexity = perplexity, check_duplicates = FALSE)
  } else{
    se_genes = rep(NA, nrow(tsne_input))
    for (i in 1:nrow(tsne_input)){
      se_genes[i] = sqrt(var(tsne_input[i,])/length(tsne_input[i,]))
    }
    decreasing_rank = order(se_genes, decreasing = TRUE)
    
    set.seed(SEED)
    
    tsne_output <- Rtsne(t(tsne_input[decreasing_rank[1:var_genes],]), dims = dimensions, perplexity = perplexity)
  }
  
  ### Determining the optimal cluster number (k) and centroid by ADPclust
  adpOUTPUT <- adpclust(tsne_output$Y, htype = "amise",centroids="auto", nclust = ncluster)
  
  ### Clustering the cells by kmeans
  tsne_kmeansOUTPUT <- kmeans(tsne_output$Y, tsne_output$Y[adpOUTPUT$centers,], adpOUTPUT$nclust)
  
  return(tsne_kmeansOUTPUT)
}


seurat_SAFE <- function(inputTags, datatype, nPC.seurat, resolution, seurat_min_cell, resolution_min, SEED){
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(methods)
  
  seuratOUTPUT <- NULL
  
  # Initialize the Seurat object with the raw data (non-normalized data)
  # Keep all genes expressed in >= 3 cells, keep all cells with >= 200 genes
  #seuratOUTPUT <- CreateSeuratObject(raw.data = inputTags, min.cells = 3, min.genes = 200, project = "single-cell clustering")
  seuratOUTPUT <- CreateSeuratObject(raw.data = inputTags, min.cells = 3, min.features = 200, project = "single-cell clustering")
  
  # Perform log-normalization, first scaling each cell to a total of 1e4 molecules (as in Macosko et al. Cell 2015)
  if (datatype == "count"){
    seuratOUTPUT = NormalizeData(object = seuratOUTPUT, normalization.method = "LogNormalize", scale.factor = 10000)
  } else if (datatype == "CPM" || datatype == "FPKM" || datatype == "RPKM" || datatype == "TPM"){
    raw.data <- GetAssayData(object = seuratOUTPUT, slot = "raw.data")
    normalized.data <- log(raw.data+1)
    colnames(x = normalized.data) <- colnames(x = raw.data)
    rownames(x = normalized.data) <- rownames(x = raw.data)
    seuratOUTPUT <- SetAssayData(object = seuratOUTPUT, assay.type = "RNA",slot = "data", new.data = normalized.data)
  }
  
  # Detection of variable genes across the single cells
  seuratOUTPUT = FindVariableGenes(object = seuratOUTPUT, mean.function = ExpMean, dispersion.function = LogVMR,
                                   x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  
  # Regress out unwanted sources of variation
  seuratOUTPUT <- ScaleData(object = seuratOUTPUT, vars.to.regress = c("nUMI"))
  
  ### Perform linear dimensional reduction
  if (nPC.seurat <= 20){
    seuratOUTPUT <- RunPCA(object = seuratOUTPUT, pc.genes = seuratOUTPUT@var.genes, do.print = FALSE)
  } else {
    seuratOUTPUT <- RunPCA(object = seuratOUTPUT, pc.genes = seuratOUTPUT@var.genes, pcs.compute = nPC.seurat, do.print = FALSE)
  }
  
  if (length(inputTags[1,]) >= seurat_min_cell){
    ### Determine statistically significant principal components
    # NOTE: This process can take a long time for big datasets, comment out for expediency.
    # More approximate techniques such as those implemented in PCElbowPlot() can be used to reduce computation time
    # Here we chooes the same number of PCs used in CIDR
    
    ### Clustering the cells by Seurat
    seuratOUTPUT <- FindClusters(object = seuratOUTPUT, reduction.type = "pca", dims.use = 1:nPC.seurat, algorithm = 3,
                                 resolution = resolution, print.output = FALSE, random.seed = SEED)
  } else {
    resolution <- resolution_min
    seuratOUTPUT <- FindClusters(object = seuratOUTPUT, reduction.type = "pca", dims.use = 1:nPC.seurat, algorithm = 3,
                                 resolution = resolution_min, print.output = FALSE, random.seed = SEED)
  }
  
  ### Complementing the missing data
  cells_dropout <- NULL
  num_genes <- colSums(inputTags > 0)
  cells_dropout <- names(num_genes[which(num_genes <= 200)])
  if (length(cells_dropout != 0)){
    seurat_output <- matrix(NA, ncol = ncol(inputTags), byrow = TRUE)
    colnames(seurat_output) <- colnames(inputTags)
    seurat_retained <- t(as.matrix(as.numeric(seuratOUTPUT@ident)))
    colnames(seurat_retained) <- colnames(seuratOUTPUT@data)
    for (i in 1:ncol(seurat_retained)){
      seurat_output[1,colnames(seurat_retained)[i]] <- seurat_retained[1,colnames(seurat_retained)[i]]
    }
  } else {
    seurat_output <- t(as.matrix(as.numeric(seuratOUTPUT@ident)))
  }
  
  return(seurat_output)
}




# running CIDR

cidrOUTPUT <- cidr(inputTags = as.matrix(data$expr), datatype = data$Unit, nPC.cidr = NULL, SEED = 1234)


# running SC3

sc3OUTPUT = SC3(inputTags = data$expr, datatype = data$Unit, gene_filter = FALSE, svm_num_cells = 5000, SEED = SEED)
cluster_number = cidrOUTPUT@nCluster

# running tSNE-Kmeans

tSNEKmeansOUTPUT = tSNE_kmeans(inputTags = data$expr, datatype = data$Unit, saver = FALSE, dimensions = 3,
                               perplexity = 10, k.min = 2, k.max = max(cluster_number), var_genes = NULL, SEED = SEED)

res = list(cidr = cidrOUTPUT, sc3 = sc3OUTPUT, tsne = tSNEKmeansOUTPUT)

if(FALSE){
  # running Seurat
  seurat_output <- seurat_SAFE(inputTags = data$expr, datatype = data$Unit, nPC.seurat = 272, resolution = 0.9,
                               seurat_min_cell = 200, resolution_min = 1.2, SEED = SEED)

  
}


