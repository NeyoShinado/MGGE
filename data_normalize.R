# cells' total count normalization & log-tranform
# X should be a N*P raw-count data matrix

data_normalize <- function(X, N, P, gt, mu_probs=0.4, cv_probs=0.3){

  totalCounts_by_gene = colSums(X)
  totalCounts_by_cell = rowSums(X)
  cell_nzero = setdiff(1:N, which(totalCounts_by_cell==0))
  gene_nzero = setdiff(1:P, which(totalCounts_by_gene==0))
  if(length(cell_nzero) != N){
    cat("### filtered", N-length(cell_nzero), "zero_count cells...\n")
    X = X[cell_nzero, ]
    totalCounts_by_cell = totalCounts_by_cell[cell_nzero]
    gt = gt[cell_nzero]
  }
  if(length(gene_nzero) != P){
    cat("### filtered", P-length(gene_nzero), "zero_count genes...\n")
    X = X[, gene_nzero]
    totalCounts_by_gene = totalCounts_by_gene[gene_nzero]
  }
  
  X = sweep(X, MARGIN = 1, 10^6/totalCounts_by_cell, FUN = "*")
  # normalize genes
  norm_gene = sqrt(colSums(X^2))
  norm_gene = t(matrix(norm_gene, length(norm_gene), N))
  X = X / norm_gene
  
  if (min(X) < 0) {
    stop("smallest read count cannot be negative!")
  }
  lg_X = log10(X + 1)
  # lg_X = log(X + 1.01)
  
  
  # bi-peak genes often correspond to relative-high miu-gene
  # find hv_genes
  count = lapply(1:ncol(lg_X), function(i) lg_X[,i])
  mu = sapply(count, mean)
  mu[is.na(mu)] = 0
  sd = sapply(count, sd)
  sd[is.na(sd)] = 0
  cv = sd/mu
  cv[is.na(cv)] = 0
  # sum(mu >= 1 & cv >= quantile(cv, 0.25), na.rm = TRUE)
  high_var_genes = which(mu >= quantile(mu, mu_probs) & cv >= quantile(cv, cv_probs))
  if(length(high_var_genes) < 500){ 
    high_var_genes = 1:ncol(lg_X)
  }
  count_hv = lg_X[, high_var_genes]
  return(list(count_hv=count_hv, gt=gt))
  
  
}

