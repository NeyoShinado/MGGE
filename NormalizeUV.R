# normalize updating vars

NormalizeUV <- function(U, V, Norm=2){
  # Object term ||X - UV||,  which X ~ N*P, noramlize U by default
  # Norm -- 2: Fro Norm
  #      -- 1: 1 Norm
  U = as.matrix(U)
  V = as.matrix(V)
  N = dim(U)[1]
  P = dim(V)[2]
  k = dim(U)[2]
  
  if(Norm != 1 && Norm != 2){
    stop("## Incorrect nrom pars for norm type identification!")
  }
  

  if(Norm == 2){
    norms = sqrt(colSums(U^2))
  }else{
    norms = sqrt(colSums(abs(U)))
  }
  
  norms[norms == 0] = 1e-10
  U = U / t(matrix(norms, k, N))
  V = V * matrix(norms, k, P)
  
  res = list(U = U, V = V)
  return(res)
}