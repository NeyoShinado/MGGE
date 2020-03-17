setwd("D:/Work/Project/Paper_algorithms/sc-RNA/scImpute-master")

library(parallel)
library(kernlab)
library(doParallel)

ans = list.files("./R")
ans = paste0("R/", ans[1:9])

sapply(ans, source)

s = Sys.time()
scimpute(count_path="inst/extdata/raw_count.csv", infile="csv", 
         outfile="csv", out_dir="./output", labeled=FALSE, drop_thre=0.5, 
         Kcluster=2, ncores=1)
e = Sys.time()
print(e-s)
