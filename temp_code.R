### NMI
clust = specc(as.matrix(res$W), 12)@.Data
nmi = NMI(clust, gt)


L = diag(colSums(S)) - S
clust = specc(L, K)@.Data
NMI(clust, gt)


H = eigen(L)$vectors
H = H[, (N-K):(N-1)]
clust = kmeans(H, K)$cluster
NMI(clust, gt)




# droprate preprose
step = 4
map = c(0, c(1:step) * (1/step))
droprate = sapply(1:dim(droprate)[2], function(i){
  temp = droprate[, i]
  for(j in 2:step){
    temp[which(temp <= map[j]) && temp > map[j-1]] = map[j]
  }
  return(temp)
})


### droprate based on group
map = unique(gt)
group_id = lapply(1:length(map), function(i){
  res = which(gt == map[i])
  return(res)
})
group_droprate = sapply(group_id, function(i){
  clust_droprate = colSums(droprate[i, ] / length(i))
  return(clust_droprate)
})
#group_droprate = as.data.frame(group_droprate)
#group_droprate = cbind(group_droprate, data.frame(group=map))

# visualize
#test = melt(group_droprate, ID = colnames(group_droprate)[1:(ncol(group_droprate)-1)], group_id="group")
#ggplot(test, aes(variable, value), group = group) + geom_line()

sort_geneid = order(group_droprate[, 1], decreasing=TRUE)
for(i in 1:length(map)){
#  png(paste0("result/figure/gene mean drop by cluster/gene mean droprate of cluster ", map[i], ".png"), width=600, height=600)
  dev.new()
  plot(group_droprate[sort_geneid, i], main=paste0("gene\'s mean droprate of cluster ", map[i]),
       xlab = paste0("gene\'s mean droprate"), ylab = paste0("sorted gene index of cluster ", map[1]))
#  dev.off()
}




### cell's droprate based on group
Nsp = 2
group = 1
sample_id = sample(group_id[[group]], Nsp, replace=FALSE)
sort_geneid = order(group_droprate[, 2], decreasing=TRUE)
for(i in 1:Nsp){
  if(i == 1){
    #png(paste0("result/figure/cell drop by cluster/mean droprate of cluster ", map[group], ".png"), width=600, height=600)
    dev.new()
    plot(group_droprate[sort_geneid, group], main=paste0("mean droprate of cluster ", map[group], " -- index based on cluster ", map[1]),
         xlab = paste0("sorted gene index of cluster ", map[1]), ylab = paste0("mean droprate of cluster ", map[group]))
    #dev.off()
  }
  
  #png(paste0("result/figure/cell drop by cluster/cell droprate of cell ", i, " from cluster ", map[group], ".png"), width=600, height=600)
  dev.new()
  plot(droprate[sample_id[i], sort_geneid], main=paste0("droprate of cell ", i, " from cluster ", map[group]),
       xlab = paste0("sorted gene index of cluster ", map[1]), ylab = paste0("droprate of cell ", i, " from cluster ", map[group]))
  #dev.off()
}




### reading weighted res
nmi = matrix(0, 5, 5)
rownames(nmi) = as.character(c(.5, 1, 2, 3, 4))
colnames(nmi) = as.character(c(.5, 1, 2, 3, 4))

for(sigmac in c(.5, 1, 2, 3, 4)){
  for(sigmag in c(.5, 1, 2, 3, 4)){
    try({
      res = readRDS(paste0("result/weighted/mouse2c_", sigmac, '_g_', sigmag,'.rds'))
      res = res$nmi[301, 1]
      nmi[as.character(sigmac), as.character(sigmag)] = res
    })
  }
}


### visualize
dev.new()
plot(res$lambda1)
dev.new()
plot(res$lambda2)
dev.new()
plot(res$J_HE)
dev.new()
plot(res$J_LE)
dev.new()
plot(res$nmi[,1])
title("nmi_H")
dev.new()
plot(res$nmi[,2])
title("nmi_WL")
dev.new()
plot(res$nmi[,3])
title("nmi_W_sp")