proj_sim <- function(ad){
  k = 1
  ft = 1
  n = length(ad)
  ad0 = ad - mean(ad) + 1/n

  admin = min(ad0)
  if(admin < 0){
    f = 1
    lambda_m = 0
    while(abs(f) > 10^-10){
      ad1 = ad0 - lambda_m
      posidx = ad1 > 0
      npos = sum(posidx)
      g = -npos
      f = sum(ad1[posidx]) - k
      lambda_m = lambda_m - f/g
      ft = ft + 1
      if(ft > 100){
        x = max(ad1, 0)
        break
      }
    }
    x = max(ad1, 0)
  }else{
    x = ad0
  }
  return(x)
}