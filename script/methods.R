lw <- function(sample){
  n_ = nrow(sample)
  p = ncol(sample)
  scov = cov(sample)
  f1 = (sum(diag(scov))/p)*diag(p)
  f2 = diag(diag(scov))
  tr_s = sum(diag(scov))
  tr_s2 = sum(diag(t(scov) %*% scov))
  tr_sdiag2 = sum(diag(diag(diag(scov)) %*% diag(diag(scov))))
  
  scorr = cov2cor(scov)
  f1_cor = (sum(diag(scorr))/p)*diag(p)
  sample_sc = sample %*% (diag(diag(f2)^(-1/2)))
  
  num = 0
  for (i in 1:n_){
    fill = sum((matrix(sample[i,],ncol=1) %*% matrix(sample[i,],nrow=1) - scov)^2)
    num = num + fill
  }
  den = (n_^2)*(tr_s2-(tr_s^2)/p)
  rou_lw = num/den
  rou_lw = min(rou_lw,1)
  
  return((1-rou_lw)*scov+rou_lw*f1)
}

oas <- function(sample){
  n_ = nrow(sample)
  p = ncol(sample)
  scov = cov(sample)
  f1 = (sum(diag(scov))/p)*diag(p)
  f2 = diag(diag(scov))
  tr_s = sum(diag(scov))
  tr_s2 = sum(diag(t(scov) %*% scov))
  tr_sdiag2 = sum(diag(diag(diag(scov)) %*% diag(diag(scov))))
  
  scorr = cov2cor(scov)
  f1_cor = (sum(diag(scorr))/p)*diag(p)
  sample_sc = sample %*% (diag(diag(f2)^(-1/2)))
  
  phi = (tr_s2-(tr_s^2)/p)/(tr_s^2+(1-2/p)*tr_s2)
  rou_oas_constant = min(1/((n_+1-2/p)*phi),1)
  return((1-rou_oas_constant)*scov+rou_oas_constant*f1)
}

oasd <- function(sample){
  n_ = nrow(sample)
  p = ncol(sample)
  scov = cov(sample)
  f1 = (sum(diag(scov))/p)*diag(p)
  f2 = diag(diag(scov))
  tr_s = sum(diag(scov))
  tr_s2 = sum(diag(t(scov) %*% scov))
  tr_sdiag2 = sum(diag(diag(diag(scov)) %*% diag(diag(scov))))
  
  scorr = cov2cor(scov)
  f1_cor = (sum(diag(scorr))/p)*diag(p)
  sample_sc = sample %*% (diag(diag(f2)^(-1/2)))
  
  phi1 = (tr_s2-tr_sdiag2)/(tr_s2+tr_s^2-2*tr_sdiag2)
  rou_oas_variant = min(1/((n_+1)*phi1),1)
  return((1-rou_oas_variant)*scov+rou_oas_variant*f2)
}

oasb <- function(sample){
  n_ = nrow(sample)
  p = ncol(sample)
  scov = cov(sample)
  f1 = (sum(diag(scov))/p)*diag(p)
  f2 = diag(diag(scov))
  tr_s = sum(diag(scov))
  tr_s2 = sum(diag(t(scov) %*% scov))
  tr_sdiag2 = sum(diag(diag(diag(scov)) %*% diag(diag(scov))))
  
  scorr = cov2cor(scov)
  f1_cor = (sum(diag(scorr))/p)*diag(p)
  sample_sc = sample %*% (diag(diag(f2)^(-1/2)))
  
  phi1 = (tr_s2-tr_sdiag2)/(tr_s2+tr_s^2-2*tr_sdiag2)
  rou_oas_variant = min(1/((n_+1)*phi1),1)
  theta = rou_oas_variant
  A = (tr_s^2)/p
  B = tr_s2-tr_sdiag2
  C = tr_sdiag2-A
  
  t1 = ((p-1)*C)/((p-1)*(A+C)-B)
  t2 = (B-(p-1)*C)/((p-1)*(A+C)-B)
  t3 = (n_*p*C)/(2*(p-1)*(A+C)-2*B)
  
  if (abs(t3/((t1+t2)*theta+1-t1))<1){
    gamma_oas = (theta-1)/(theta)
  }else{
    gamma_oas = (theta*(t3-t2)-1)/(theta*(t1+t3))
  }
  
  return((1-theta)*scov+theta*(gamma_oas*f2+(1-gamma_oas)*f1))
}

shafer <- function(sample){
  n_ = nrow(sample)
  p = ncol(sample)
  scov = cov(sample)
  f1 = (sum(diag(scov))/p)*diag(p)
  f2 = diag(diag(scov))
  tr_s = sum(diag(scov))
  tr_s2 = sum(diag(t(scov) %*% scov))
  tr_sdiag2 = sum(diag(diag(diag(scov)) %*% diag(diag(scov))))
  
  scorr = cov2cor(scov)
  f1_cor = (sum(diag(scorr))/p)*diag(p)
  sample_sc = sample %*% (diag(diag(f2)^(-1/2)))
  
  v <- (1/(n_*(n_ - 1))) * (crossprod(sample_sc^2) - (1/n_) * (crossprod(sample_sc))^2)
  diag(v) <- 0
  d <- (scorr - diag(diag(scorr)))^2
  rou_ss <- sum(v)/sum(d)
  rou_ss <- max(min(rou_ss, 1), 0)
  scorr_ss = (1-rou_ss)*scorr+rou_ss*f1_cor
  # SS is always based on corr=TRUE
  return((f2^(1/2)) %*% scorr_ss %*% (f2^(1/2)))
}