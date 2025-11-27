lw <- function(sample){
  n_ = nrow(sample)
  p = ncol(sample)
  scov = cov(sample)
  f1 = (sum(diag(scov))/p)*diag(p)
  f2 = diag(diag(scov))
  tr_s = sum(diag(scov))
  tr_s2 = sum(diag(t(scov) %*% scov))
  tr_sdiag2 = sum(diag(diag(diag(scov)) %*% diag(diag(scov))))
  
  num = 0
  for (i in 1:n_){
    fill = sum((matrix(sample[i,],ncol=1) %*% matrix(sample[i,],nrow=1) - scov)^2)
    num = num + fill
  }
  den = (n_^2)*(tr_s2-(tr_s^2)/p)
  rou_lw = num/den
  rou_lw = min(rou_lw,1)
  
  return(list(coef=rou_lw, covest=(1-rou_lw)*scov+rou_lw*f1))
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
  
  phi = (tr_s2-(tr_s^2)/p)/(tr_s^2+(1-2/p)*tr_s2)
  rou_oas_constant = min(1/((n_+1-2/p)*phi),1)
  
  return(list(coef=rou_oas_constant, covest=(1-rou_oas_constant)*scov+rou_oas_constant*f1))
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
  
  phi1 = (tr_s2-tr_sdiag2)/(tr_s2+tr_s^2-2*tr_sdiag2)
  rou_oas_variant = min(1/((n_+1)*phi1),1)
  
  return(list(coef=rou_oas_variant, covest=(1-rou_oas_variant)*scov+rou_oas_variant*f2))
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
  
  return(list(coef=gamma_oas, covest=(1-theta)*scov+theta*(gamma_oas*f2+(1-gamma_oas)*f1)))
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
  
  sample_sc = sample %*% (diag(diag(f2)^(-1/2)))
  v <- (1/(n_*(n_ - 1))) * (crossprod(sample_sc^2) - (1/n_) * (crossprod(sample_sc))^2)
  diag(v) <- 0
  d <- (scorr - diag(diag(scorr)))^2
  rou_ss <- sum(v)/sum(d)
  rou_ss <- max(min(rou_ss, 1), 0)
  scorr_ss = (1-rou_ss)*scorr+rou_ss*diag(p)
  # SS is always based on corr=TRUE
  
  return(list(coef=rou_ss, covest=(f2^(1/2)) %*% scorr_ss %*% (f2^(1/2))))
}

od <- function(sample, cov){
  n_ = nrow(sample)
  p = ncol(sample)
  scov = cov(sample)
  f2 = diag(diag(scov))
  
  tr_t = sum(diag(cov))
  tr_t2 = sum(diag(t(cov) %*% cov))
  tr_tdiag2 = sum(diag(diag(diag(cov)) %*% diag(diag(cov))))
  
  P = tr_tdiag2
  R = tr_t2
  Q = tr_t^2
  
  rou_o_variant = (Q+R-2*P)/((n_+1)*R+Q-(n_+2)*P)

  return(list(coef=rou_o_variant, covest=(1-rou_o_variant)*scov+rou_o_variant*f2))
}

ob <- function(sample, cov){
  n_ = nrow(sample)
  p = ncol(sample)
  scov = cov(sample)
  f1 = (sum(diag(scov))/p)*diag(p)
  f2 = diag(diag(scov))
  
  tr_t = sum(diag(cov))
  tr_t2 = sum(diag(t(cov) %*% cov))
  tr_tdiag2 = sum(diag(diag(diag(cov)) %*% diag(diag(cov))))
  
  P = tr_tdiag2
  R = tr_t2
  Q = tr_t^2
  
  rou_b_o = (Q+R-2*P)/((n_+1)*R+Q-(n_+2)*P)
  gamma_b_o = 1-(1/rou_b_o)*(((2/n_)*P-(2/(n_*p))*R)/(((n_+2)/n_)*P-(2/(n_*p))*R-(1/p)*Q))
  
  scov_o_both = (1-rou_b_o)*scov+rou_b_o*(gamma_b_o*f2+(1-gamma_b_o)*f1)
  
  return(list(coef=gamma_b_o, covest=(1-rou_b_o)*scov+rou_b_o*(gamma_b_o*f2+(1-gamma_b_o)*f1)))
}