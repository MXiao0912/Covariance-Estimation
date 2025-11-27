sim_full <- function(sd_, r_, n_, iter, corr_ind=FALSE, inv_=FALSE, dist="normal"){
  
  corr = matrix(0, ncol=p, nrow=p)
  for (i in 1:p){
    for (j in 1:p){
      corr[i,j]=r_^abs(i-j)
    }
  }
  
  sd = diag(c(rep(small,p/2),rep(sd_,p/2)))
  cov = sd %*% corr %*% sd
  
  sim_err <- function(iter_){
    
    # Generate samples from n(0,cov)
    if (dist=="normal"){
      sample = mvrnorm(n_, rep(0,p), cov)
    }else if (dist=="t"){
      sample=rmvt(n = n_, sigma = (3/5)*cov, df = 5, delta = rep(0,p))
    }else if (dist=="uniform"){
      sample <- matrix(runif(n_*p, min=-sqrt(3), max=sqrt(3)), ncol=p)
      L <- chol(cov)
      sample = sample %*% L
    }
    
    # Estimate the covariance matrix
    scov = (t(sample) %*% sample)/nrow(sample)
    f1 = (sum(diag(scov))/p)*diag(p)
    f2 = diag(diag(scov))
    tr_s = sum(diag(scov))
    tr_s2 = sum(diag(t(scov) %*% scov))
    tr_sdiag2 = sum(diag(diag(diag(scov)) %*% diag(diag(scov))))
    
    tr_t = sum(diag(cov))
    tr_t2 = sum(diag(t(cov) %*% cov))
    tr_tdiag2 = sum(diag(diag(diag(cov)) %*% diag(diag(cov))))
    
    scorr = cov2cor(scov)
    f1_cor = (sum(diag(scorr))/p)*diag(p)
    sample_sc = sample %*% (diag(diag(f2)^(-1/2)))
    
    if(corr_ind){
      tr_s_cor = sum(diag(scorr))
      tr_s2_cor = sum(diag(t(scorr) %*% scorr))
    }
    
    # LW
    num = 0
    if (!corr_ind){
      for (i in 1:n_){
        fill = sum((matrix(sample[i,],ncol=1) %*% matrix(sample[i,],nrow=1) - scov)^2)
        num = num + fill
      }
      den = (n_^2)*(tr_s2-(tr_s^2)/p)
    }else{
      for (i in 1:n_){
        fill = sum((matrix(sample_sc[i,],ncol=1) %*% matrix(sample_sc[i,],nrow=1) - scorr)^2)
        num = num + fill
      }
      den = (n_^2)*(tr_s2_cor-(tr_s_cor^2)/p)
    }
    rou_lw = num/den
    rou_lw = min(rou_lw,1)
    
    if (corr_ind){
      scorr_lw = (1-rou_lw)*scorr+rou_lw*f1_cor
      scov_lw = (f2^(1/2)) %*% scorr_lw %*% (f2^(1/2))
    }else{
      scov_lw = (1-rou_lw)*scov+rou_lw*f1
    }
    
    
    # OAS
    if (!corr_ind){
      phi = (tr_s2-(tr_s^2)/p)/(tr_s^2+(1-2/p)*tr_s2)
    }else{
      phi = (tr_s2_cor-(tr_s_cor^2)/p)/(tr_s_cor^2+(1-2/p)*tr_s2_cor)
    }
    rou_oas_constant = min(1/((n_+1-2/p)*phi),1)
    
    if (corr_ind){
      scorr_oas_constant = (1-rou_oas_constant)*scorr+rou_oas_constant*f1_cor
      scov_oas_constant = (f2^(1/2)) %*% scorr_oas_constant %*% (f2^(1/2))
    }else{
      scov_oas_constant = (1-rou_oas_constant)*scov+rou_oas_constant*f1
    }
    
    
    # OASD (designed for diagonal target,so don't apply it to corr=TRUE)
    phi1 = (tr_s2-tr_sdiag2)/(tr_s2+tr_s^2-2*tr_sdiag2)
    rou_oas_variant = min(1/((n_+1)*phi1),1)
    scov_oas_variant = (1-rou_oas_variant)*scov+rou_oas_variant*f2
    
    
    # OASB
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
    
    scov_oas_both = (1-theta)*scov+theta*(gamma_oas*f2+(1-gamma_oas)*f1)
    
    # SS (designed for diagonal target,so don't apply it to corr=TRUE)
    v <- (1/(n_*(n_ - 1))) * (crossprod(sample_sc^2) - (1/n_) * (crossprod(sample_sc))^2)
    diag(v) <- 0
    d <- (scorr - diag(diag(scorr)))^2
    rou_ss <- sum(v)/sum(d)
    rou_ss <- max(min(rou_ss, 1), 0)
    scorr_ss = (1-rou_ss)*scorr+rou_ss*f1_cor
    # SS is always based on corr=TRUE
    scov_ss = (f2^(1/2)) %*% scorr_ss %*% (f2^(1/2))
    
    # # ORACLE (the best oracle is the one that uses the diagonal target directly because the weight optimizes the evaluation criteria)
    # # (designed for diagonal target,so don't apply it to corr=TRUE)
    # rou_o_variant = ((1/n_)*tr_t2-(2/n_)*tr_tdiag2+(1/n_)*tr_t^2)/(((n_+1)/n_)*tr_t2+(1/n_)*tr_t^2-((n_+2)/n_)*tr_tdiag2)
    # scov_o_variant = (1-rou_o_variant)*scov+rou_o_variant*f2
    
    # Method 9: Oracle_b
    P = tr_tdiag2
    R = tr_t2
    Q = tr_t^2
    
    rou_b_o = (Q+R-2*P)/((n_+1)*R+Q-(n_+2)*P)
    gamma_b_o = 1-(1/rou_b_o)*(((2/n_)*P-(2/(n_*p))*R)/(((n_+2)/n_)*P-(2/(n_*p))*R-(1/p)*Q))
    
    scov_o_both = (1-rou_b_o)*scov+rou_b_o*(gamma_b_o*f2+(1-gamma_b_o)*f1)
    
    # OD (designed for diagonal target,so don't apply it to corr=TRUE)
    rou_od = ((1/n_)*tr_s2-(2/n_)*tr_sdiag2+(1/n_)*tr_s^2)/(((n_+1)/n_)*tr_s2+(1/n_)*tr_s^2-((n_+2)/n_)*tr_sdiag2)
    scov_od = (1-rou_od)*scov+rou_od*f2
    
    # Collect errors
    if (!inv_){
      err = data.frame(matrix(nrow=1,ncol=8))
      colnames(err) = c('OASD','OAS', 'SS','LW','s','OD','Oracle', 'OASB')
      err[1,'LW'] = get_err(scov_lw,cov)
      err[1,'SS'] = get_err(scov_ss,cov)
      err[1,'OAS'] = get_err(scov_oas_constant,cov)
      err[1,'OASD'] = get_err(scov_oas_variant,cov)
      err[1,'OASB'] = get_err(scov_oas_both,cov)
      err[1,'s']=sum((scov-cov)^2)
      # err[1,'Oracle'] = get_err(scov_o_variant,cov)
      err[1,'Oracle'] = get_err(scov_o_both,cov)
      err[1,'OD'] = get_err(scov_od,cov)
    }else{
      err = data.frame(matrix(nrow=1,ncol=9))
      colnames(err) = c('OASD','OAS', 'SS','LW','s','OD','MP','Oracle','OASB')
      err[1,'LW'] = get_err(solve(scov_lw),solve(cov))
      err[1,'SS'] = get_err(solve(scov_ss),solve(cov))
      err[1,'OAS'] = get_err(solve(scov_oas_constant),solve(cov))
      err[1,'OASD'] = get_err(solve(scov_oas_variant),solve(cov))
      err[1,'OASB'] = get_err(solve(scov_oas_both),solve(cov))
      err[1,'s']=NA
      # err[1,'Oracle'] = get_err(solve(scov_o_variant),solve(cov))
      err[1,'Oracle'] = get_err(solve(scov_o_both),solve(cov))
      err[1,'MP'] = get_err(ginv(scov),solve(cov))
      err[1,'OD'] = get_err(ginv(scov_od),solve(cov))
    }
    
    coef_ = data.frame(matrix(nrow=1,ncol=8))
    colnames(coef_) = c('OASD','OAS', 'SS', 'LW','OD','Oracle', 'Oracle_g','OASB')
    coef_[1,'LW'] = rou_lw
    coef_[1,'SS'] = rou_ss
    coef_[1,'OAS'] = rou_oas_constant
    coef_[1,'OASD'] = rou_oas_variant
    coef_[1,'Oracle'] = rou_b_o
    coef_[1,'Oracle_g'] = gamma_b_o
    coef_[1,'OD'] = rou_od
    coef_[1,'OASB'] = gamma_oas

    return(list(err,coef_))
  }
  
  res = map(seq(1,iter), sim_err)
  res_err = map_dfr(res, function(x){x[[1]]}) %>% 
    mutate(sd = sd_, r=r_, n=n_)
  res_coef = map_dfr(res, function(x){x[[2]]}) %>% 
    mutate(sd = sd_, r=r_, n=n_)
  
  return(list(res_err, res_coef))
}