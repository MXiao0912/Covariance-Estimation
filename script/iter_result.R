sim_full_rho <- function(cov, n_, p, iter, dist="normal"){

  
  sim_err <- function(iter){
    
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
    
  
    # OASD (designed for diagonal target,so don't apply it to corr=TRUE)
    phi1 = (tr_s2-tr_sdiag2)/(tr_s2+tr_s^2-2*tr_sdiag2)
    rou_oas_variant = min(1/((n_+1)*phi1),1)
   
    # ORACLE (the best oracle is the one that uses the diagonal target directly because the weight optimizes the evaluation criteria)
    # (designed for diagonal target,so don't apply it to corr=TRUE)
    rou_o_variant = ((1/n_)*tr_t2-(2/n_)*tr_tdiag2+(1/n_)*tr_t^2)/(((n_+1)/n_)*tr_t2+(1/n_)*tr_t^2-((n_+2)/n_)*tr_tdiag2)
   
    
    coef_ = data.frame(matrix(nrow=1,ncol=2))
    colnames(coef_) = c('OASD','Oracle')
    coef_[1,'OASD_phi'] = phi1
    coef_[1,'OASD_rho'] = rou_oas_variant
    coef_[1,'Oracle'] = rou_o_variant
    
    return(coef_)
  }
  
  res = map_dfr(seq(1,iter), sim_err)
  
  # res = res %>% mutate(err = (OASD-Oracle)^2)
  # 
  # res = mean(res$err)
  
  return(list(mean(res$OASD), sd(res$OASD), mean(res$Oracle), sd(res$Oracle)))
  # return(res)
} 