source("script/methods.R")

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
    
    scov = cov(sample)
    
    if (corr_ind==FALSE){
      # LW
      lw_res = lw(sample)
      scov_lw = lw_res$covest
      rou_lw = lw_res$coef
      
      # OAS
      oas_res = oas(sample)
      scov_oas = oas_res$covest
      rou_oas = oas_res$coef  
      
    }else{
      
      f2 = diag(diag(cov(sample)))
      sample_sc = sample %*% (diag(diag(f2)^(-1/2)))
      
      # LW
      lw_res = lw(sample_sc)
      scov_lw = (f2^(1/2)) %*% lw_res$covest %*% (f2^(1/2))
      rou_lw = lw_res$coef
      
      # OAS
      oas_res = oas(sample_sc)
      scov_oas = (f2^(1/2)) %*% oas_res$covest %*% (f2^(1/2))
      rou_oas = oas_res$coef  
      
    }
  
    
    # OASD (designed for diagonal target,so don't apply it to corr=TRUE)
    oasd_res = oasd(sample)
    scov_oasd = oasd_res$covest
    rou_oasd = oasd_res$coef
    
    # OASB
    oasb_res = oasb(sample)
    scov_oasb = oasb_res$covest
    gamma_oasb = oasb_res$coef
    
    # SS (designed for diagonal target,so don't apply it to corr=TRUE)
    ss_res = shafer(sample)
    scov_ss = ss_res$covest
    rou_ss = ss_res$coef
    
    # OB
    ob_res = ob(sample,cov)
    scov_ob = ob_res$covest
    gamma_ob = ob_res$coef
    
    # OD 
    od_res = od(sample,cov)
    scov_od = od_res$covest
    rou_od = od_res$coef
    
    
    # Collect errors
    if (!inv_){
      
      err = data.frame(matrix(nrow=1,ncol=8))
      colnames(err) = c('OASD','OAS', 'SS','LW','s','OD','OB', 'OASB')
      err[1,'LW'] = get_err(scov_lw,cov)
      err[1,'SS'] = get_err(scov_ss,cov)
      err[1,'OAS'] = get_err(scov_oas,cov)
      err[1,'OASD'] = get_err(scov_oasd,cov)
      err[1,'OASB'] = get_err(scov_oasb,cov)
      err[1,'s']=sum((scov-cov)^2)
      err[1,'OB'] = get_err(scov_ob,cov)
      err[1,'OD'] = get_err(scov_od,cov)
      
    }else{
      
      err = data.frame(matrix(nrow=1,ncol=8))
      colnames(err) = c('OASD','OAS', 'SS','LW','OD','MP','OB','OASB')
      err[1,'LW'] = get_err(solve(scov_lw),solve(cov))
      err[1,'SS'] = get_err(solve(scov_ss),solve(cov))
      err[1,'OAS'] = get_err(solve(scov_oas),solve(cov))
      err[1,'OASD'] = get_err(solve(scov_oasd),solve(cov))
      err[1,'OASB'] = get_err(solve(scov_oasb),solve(cov))
      err[1,'OB'] = get_err(solve(scov_ob),solve(cov))
      err[1,'MP'] = get_err(ginv(scov),solve(cov))
      err[1,'OD'] = get_err(solve(scov_od),solve(cov))
      
    }
    
    coef_ = data.frame(matrix(nrow=1,ncol=7))
    colnames(coef_) = c('OASD','OAS', 'SS', 'LW','OD','OB','OASB')
    coef_[1,'LW'] = rou_lw
    coef_[1,'SS'] = rou_ss
    coef_[1,'OAS'] = rou_oas
    coef_[1,'OASD'] = rou_oasd
    coef_[1,'OB'] = gamma_ob
    coef_[1,'OD'] = rou_od
    coef_[1,'OASB'] = gamma_oasb
    
    return(list(err,coef_))
  }
  
  res = map(seq(1,iter), sim_err)
  res_err = map_dfr(res, function(x){x[[1]]}) %>% 
    mutate(sd = sd_, r=r_, n=n_)
  res_coef = map_dfr(res, function(x){x[[2]]}) %>% 
    mutate(sd = sd_, r=r_, n=n_)
  
  return(list(res_err, res_coef))
}