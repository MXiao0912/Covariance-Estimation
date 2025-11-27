r_new_list = c(0.3,0.5,0.7,0.9)

set.seed(1234)

for (r in r_new_list){
  # True covariance matrix
  corr = matrix(0, ncol=p, nrow=p)
  for (i in 1:p){
    for (j in 1:p){
      corr[i,j]=r^abs(i-j)
    }
  }
  
  sd_s = 1
  sd_b = 100
  sd = diag(c(rep(sd_s, p/2), rep(sd_b, p/2)))
  # sd = diag(rlnorm(p, meanlog = 0, sdlog=1))
  cov = sd %*% corr %*% sd
  
  simulation <- function(n){
    
    sim_err <- function(iter){
      err = data.frame(matrix(nrow=1,ncol=5))
      colnames(err) = c('oas_v','oas_c', 'ss', 'rblw','lw')
      
      # Generate samples from N(0,cov)
      sample = mvrnorm(n, rep(0,p), cov)
      
      # Estimate the covariance matrix using different methods of combining the
      # following two matrices
      scov = (t(sample) %*% sample)/nrow(sample)
      f1 = (sum(diag(scov))/p)*diag(p)
      f2 = diag(diag(scov))
      tr_s = sum(diag(scov))
      tr_s2 = sum(diag(t(scov) %*% scov))
      tr_t = sum(diag(cov))
      tr_t2 = sum(diag(t(cov) %*% cov))
      tr_tdiag2 = sum(diag(diag(diag(cov)) %*% diag(diag(cov))))
      tr_sdiag2 = sum(diag(diag(diag(scov)) %*% diag(diag(scov))))
      
      scorr = cov2cor(scov)
      tr_s_cor = sum(diag(scorr))
      tr_s2_cor = sum(diag(t(scorr) %*% scorr))
      f1_cor = (sum(diag(scorr))/p)*diag(p)
      sample_sc = sample %*% (diag(diag(f2)^(-1/2)))
      
      # Method 2&3: LW & RBLW
      num = 0
      for (i in 1:n){
        fill = sum((matrix(sample_sc[i,],ncol=1) %*% matrix(sample_sc[i,],nrow=1) - scorr)^2)
        num = num + fill
      }
      den = (n^2)*(tr_s2_cor-(tr_s_cor^2)/p)
      rou_lw = num/den
      rou_lw = min(rou_lw,1)
      scorr_lw = (1-rou_lw)*scorr+rou_lw*f1_cor
      scov_lw = (f2^(1/2)) %*% scorr_lw %*% (f2^(1/2))
      
      rou_rblw = (((n-2)/n)*tr_s2_cor+tr_s_cor^2)/((n+2)*(tr_s2_cor-(tr_s_cor^2)/p))
      rou_rblw = min(rou_rblw,1)
      scorr_rblw = (1-rou_rblw)*scorr+rou_rblw*f1_cor
      scov_rblw = (f2^(1/2)) %*% scorr_rblw %*% (f2^(1/2))
      
      # Method 4: oas_c
      phi = (tr_s2_cor-(tr_s_cor^2)/p)/(tr_s_cor^2+(1-2/p)*tr_s2_cor)
      rou_oas_constant = min(1/((n+1-2/p)*phi),1)
      scorr_oas_constant = (1-rou_oas_constant)*scorr+rou_oas_constant*f1_cor
      scov_oas_constant = (f2^(1/2)) %*% scorr_oas_constant %*% (f2^(1/2))
      
      # Method 5: oas_v
      phi1 = (tr_s2-tr_sdiag2)/(tr_s2+tr_s^2-2*tr_sdiag2)
      rou_oas_variant = min(1/((n+1)*phi1),1)
      scov_oas_variant = (1-rou_oas_variant)*scov+rou_oas_variant*f2
      # phi1 = (tr_s2_cor-tr_sdiag2_cor)/(tr_s2_cor+tr_s_cor^2-2*tr_sdiag2_cor)
      # rou_oas_variant = min(1/((n+1)*phi1),1)
      # scorr_oas_variant = (1-rou_oas_variant)*scorr+rou_oas_variant*f1_cor
      # scov_oas_variant = (f2^(1/2)) %*% scorr_oas_variant %*% (f2^(1/2))
      
      # Method 8: Shafer(applied to correlation matrix matrix)
      v <- (1/(n*(n - 1))) * (crossprod(sample_sc^2) - (1/n) * (crossprod(sample_sc))^2)
      diag(v) <- 0
      d <- (scorr - diag(diag(scorr)))^2
      rou_ss <- sum(v)/sum(d)
      rou_ss <- max(min(rou_ss, 1), 0)
      scorr_ss = (1-rou_ss)*scorr+rou_ss*f1_cor
      scov_ss = (f2^(1/2)) %*% scorr_ss %*% (f2^(1/2))
      
      # Collect errors
      err[1,'lw'] = get_err(scov_lw,cov)
      err[1,'rblw'] = get_err(scov_rblw,cov)
      err[1,'ss'] = get_err(scov_ss,cov)
      err[1,'oas_c'] = get_err(scov_oas_constant,cov)
      err[1,'oas_v'] = get_err(scov_oas_variant,cov)
      
      return(err)
    }
    res = map_dfr(seq(1,1000), sim_err) %>% mutate(sample_size = n)
    return(res)
  }
  
  res_sample_size = map_dfr(seq(n_min,n_max),simulation,.progress=TRUE)
  saveRDS(res_sample_size, glue('100r={r}.rds'))
}

