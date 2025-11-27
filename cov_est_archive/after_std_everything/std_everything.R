library(MASS)
library(purrr)
library(dplyr)
library(glue)
library(tidyverse)
r_new_list = c(0.3,0.5,0.7,0.9)
p = 100
n_min = 6
n_max = 30

set.seed(1234)
shrink.estim <- function(res)
{
  n<-nrow(res)
  covm <- cov(res)
  tar <- diag(diag(covm))
  corm <- stats::cov2cor(covm)
  xs <- scale(res, center = FALSE, scale = sqrt(diag(covm)))
  xs <- xs[stats::complete.cases(xs),]
  v <- (n/(n - 1)^3) * (crossprod(xs^2) - 1/(n-1) * (crossprod(xs))^2) # corrected
  diag(v) <- 0
  corapn <- stats::cov2cor(tar)
  d <- (corm - corapn)^2
  lambda <- sum(v)/sum(d)
  lambda <- max(min(lambda, 1), 0)
  W <- lambda * tar + (1 - lambda) * covm
  return(W)
}

for (r in r_new_list){
  # True covariance matrix
  corr = matrix(0, ncol=p, nrow=p)
  for (i in 1:p){
    for (j in 1:p){
      corr[i,j]=r^abs(i-j)
    }
  }
  
  sd_s = 1
  sd_b = 10
  sd = diag(c(rep(sd_s, p/2), rep(sd_b, p/2)))
  cov = sd %*% corr %*% sd
  
  simulation <- function(n){
    
    sim_err <- function(iter){
      err = data.frame(matrix(nrow=1,ncol=11))
      colnames(err) = c('corr_oracle_c','corr_oracle_v','corr_oas_c', 'corr_oas_v','oracle_c',
                        'oracle_v','lw','rblw','oas_c','oas_v','shafer')
      
      # Generate samples from N(0,cov)
      sample = mvrnorm(n, rep(0,p), cov)
      
      # Estimate the covariance matrix using different methods of combining the
      # following two matrices
      scov = (t(sample) %*% sample)/nrow(sample)
      # scov = cov(sample)
      f1 = (sum(diag(scov))/p)*diag(p)
      f2 = diag(diag(scov))
      tr_s = sum(diag(scov))
      tr_s2 = sum(diag(t(scov) %*% scov))
      tr_t = sum(diag(cov))
      tr_t2 = sum(diag(t(cov) %*% cov))
      tr_tdiag2 = sum(diag(diag(diag(cov)) %*% diag(diag(cov))))
      tr_sdiag2 = sum(diag(diag(diag(scov)) %*% diag(diag(scov))))
      
      ## Corr Methods
      scorr = cov2cor(scov)
      corr = cov2cor(cov)
      tr_t2_cor = sum(diag(t(corr) %*% corr))
      tr_t_cor = sum(diag(corr))
      tr_s_cor = sum(diag(scorr))
      tr_s2_cor = sum(diag(t(scorr) %*% scorr))
      tr_tdiag2_cor = sum(diag(diag(diag(corr)) %*% diag(diag(corr))))
      tr_sdiag2_cor = sum(diag(diag(diag(scorr)) %*% diag(diag(scorr))))
      f1_cor = (sum(diag(scorr))/p)*diag(p)
      f2_cor = diag(diag(scorr))
      
      ## Method1: Benchmark (oracle with constant and variant variance)
      rou_o_constant = ((1-2/p)*tr_t2 + tr_t^2)/((n+1-2/p)*tr_t2+(1-n/p)*tr_t^2)
      scov_o_constant = (1-rou_o_constant)*scov+rou_o_constant*f1
      rou_o_variant = ((1/n)*tr_t2-(2/n)*tr_tdiag2+(1/n)*tr_t^2)/(((n+1)/n)*tr_t2+(1/n)*tr_t^2-((n+2)/n)*tr_tdiag2)
      scov_o_variant = (1-rou_o_variant)*scov+rou_o_variant*f2
      
      ## Method2&3: LW & RBLW (trying to make the oracle under constant variance estimable)
      num = 0
      for (i in 1:n){
        fill = sum((matrix(sample[i,],ncol=1) %*% matrix(sample[i,],nrow=1) - scov)^2)
        num = num + fill
      }
      den = (n^2)*(tr_s2-(tr_s^2)/p)
      rou_lw = num/den
      rou_lw = min(rou_lw,1)
      scov_lw = (1-rou_lw)*scov+rou_lw*f1
      
      rou_rblw = (((n-2)/n)*tr_s2+tr_s^2)/((n+2)*(tr_s2-(tr_s^2)/p))
      rou_rblw = min(rou_rblw,1)
      scov_rblw = (1-rou_rblw)*scov+rou_rblw*f1
      
      ## Method4: OAS with constant variance
      phi = (tr_s2-(tr_s^2)/p)/(tr_s^2+(1-2/p)*tr_s2)
      rou_oas_constant = min(1/((n+1-2/p)*phi),1)
      scov_oas_constant = (1-rou_oas_constant)*scov+rou_oas_constant*f1
      
      ## Method5: OAS with variant variance
      phi1 = (tr_s2-tr_sdiag2)/(tr_s2+tr_s^2-2*tr_sdiag2)
      rou_oas_variant = min(1/((n+1)*phi1),1)
      scov_oas_variant = (1-rou_oas_variant)*scov+rou_oas_variant*f2
      
      # Method6: Oracle_c with constant variance applied to correlation estimation
      rou_o_constant_cor = ((1-2/p)*tr_t2_cor + tr_t_cor^2)/((n+1-2/p)*tr_t2_cor+(1-n/p)*tr_t_cor^2)
      scorr_o_constant_cor = (1-rou_o_constant_cor)*scorr+rou_o_constant_cor*f1_cor
      scov_o_constant_cor = diag(sqrt(diag(scov))) %*% scorr_o_constant_cor %*% diag(sqrt(diag(scov)))
      
      # Method7: OAS with constant variance applied to correlation estimation
      phi_cor = (tr_s2_cor-(tr_s_cor^2)/p)/(tr_s_cor^2+(1-2/p)*tr_s2_cor)
      rou_oas_constant_cor = min(1/((n+1-2/p)*phi_cor),1)
      scorr_oas_constant = (1-rou_oas_constant_cor)*scorr+rou_oas_constant_cor*f1_cor
      scov_oas_constant_cor = diag(sqrt(diag(scov))) %*% scorr_oas_constant %*% diag(sqrt(diag(scov)))
      
      # Oracle_v applied to correlation estimation
      rou_v_cor = (tr_t2_cor-2*tr_tdiag2_cor+tr_t_cor^2)/((n+1)*tr_t2_cor+tr_t_cor^2-(n+2)*tr_tdiag2_cor)
      scorr_v_cor = (1-rou_v_cor)*scorr+rou_v_cor*f2_cor
      scov_v_cor = diag(sqrt(diag(scov))) %*% scorr_v_cor %*% diag(sqrt(diag(scov)))
      
      # OAS_v applied to correlation estimation
      phi_cor_v = (tr_s2_cor-tr_sdiag2_cor)/(tr_s2_cor+tr_s_cor^2-2*tr_sdiag2_cor)
      rou_oas_v_cor = min(1/((n+1)*phi_cor_v),1)
      scorr_oas_v_cor = (1-rou_oas_v_cor)*scorr+rou_oas_v_cor*f2_cor
      scov_oas_v_cor = diag(sqrt(diag(scov))) %*% scorr_oas_v_cor %*% diag(sqrt(diag(scov)))
      
      # Method8: Shafer
      scov_shafer = shrink.estim(sample)
      
      
      # Collect errors
      err[1,'corr_oracle_c'] = get_err(scov_o_constant_cor,cov)
      err[1,'corr_oracle_v'] = get_err(scov_v_cor,cov)
      err[1,'corr_oas_c'] = get_err(scov_oas_constant_cor,cov)
      err[1,'corr_oas_v'] = get_err(scov_oas_v_cor,cov)
      
      err[1,'oracle_c'] = get_err(scov_o_constant,cov)
      err[1,'oracle_v'] = get_err(scov_o_variant,cov)
      err[1,'lw'] = get_err(scov_lw,cov)
      err[1,'rblw'] = get_err(scov_rblw,cov)
      err[1,'oas_c'] = get_err(scov_oas_constant,cov)
      err[1,'oas_v'] = get_err(scov_oas_variant,cov)
      err[1, 'shafer'] = get_err(scov_shafer, cov)
      
      return(err)
    }
    res = map_dfr(seq(1,1000), sim_err) %>% mutate(sample_size = n)
    return(res)
  }
  
  res_sample_size = map_dfr(seq(n_min,n_max),simulation,.progress=TRUE)
  saveRDS(res_sample_size, glue('r={r}.rds'))
}

