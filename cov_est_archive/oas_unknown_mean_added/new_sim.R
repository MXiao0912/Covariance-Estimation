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
      err = data.frame(matrix(nrow=1,ncol=2))
      colnames(err) = c('oas_v_unknown_mean','shafer')
      
      # Generate samples from N(0,cov)
      sample = mvrnorm(n, rep(0,p), cov)
      
      # Estimate the covariance matrix using different methods of combining the
      # following two matrices
      # scov = (t(sample) %*% sample)/nrow(sample)
      scov = cov(sample)
      f1 = (sum(diag(scov))/p)*diag(p)
      f2 = diag(diag(scov))
      tr_s = sum(diag(scov))
      tr_s2 = sum(diag(t(scov) %*% scov))
      tr_t = sum(diag(cov))
      tr_t2 = sum(diag(t(cov) %*% cov))
      tr_tdiag2 = sum(diag(diag(diag(cov)) %*% diag(diag(cov))))
      tr_sdiag2 = sum(diag(diag(diag(scov)) %*% diag(diag(scov))))
      
      
      # OAS_v but with unknown mean
      phi_unknown = (tr_s2-tr_sdiag2)/((-2)*tr_sdiag2+tr_s2+tr_s^2)
      rou_unknown = min(1/(n*phi_unknown),1)
      scov_oas_unknown_mean = (1-rou_unknown)*scov+rou_unknown*f2
      
     
      # Method8: Shafer
      scov_shafer = shrink.estim(sample)
      
      
      # Collect errors
      err[1,'oas_v_unknown_mean'] = get_err(scov_oas_unknown_mean,cov)
      err[1, 'shafer'] = get_err(scov_shafer, cov)
      
      return(err)
    }
    res = map_dfr(seq(1,1000), sim_err) %>% mutate(sample_size = n)
    return(res)
  }
  
  res_sample_size = map_dfr(seq(n_min,n_max),simulation,.progress=TRUE)
  saveRDS(res_sample_size, glue('r={r}.rds'))
}
