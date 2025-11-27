# Packages
library(MASS)
library(purrr)
library(dplyr)
library(glue)
library(tidyverse)

# Create covariance matrix
p = 100
n_min = 6
n_max = 30
r_list = c(0, 0.1, 0.5, 0.9)
r_new_list = c(0.05,0.2,0.3,0.4,0.6,0.7,0.8)


# corr = matrix(0, ncol=p, nrow=p)
# for (i in 1:p){
#   for (j in 1:p){
#     corr[i,j]=r^abs(i-j)
#   }
# }
# 
# sd_s = 1
# sd_b = sqrt(10)
# sd = diag(c(rep(sd_s, p/2), rep(sd_b, p/2)))
# cov = sd %*% corr %*% sd
# 
# # Generate samples from N(0,cov)
# sample = mvrnorm(n, rep(0,p), cov)
# 
# # Estimate the covariance matrix using different methods of combining the
# # following two matrices
# scov = cov(sample)
# f1 = (sum(diag(scov))/p)*diag(p)
# f2 = diag(diag(scov))
# tr_s = sum(diag(scov))
# tr_s2 = sum(diag(t(scov) %*% scov))
# tr_t = sum(diag(cov))
# tr_t2 = sum(diag(t(cov) %*% cov))
# tr_tdiag2 = sum(diag(diag(diag(cov)) %*% diag(diag(cov))))
# tr_sdiag2 = sum(diag(diag(diag(scov)) %*% diag(diag(scov))))
# 
# ## Method1: Benchmark (oracle with constant and variant variance)
# rou_o_constant = ((1-2/p)*tr_t2 + tr_t^2)/((n+1-2/p)*tr_t2+(1-n/p)*tr_t^2)
# scov_o_constant = (1-rou_o_constant)*scov+rou_o_constant*f1
# rou_o_variant = ((1/n)*tr_t2-(2/n)*tr_tdiag2+(1/n)*tr_t^2)/(((n+1)/n)*tr_t2+(1/n)*tr_t^2-((n+2)/n)*tr_tdiag2)
# scov_o_variant = (1-rou_o_variant)*scov+rou_o_variant*f2
# 
# ## Method2&3: LW & RBLW (trying to make the oracle under constant variance estimable)
# num = 0
# for (i in 1:n){
#   fill = sum((t(sample[i,]) %*% sample[i,] - scov)^2)
#   num = num + fill
# }
# den = n^2(tr_s2-(tr_s^2)/p)
# rou_lw = num/den
# rou_lw = min(rou_lw,1)
# scov_lw = (1-rou_lw)*scov+rou_lw*f1
# 
# rou_rblw = (((n-2)/n)*tr_s2+tr_s^2)/((n+2)*(tr_s2-(tr_s^2)/p))
# rou_rblw = min(rou_rblw,1)
# scov_rblw = (1-rou_rblw)*scov+rou_rblw*f1
# 
# ## Method4: OAS with constant variance
# phi = (tr_s2-(tr_s^2)/p)/(tr_s2+(1-2/p)*tr_s^2)
# rou_oas_constant = min(1/((n+1-2/p)*phi),1)
# scov_oas_constant = (1-rou_oas_constant)*scov+rou_oas_constant*f1
# 
# ## Method5: OAS with variant variance
# phi1 = (tr_s2-tr_sdiag2)/(tr_s2+tr_s^2-2*tr_sdiag2)
# rou_oas_variant = min(1/((n+1)*phi1),1)
# scov_oas_variant = (1-rou_oas_variant)*scov+rou_oas_variant*f2

get_err <- function(cov_est, cov_t){
  error = sum((cov_est-cov_t)^2)/(100^2)
  return(error)
}

# set.seed(1234)
# 
# for (r in r_new_list){
#   # True covariance matrix
#   corr = matrix(0, ncol=p, nrow=p)
#   for (i in 1:p){
#     for (j in 1:p){
#       corr[i,j]=r^abs(i-j)
#     }
#   }
#   
#   sd_s = 1
#   sd_b = sqrt(10)
#   sd = diag(c(rep(sd_s, p/2), rep(sd_b, p/2)))
#   cov = sd %*% corr %*% sd
#   
#   simulation <- function(n){
#     
#     sim_err <- function(iter){
#       err = data.frame(matrix(nrow=1,ncol=6))
#       colnames(err) = c('oracle_c','oracle_v','lw','rblw','oas_c','oas_v')
#       
#       # Generate samples from N(0,cov)
#       sample = mvrnorm(n, rep(0,p), cov)
#       
#       # Estimate the covariance matrix using different methods of combining the
#       # following two matrices
#       scov = cov(sample)
#       f1 = (sum(diag(scov))/p)*diag(p)
#       f2 = diag(diag(scov))
#       tr_s = sum(diag(scov))
#       tr_s2 = sum(diag(t(scov) %*% scov))
#       tr_t = sum(diag(cov))
#       tr_t2 = sum(diag(t(cov) %*% cov))
#       tr_tdiag2 = sum(diag(diag(diag(cov)) %*% diag(diag(cov))))
#       tr_sdiag2 = sum(diag(diag(diag(scov)) %*% diag(diag(scov))))
#       
#       ## Method1: Benchmark (oracle with constant and variant variance)
#       rou_o_constant = ((1-2/p)*tr_t2 + tr_t^2)/((n+1-2/p)*tr_t2+(1-n/p)*tr_t^2)
#       scov_o_constant = (1-rou_o_constant)*scov+rou_o_constant*f1
#       rou_o_variant = ((1/n)*tr_t2-(2/n)*tr_tdiag2+(1/n)*tr_t^2)/(((n+1)/n)*tr_t2+(1/n)*tr_t^2-((n+2)/n)*tr_tdiag2)
#       scov_o_variant = (1-rou_o_variant)*scov+rou_o_variant*f2
#       
#       ## Method2&3: LW & RBLW (trying to make the oracle under constant variance estimable)
#       num = 0
#       for (i in 1:n){
#         fill = sum((matrix(sample[i,],ncol=1) %*% matrix(sample[i,],nrow=1) - scov)^2)
#         num = num + fill
#       }
#       den = (n^2)*(tr_s2-(tr_s^2)/p)
#       rou_lw = num/den
#       rou_lw = min(rou_lw,1)
#       scov_lw = (1-rou_lw)*scov+rou_lw*f1
#       
#       rou_rblw = (((n-2)/n)*tr_s2+tr_s^2)/((n+2)*(tr_s2-(tr_s^2)/p))
#       rou_rblw = min(rou_rblw,1)
#       scov_rblw = (1-rou_rblw)*scov+rou_rblw*f1
#       
#       ## Method4: OAS with constant variance
#       phi = (tr_s2-(tr_s^2)/p)/(tr_s^2+(1-2/p)*tr_s2)
#       rou_oas_constant = min(1/((n+1-2/p)*phi),1)
#       scov_oas_constant = (1-rou_oas_constant)*scov+rou_oas_constant*f1
#       
#       ## Method5: OAS with variant variance
#       phi1 = (tr_s2-tr_sdiag2)/(tr_s2+tr_s^2-2*tr_sdiag2)
#       rou_oas_variant = min(1/((n+1)*phi1),1)
#       scov_oas_variant = (1-rou_oas_variant)*scov+rou_oas_variant*f2
#       
#       # Collect errors
#       err[1,'oracle_c'] = get_err(scov_o_constant,cov)
#       err[1,'oracle_v'] = get_err(scov_o_variant,cov)
#       err[1,'lw'] = get_err(scov_lw,cov)
#       err[1,'rblw'] = get_err(scov_rblw,cov)
#       err[1,'oas_c'] = get_err(scov_oas_constant,cov)
#       err[1,'oas_v'] = get_err(scov_oas_variant,cov)
#       
#       return(err)
#     }
#     res = map_dfr(seq(1,1000), sim_err) %>% mutate(sample_size = n)
#     return(res)
#   }
#   
#   res_sample_size = map_dfr(seq(n_min,n_max),simulation,.progress=TRUE)
#   saveRDS(res_sample_size, glue('r={r}.rds'))
# }
# 
# # Evaluate results
for (r in r_new_list[1]){
  data = readRDS(glue('r={r}.rds')) %>% 
    group_by(sample_size) %>% 
    summarise(across(everything(), function(x){(x^2)/length(x)})) %>%
    pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
  ggplot(data, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
    labs(title = glue('MSE for r={r}'))
  ggsave(glue('r={r}.png'))
}

for (r in r_new_list){
  df = readRDS(glue('r={r}.rds')) %>% 
    group_by(sample_size) %>% 
    summarise(across(everything(), mean)) %>%
    select(sample_size, test, oas_v, corr_oas, oracle_v, shafer) %>% 
    pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
  ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
    labs(title = glue('MSE for r={r}'))
  ggsave(glue('r={r}.png'))
}


for (r in r_new_list){
  df = readRDS(glue('r={r}_rou.rds')) %>% 
    group_by(sample_size) %>% 
    summarise(across(everything(), mean)) %>%
    pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'rou')
  ggplot(df, aes(x= sample_size, y=rou, group=methods, color=methods))+geom_line()+
    labs(title = glue('average rou for r={r}'))
  ggsave(glue('r={r}_rou.png'))
}

library(gridExtra)
library(ggpubr)
df = readRDS(glue('r=0.9_offdiag.rds')) %>% 
  group_by(sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
p4 = ggplot(df, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
  labs(title = glue('MSE for r=0.9'))

ggarrange(p1, p2, p3, p4, nrow=2, ncol = 2, common.legend = TRUE, legend="bottom")

# 

# 
# # Examine error contributions from diagonal
# p = 100
# n_min = 6
# n_max = 30
# get_err_diag <- function(cov_est, cov_t){
#   error = sum(diag((cov_est-cov_t)^2))/(100)
#   return(error)
# }


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

shrink.estim_coef <- function(res)
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
  return(lambda)
}

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
  sd_b = 10
  sd = diag(c(rep(sd_s, p/2), rep(sd_b, p/2)))
  cov = sd %*% corr %*% sd
  
  simulation <- function(n){
    
    sim_err <- function(iter){
      err = data.frame(matrix(nrow=1,ncol=9))
      colnames(err) = c('oracle_c','oracle_v','lw','rblw','oas_c','oas_v','corr_oracle','corr_oas', 'shafer')
      
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
      
      ## Corr Methods
      scorr = cov2cor(scov)
      corr = cov2cor(cov)
      tr_t2_cor = sum(diag(t(corr) %*% corr))
      tr_t_cor = sum(diag(corr))
      tr_s_cor = sum(diag(scorr))
      tr_s2_cor = sum(diag(t(scorr) %*% scorr))
      f1_cor = (sum(diag(scorr))/p)*diag(p)
      
      # Method6: Oracle with constant variance applied to correlation estimation
      rou_o_constant_cor = ((1-2/p)*tr_t2_cor + tr_t_cor^2)/((n+1-2/p)*tr_t2_cor+(1-n/p)*tr_t_cor^2)
      scorr_o_constant_cor = (1-rou_o_constant_cor)*scorr+rou_o_constant_cor*f1_cor
      scov_o_constant_cor = diag(sqrt(diag(scov))) %*% scorr_o_constant_cor %*% diag(sqrt(diag(scov)))
      
      # Method7: OAS with constant variance applied to correlation estimation
      phi_cor = (tr_s2_cor-(tr_s_cor^2)/p)/(tr_s_cor^2+(1-2/p)*tr_s2_cor)
      rou_oas_constant_cor = min(1/((n+1-2/p)*phi_cor),1)
      scorr_oas_constant = (1-rou_oas_constant_cor)*scorr+rou_oas_constant_cor*f1_cor
      scov_oas_constant_cor = diag(sqrt(diag(scov))) %*% scorr_oas_constant %*% diag(sqrt(diag(scov)))
      
      # Method8: Shafer
      # scov_shafer = shrink.estim(sample)
      v <- (1/(n*(n - 1))) * (crossprod(sample^2) - (1/n) * (crossprod(sample))^2)
      diag(v) <- 0
      d <- (scov - f2)^2
      rou_ss <- sum(v)/sum(d)
      rou_ss <- max(min(rou_ss, 1), 0)
      scov_shafer = (1-rou_ss)*scov+rou_ss*f2
      
      # # Method9: test for diagonality before choosing between OAS or OAS_diag
      # gamma_2 = ((n^2)/((n-1)*(n+2)))*(1/p)*(tr_s2_cor-(1/n)*tr_s_cor^2)
      # gamma_1 = tr_s_cor/p
      # test_stat = (n/2)*((gamma_2/(gamma_1^2))-1)
      # if (test_stat>2){
      #   scov_test = scov_oas_variant
      # }else{
      #   scov_test = scov_oas_constant_cor
      # }
      
      # Collect errors
      err[1,'oracle_c'] = get_err(scov_o_constant,cov)
      err[1,'oracle_v'] = get_err(scov_o_variant,cov)
      err[1,'lw'] = get_err(scov_lw,cov)
      err[1,'rblw'] = get_err(scov_rblw,cov)
      err[1,'oas_c'] = get_err(scov_oas_constant,cov)
      err[1,'oas_v'] = get_err(scov_oas_variant,cov)
      err[1,'corr_oracle'] = get_err(scov_o_constant_cor,cov)
      err[1,'corr_oas'] = get_err(scov_oas_constant_cor,cov)
      err[1, 'shafer'] = get_err(scov_shafer, cov)
      # err[1, 'test'] = get_err(scov_test, cov)
      
      return(err)
    }
    res = map_dfr(seq(1,1000), sim_err) %>% mutate(sample_size = n)
    return(res)
  }
  
  res_sample_size = map_dfr(seq(n_min,n_max),simulation,.progress=TRUE)
  saveRDS(res_sample_size, glue('testr={r}.rds'))
}


# added both:
set.seed(1234)
p = 100
n_min = 6
n_max = 30
r_new_list = c(0.3,0.5,0.7,0.9)
for (r in r_new_list){
  # True covariance matrix
  corr = matrix(0, ncol=p, nrow=p)
  for (i in 1:p){
    for (j in 1:p){
      corr[i,j]=r^abs(i-j)
    }
  }
  
  # sd_s = 1
  # sd_b = 10
  # sd = diag(c(rep(sd_s, p/2), rep(sd_b, p/2)))
  # sd = diag(rep(1,p))
  sd = diag(rlnorm(p, meanlog = 0, sdlog=0.1))
  cov = sd %*% corr %*% sd
  
  simulation <- function(n){
    
    sim_err <- function(iter){
      # rous = data.frame(matrix(nrow=1,ncol=9))
      # err = data.frame(matrix(nrow=1,ncol=8))
      err = data.frame(matrix(nrow=1,ncol=2))
      # colnames(err) = c('oracle_v','oas_v','oracle_c','oas_c', 'shafer', 'rblw','lw','oracle_b')
      colnames(err) = c('oracle_b','oas_b')
      
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
      
      # ## Method1: Benchmark (oracle with constant and variant variance)
      # rou_o_constant = ((1-2/p)*tr_t2 + tr_t^2)/((n+1-2/p)*tr_t2+(1-n/p)*tr_t^2)
      # scov_o_constant = (1-rou_o_constant)*scov+rou_o_constant*f1
      # rou_o_variant = ((1/n)*tr_t2-(2/n)*tr_tdiag2+(1/n)*tr_t^2)/(((n+1)/n)*tr_t2+(1/n)*tr_t^2-((n+2)/n)*tr_tdiag2)
      # scov_o_variant = (1-rou_o_variant)*scov+rou_o_variant*f2
      # # 
      # # ## Method2&3: LW & RBLW (trying to make the oracle under constant variance estimable)
      # num = 0
      # for (i in 1:n){
      #   fill = sum((matrix(sample[i,],ncol=1) %*% matrix(sample[i,],nrow=1) - scov)^2)
      #   num = num + fill
      # }
      # den = (n^2)*(tr_s2-(tr_s^2)/p)
      # rou_lw = num/den
      # rou_lw = min(rou_lw,1)
      # scov_lw = (1-rou_lw)*scov+rou_lw*f1
      # 
      # rou_rblw = (((n-2)/n)*tr_s2+tr_s^2)/((n+2)*(tr_s2-(tr_s^2)/p))
      # rou_rblw = min(rou_rblw,1)
      # scov_rblw = (1-rou_rblw)*scov+rou_rblw*f1
      # # 
      # # ## Method4: OAS with constant variance
      # phi = (tr_s2-(tr_s^2)/p)/(tr_s^2+(1-2/p)*tr_s2)
      # rou_oas_constant = min(1/((n+1-2/p)*phi),1)
      # scov_oas_constant = (1-rou_oas_constant)*scov+rou_oas_constant*f1
      # 
      # ## Method5: OAS with variant variance
      # phi1 = (tr_s2-tr_sdiag2)/(tr_s2+tr_s^2-2*tr_sdiag2)
      # rou_oas_variant = min(1/((n+1)*phi1),1)
      # scov_oas_variant = (1-rou_oas_variant)*scov+rou_oas_variant*f2
      
      ## Corr Methods
      # scorr = cov2cor(scov)
      # corr = cov2cor(cov)
      # tr_t2_cor = sum(diag(t(corr) %*% corr))
      # tr_t_cor = sum(diag(corr))
      # tr_s_cor = sum(diag(scorr))
      # tr_s2_cor = sum(diag(t(scorr) %*% scorr))
      # f1_cor = (sum(diag(scorr))/p)*diag(p)
      
      # Method6: Oracle with constant variance applied to correlation estimation
      # rou_o_constant_cor = ((1-2/p)*tr_t2_cor + tr_t_cor^2)/((n+1-2/p)*tr_t2_cor+(1-n/p)*tr_t_cor^2)
      # scorr_o_constant_cor = (1-rou_o_constant_cor)*scorr+rou_o_constant_cor*f1_cor
      # scov_o_constant_cor = diag(sqrt(diag(scov))) %*% scorr_o_constant_cor %*% diag(sqrt(diag(scov)))
      # 
      # Method7: OAS with constant variance applied to correlation estimation
      # phi_cor = (tr_s2_cor-(tr_s_cor^2)/p)/(tr_s_cor^2+(1-2/p)*tr_s2_cor)
      # rou_oas_constant_cor = min(1/((n+1-2/p)*phi_cor),1)
      # scorr_oas_constant = (1-rou_oas_constant_cor)*scorr+rou_oas_constant_cor*f1_cor
      # scov_oas_constant_cor = diag(sqrt(diag(scov))) %*% scorr_oas_constant %*% diag(sqrt(diag(scov)))
      # 
      # Method8: Shafer
      # rou_shafer = shrink.estim_coef(sample)
      # scov_shafer = (1-rou_shafer)*scov+rou_shafer*f2
      
      # Method9: test for diagonality before choosing between OAS or OAS_diag
      # gamma_2 = ((n^2)/((n-1)*(n+2)))*(1/p)*(tr_s2_cor-(1/n)*tr_s_cor^2)
      # gamma_1 = tr_s_cor/p
      # test_stat = (n/2)*((gamma_2/(gamma_1^2))-1)
      # if (test_stat>2){
      #   scov_test = scov_oas_variant
      # }else{
      #   scov_test = scov_oas_constant_cor
      # }
      
      # Method10: Different shrinkage for diag and offdiag
      # A = (tr_s^2)/p
      # B = tr_s2-tr_sdiag2
      # C = tr_sdiag2-A
      
      P = tr_tdiag2
      R = tr_t2
      Q = tr_t^2
      
      # if (A>max((p/((1+n)*(p-1)))*C+(2/(p-1))*B, (n/(p-1))*B)){
      #   rou_b = 1
      #   gamma_b = 0
      # } else if ((A>(C-B)/(p-1) & A<((n+2)/(n*(p-1)))*B-(p/(n*(p-1)))*C)|(A>((n+2)/(n*(p-1)))*B+(p/(n*(p-1)))*C & A<(n/(p-1))*B)){
      #   rou_b = ((p-1)*A+B)/((1+n)*B)
      #   gamma_b = ((p-1)*A-n*B)/((p-1)*A+B)
      # }
      
      rou_b_o = (Q+R-2*P)/((n+1)*R+Q-(n+2)*P)
      gamma_b_o = 1-(1/rou_b_o)*(((2/n)*P-(2/(n*p))*R)/(((n+2)/n)*P-(2/(n*p))*R-(1/p)*Q))
      # gamma_b_o = min(max(gamma_b_o, 0),1)
      # print(gamma_b_o)
      
      scov_both_o = (1-rou_b_o)*scov+rou_b_o*(gamma_b_o*f2+(1-gamma_b_o)*f1)
      # scov_both = (1-rou_b)*scov+rou_b*(gamma_b*f2+(1-gamma_b)*f1)
      
      A = (tr_s^2)/p
      B = tr_s2-tr_sdiag2
      C = tr_sdiag2-A
      
      t1 = ((p-1)*C)/((p-1)*(A+C)-B)
      t2 = (B-(p-1)*C)/((p-1)*(A+C)-B)
      t3 = (n*p*C)/(2*(p-1)*(A+C)-2*B)
      
      u2 = C/(B-C+(p-1)*A)
      u = B/(B-C+(p-1)*A)
      u1 = u-u2
      
      
      if (A>(n/(p-1))*max(c((p/2)*C, B))){
        rou_b = 1
        gamma_b = 0
      } else if((A>(n*B/(p-1))) & (A<= ((n*p)/(2*(p-1)))*C)){
        rou_b = 1
        gamma_b = (t3-t2-1)/(t1+t3)
      } else if(((A>(((n+1)*n*p*C)/(2*(n+2)*(p-1))+(n*B/((n+2)*(p-1))))) & (A<=n*B/(p-1)))|(A<(n*B/((n+2)*(p-1)))-((n+1)*n*p*C)/(2*(n+2)*(p-1)))){
        rou_b = (1+u2)/((n+1)*u)
        gamma_b = 1-((n+1)*u)/(1+u2)
      } else{
        rou_b = (t1+t3+u2)/((u1+n*u)*(t1+t3)+(t3-t2)*u2)
        gamma_b = (t3-t2-u1-n*u)/(t1+t3+u2)
      }
      scov_both_oas = (1-rou_b)*scov+rou_b*(gamma_b*f2+(1-gamma_b)*f1)
      
      
      # Collect errors
      # err[1,'oracle_c'] = get_err(scov_o_constant,cov)
      # err[1,'oracle_v'] = get_err(scov_o_variant,cov)
      # # rous[1,'oracle_v'] = rou_o_variant
      # err[1,'lw'] = get_err(scov_lw,cov)
      # err[1,'rblw'] = get_err(scov_rblw,cov)
      # err[1,'shafer'] = get_err(scov_shafer,cov)
      # err[1,'oas_c'] = get_err(scov_oas_constant,cov)
      # err[1,'oas_v'] = get_err(scov_oas_variant,cov)
      err[1,'oracle_b'] = get_err(scov_both_o,cov)
      err[1,'oas_b'] = get_err(scov_both_oas,cov)
      # err[1,'oas_b'] = get_err(scov_both,cov)
      # rous[1,'oas_v'] = rou_oas_variant
      # rous[1,'corr_oracle'] = rou_o_constant_cor
      # rous[1,'corr_oas'] = rou_oas_constant_cor
      # rous[1, 'shafer'] = rou_shafer
      # err[1, 'test'] = get_err(scov_test, cov)
      
      return(err)
    }
    res = map_dfr(seq(1,1000), sim_err) %>% mutate(sample_size = n)
    return(res)
  }
  
  res_sample_size = map_dfr(seq(n_min,n_max),simulation,.progress=TRUE)
  saveRDS(res_sample_size, glue('r={r}.rds'))
}



# # Examine off-diagonal errors
# get_err_offdiag <- function(cov_est, cov_t){
#   error_diag = sum(diag((cov_est-cov_t)^2))
#   error_off_diag = (sum((cov_est-cov_t)^2)-error_diag)/(100^2-100)
#   return(error_off_diag)
# }
# 
# set.seed(1234)
# 
# for (r in r_new_list){
#   # True covariance matrix
#   corr = matrix(0, ncol=p, nrow=p)
#   for (i in 1:p){
#     for (j in 1:p){
#       corr[i,j]=r^abs(i-j)
#     }
#   }
#   
#   sd_s = 1
#   sd_b = sqrt(10)
#   sd = diag(c(rep(sd_s, p/2), rep(sd_b, p/2)))
#   cov = sd %*% corr %*% sd
#   
#   simulation <- function(n){
#     
#     sim_err <- function(iter){
#       err = data.frame(matrix(nrow=1,ncol=6))
#       colnames(err) = c('oracle_c','oracle_v','lw','rblw','oas_c','oas_v')
#       
#       # Generate samples from N(0,cov)
#       sample = mvrnorm(n, rep(0,p), cov)
#       
#       # Estimate the covariance matrix using different methods of combining the
#       # following two matrices
#       scov = cov(sample)
#       f1 = (sum(diag(scov))/p)*diag(p)
#       f2 = diag(diag(scov))
#       tr_s = sum(diag(scov))
#       tr_s2 = sum(diag(t(scov) %*% scov))
#       tr_t = sum(diag(cov))
#       tr_t2 = sum(diag(t(cov) %*% cov))
#       tr_tdiag2 = sum(diag(diag(diag(cov)) %*% diag(diag(cov))))
#       tr_sdiag2 = sum(diag(diag(diag(scov)) %*% diag(diag(scov))))
#       
#       ## Method1: Benchmark (oracle with constant and variant variance)
#       rou_o_constant = ((1-2/p)*tr_t2 + tr_t^2)/((n+1-2/p)*tr_t2+(1-n/p)*tr_t^2)
#       scov_o_constant = (1-rou_o_constant)*scov+rou_o_constant*f1
#       rou_o_variant = ((1/n)*tr_t2-(2/n)*tr_tdiag2+(1/n)*tr_t^2)/(((n+1)/n)*tr_t2+(1/n)*tr_t^2-((n+2)/n)*tr_tdiag2)
#       scov_o_variant = (1-rou_o_variant)*scov+rou_o_variant*f2
#       
#       ## Method2&3: LW & RBLW (trying to make the oracle under constant variance estimable)
#       num = 0
#       for (i in 1:n){
#         fill = sum((matrix(sample[i,],ncol=1) %*% matrix(sample[i,],nrow=1) - scov)^2)
#         num = num + fill
#       }
#       den = (n^2)*(tr_s2-(tr_s^2)/p)
#       rou_lw = num/den
#       rou_lw = min(rou_lw,1)
#       scov_lw = (1-rou_lw)*scov+rou_lw*f1
#       
#       rou_rblw = (((n-2)/n)*tr_s2+tr_s^2)/((n+2)*(tr_s2-(tr_s^2)/p))
#       rou_rblw = min(rou_rblw,1)
#       scov_rblw = (1-rou_rblw)*scov+rou_rblw*f1
#       
#       ## Method4: OAS with constant variance
#       phi = (tr_s2-(tr_s^2)/p)/(tr_s^2+(1-2/p)*tr_s2)
#       rou_oas_constant = min(1/((n+1-2/p)*phi),1)
#       scov_oas_constant = (1-rou_oas_constant)*scov+rou_oas_constant*f1
#       
#       ## Method5: OAS with variant variance
#       phi1 = (tr_s2-tr_sdiag2)/(tr_s2+tr_s^2-2*tr_sdiag2)
#       rou_oas_variant = min(1/((n+1)*phi1),1)
#       scov_oas_variant = (1-rou_oas_variant)*scov+rou_oas_variant*f2
#       
#       # Collect errors
#       err[1,'oracle_c'] = get_err_offdiag(scov_o_constant,cov)
#       err[1,'oracle_v'] = get_err_offdiag(scov_o_variant,cov)
#       err[1,'lw'] = get_err_offdiag(scov_lw,cov)
#       err[1,'rblw'] = get_err_offdiag(scov_rblw,cov)
#       err[1,'oas_c'] = get_err_offdiag(scov_oas_constant,cov)
#       err[1,'oas_v'] = get_err_offdiag(scov_oas_variant,cov)
#       
#       return(err)
#     }
#     res = map_dfr(seq(1,1000), sim_err) %>% mutate(sample_size = n)
#     return(res)
#   }
#   
#   res_sample_size = map_dfr(seq(n_min,n_max),simulation,.progress=TRUE)
#   saveRDS(res_sample_size, glue('r={r}_offdiag.rds'))
# }
# 
# # Evaluate results
# for (r in r_new_list){
#   data = readRDS(glue('r={r}_diag.rds'))%>% group_by(sample_size) %>% summarise(across(everything(),mean)) %>% 
#     pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
#   ggplot(data, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
#     labs(title = glue('Diagonal MSE for r={r}'))
#   ggsave(glue('r={r}_diag.png'))
# }
# 
# for (r in r_new_list){
#   data = readRDS(glue('r={r}_offdiag.rds'))%>% group_by(sample_size) %>% summarise(across(everything(),mean)) %>% 
#     pivot_longer(cols = -sample_size, names_to = 'methods','values_to' = 'MSE')
#   ggplot(data, aes(x= sample_size, y=MSE, group=methods, color=methods))+geom_line()+
#     labs(title = glue('Off-diagonal MSE for r={r}'))
#   ggsave(glue('r={r}_offdiag.png'))
# }