# simulate the evolution of shrinkage coefficient

set.seed(1234)
p = 100
n = 30
r = 0.5
corr = matrix(0, ncol=p, nrow=p)
for (i in 1:p){
  for (j in 1:p){
    corr[i,j]=r^abs(i-j)
  }
}

sd = diag(rlnorm(p, meanlog = 0, sdlog=0.5))
cov = sd %*% corr %*% sd

sim_err <- function(iter){

    err = data.frame(matrix(nrow=100,ncol=1))
    
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
    
    # ## Method1: Benchmark (oracle with constant and variant variance)
    rou_o_variant = ((1/n)*tr_t2-(2/n)*tr_tdiag2+(1/n)*tr_t^2)/(((n+1)/n)*tr_t2+(1/n)*tr_t^2-((n+2)/n)*tr_tdiag2)

    # ## Method5: OAS with variant variance
    phi1 = (tr_s2-tr_sdiag2)/(tr_s2+tr_s^2-2*tr_sdiag2)
    rou_oas_variant = min(1/((n+1)*phi1),1)
    err[100,1] = (rou_o_variant-rou_oas_variant)^2
    
    rou_init = ((1/n)*tr_s2-(2/n)*tr_sdiag2+(1/n)*tr_s^2)/(((n+1)/n)*tr_s2+(1/n)*tr_s^2-((n+2)/n)*tr_sdiag2)
    err[1,1] = (rou_init-rou_o_variant)^2
    rou_t = rou_init
    for (it in 2:99){
      sigma_t = (1-rou_t)*scov+rou_t*f2
      rou_t = ((1/n)*sum(diag(sigma_t %*% scov))-(2/n)*tr_sdiag2+(1/n)*tr_s^2)/(((n+1)/n)*sum(diag(sigma_t %*% scov))+(1/n)*tr_s^2-((n+2)/n)*tr_sdiag2)
      err[it,1] = (rou_t-rou_o_variant)^2
    }
    return(err)
  }
res = map_dfc(seq(1,1000), sim_err,.progress=TRUE)
res_agg = rowSums(res)

# simulate for the separate shrinkage case
set.seed(1234)
p = 100
n = 30
r = 0.5
corr = matrix(0, ncol=p, nrow=p)
for (i in 1:p){
  for (j in 1:p){
    corr[i,j]=r^abs(i-j)
  }
}

sd = diag(rlnorm(p, meanlog = 0, sdlog=0.5))
cov = sd %*% corr %*% sd

sim_err <- function(iter){
  
  err_rou = data.frame(matrix(nrow=30,ncol=1))
  err_gamma = data.frame(matrix(nrow=30,ncol=1))
  
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
  
  P = tr_tdiag2
  R = tr_t2
  Q = tr_t^2
  
  P1 = tr_sdiag2
  R1 = tr_s2
  Q1 = tr_s^2
  
  # ## Method1: Benchmark (oracle with constant and variant variance)
  rou_b_o = (Q+R-2*P)/((n+1)*R+Q-(n+2)*P)
  gamma_b_o = 1-(1/rou_b_o)*(((2/n)*P-(2/(n*p))*R)/(((n+2)/n)*P-(2/(n*p))*R-(1/p)*Q))
  
  # ## Method5: OAS with variant variance
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
  
  err_rou[30,1] = (rou_b_o-rou_b)^2
  err_gamma[30,1] = (gamma_b_o-gamma_b)^2
  
  rou_init = (Q1+R1-2*P1)/((n+1)*R1+Q1-(n+2)*P1)
  err_rou[1,1] = (rou_init-rou_b_o)^2
  rou_t = rou_init
  for (it in 2:29){
    sigma_t = (1-rou_t)*scov+rou_t*(gamma_b*f2+(1-gamma_b)*f1)
    rou_t = (sum(diag(sigma_t))^2+sum(diag(sigma_t %*% scov))-2*sum(diag(diag(diag(sigma_t)) %*% diag(diag(scov)))))/((n+1)*sum(diag(sigma_t %*% scov))+sum(diag(sigma_t))^2-(n+2)*sum(diag(diag(diag(sigma_t)) %*% diag(diag(scov)))))
    err_rou[it,1] = (rou_t-rou_b_o)^2
  }

  gamma_init = 1-(1/rou_b)*(((2/n)*P1-(2/(n*p))*R1)/(((n+2)/n)*P1-(2/(n*p))*R1-(1/p)*Q1))
  err_gamma[1,1] = (gamma_init-gamma_b_o)^2
  gamma_t = gamma_init
  for (it in 2:29){
    sigma_t = (1-rou_b)*scov+rou_b*(gamma_t*f2+(1-gamma_t)*f1)
    gamma_t = 1-(1/rou_b)*(2*p*sum(diag(diag(diag(sigma_t)) %*% diag(diag(scov))))-2*sum(diag(sigma_t %*% scov)))/(p*(n+2)*sum(diag(diag(diag(sigma_t)) %*% diag(diag(scov))))-2*sum(diag(sigma_t %*% scov))-n*sum(diag(sigma_t))^2)
    err_gamma[it,1] = (gamma_t-gamma_b_o)^2
  }
  # gamma_init = 1-(1/rou_init)*(((2/n)*P1-(2/(n*p))*R1)/(((n+2)/n)*P1-(2/(n*p))*R1-(1/p)*Q1))
  # err_gamma[1,1] = (gamma_init-gamma_b_o)^2
  # gamma_t = gamma_init
  # for (it in 2:29){
  #   sigma_t = (1-rou_t)*scov+rou_t*(gamma_t*f2+(1-gamma_t)*f1)
  #   rou_t = (sum(diag(sigma_t))^2+sum(diag(sigma_t %*% scov))-2*sum(diag(diag(diag(sigma_t)) %*% diag(diag(scov)))))/((n+1)*sum(diag(sigma_t %*% scov))+sum(diag(sigma_t))^2-(n+2)*sum(diag(diag(diag(sigma_t)) %*% diag(diag(scov)))))
  #   gamma_t = 1-(1/rou_t)*(2*p*sum(diag(diag(diag(sigma_t)) %*% diag(diag(scov))))-2*sum(diag(sigma_t %*% scov)))/(p*(n+2)*sum(diag(diag(diag(sigma_t)) %*% diag(diag(scov))))-2*sum(diag(sigma_t %*% scov))-n*sum(diag(sigma_t))^2)
  #   err_gamma[it,1] = (gamma_t-gamma_b_o)^2
  #   err_rou[it,1] = (rou_t-rou_b_o)^2
  # }
  
  return(list(err_rou, err_gamma))
}
res = map(seq(1,1000), sim_err,.progress=TRUE)

res_rou = map_dfc(res, function(n){n[[1]]})
res_gamma = map_dfc(res, function(n){n[[2]]})
res_rou = rowSums(res_rou)
res_gamma = rowSums(res_gamma)