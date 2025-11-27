library(Matrix)
library(latex2exp)

p = 1000
n = 500
r = 0.3
sdlog=2

corr = matrix(0, ncol=p, nrow=p)
for (i in 1:p){
  for (j in 1:p){
    corr[i,j]=r^abs(i-j)
  }
}
sd = diag(rlnorm(p, meanlog = 0, sdlog=sdlog))
cov = sd %*% corr %*% sd
# rhos_ls = vector(mode='list', length=1000)
rhooas_ls = vector(mode='list', length=1000)
for (iter in 1:1000){
  print(iter)
  sample = mvrnorm(n, rep(0,p), cov)
  S = (t(sample) %*% sample)/n
  tr_s = sum(diag(S))
  tr_s2 = sum(diag(t(S) %*% S))
  tr_sdiag2 = sum(diag(diag(diag(S)) %*% diag(diag(S))))
  # rhos = (tr_s2-2*tr_sdiag2+tr_s^2)/((n+1)*tr_s2+tr_s^2-(n+2)*tr_sdiag2)  
  # rhos_ls[[iter]] = rhos
  phi = (tr_s2-tr_sdiag2)/(tr_s2+tr_s^2-2*tr_sdiag2)
  rhooas = min(1/((n+1)*phi),1)
  rhooas_ls[[iter]] = rhooas
}

tr_t = sum(diag(cov))
tr_t2 = sum(diag(t(cov) %*% cov))
tr_tdiag2 = sum(diag(diag(diag(cov)) %*% diag(diag(cov))))
rho_oracle = (tr_t2-2*tr_tdiag2+tr_t^2)/((n+1)*tr_t2+tr_t^2-(n+2)*tr_tdiag2) 
plot(unlist(rhooas_ls), ylim=c(0,1), main=glue('n={n}, r={r},sdlog={sdlog},OAS_diag'))
abline(h=rho_oracle, col='red')
abline(h=mean(unlist(rhooas_ls)),col='blue')

rhos_ls = vector(mode='list', length=1000)
rhooas_ls = vector(mode='list', length=1000)
rhooasc_ls = vector(mode='list', length=1000)
rhosc_ls = vector(mode='list', length=1000)
for (iter in 1:1000){
  sample = mvrnorm(n, rep(0,100), cov)
  S = (t(sample) %*% sample)/n
  tr_s = sum(diag(S))
  tr_s2 = sum(diag(t(S) %*% S))
  tr_sdiag2 = sum(diag(diag(diag(S)) %*% diag(diag(S))))
  rhos = (tr_s2-2*tr_sdiag2+tr_s^2)/((n+1)*tr_s2+tr_s^2-(n+2)*tr_sdiag2)  
  phi = (tr_s2-tr_sdiag2)/(tr_s2+tr_s^2-2*tr_sdiag2)
  rhooas = min(1/((n+1)*phi),1)
  rhosc = ((1-2/p)*tr_s2 + tr_s^2)/((n+1-2/p)*tr_s2+(1-n/p)*tr_s^2)
  phi1 = (tr_s2-(tr_s^2)/p)/(tr_s^2+(1-2/p)*tr_s2)
  rhooasc = min(1/((n+1-2/p)*phi1),1)
  rhos_ls[[iter]] = rhos
  rhooas_ls[[iter]] = rhooas
  rhooasc_ls[[iter]] = rhooasc
  rhosc_ls[[iter]] = rhosc
}

