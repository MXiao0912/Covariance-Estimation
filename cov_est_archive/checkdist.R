# check derivation of the distribution of sum of squared sample variances

library(Matrix)
p = 300
n = 200
cov = diag(rep(1,p))

ans_ls = vector(mode='list', length=1000)
for (iter in 1:1000){
  print(iter)
  sample = mvrnorm(n, rep(0,p), cov)
  S = (t(sample) %*% sample)/n
  tr_sdiag2 = sum(diag(diag(diag(S)) %*% diag(diag(S))))
  ans_ls[iter] = tr_sdiag2
}