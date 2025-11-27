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

## Corr Methods
scorr = cov2cor(scov)
corr = cov2cor(cov)
tr_t2_cor = sum(diag(t(corr) %*% corr))
tr_t_cor = sum(diag(corr))
tr_s_cor = sum(diag(scorr))
tr_s2_cor = sum(diag(t(scorr) %*% scorr))
f1_cor = (sum(diag(scorr))/p)*diag(p)

# gamma_2 = ((n^2)/((n-1)*(n+2)))*(1/p)*(tr_s2-(1/n)*tr_s^2)
# gamma_1 = tr_s/p
# test_stat = (n/2)*(gamma_2/(gamma_1^2)-1)

gamma_2 = ((n^2)/((n-1)*(n+2)))*(1/p)*(tr_s2_cor-(1/n)*tr_s_cor^2)
gamma_1 = tr_s_cor/p
test_stat = (n/2)*(gamma_2/(gamma_1^2)-1)

U = (1/p)*sum(diag(((scorr*p/tr_s_cor)-diag(p)) %*% ((scorr*p/tr_s_cor)-diag(p))))
test_stat = n*U-p