par(mfrow=c(2,2))
r = 0.3
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

eval = eigen(cov)$values
plot(seq_along(eval), eval, xlab = '', ylab = "eigenvalues", main = 'r=0.3')
lines(seq_along(eval),rep(mean(eval),100), type='l',col='red')
  

r = 0.5
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

eval = eigen(cov)$values
plot(seq_along(eval), eval, xlab = '', ylab = "eigenvalues", main = 'r=0.5')
lines(seq_along(eval),rep(mean(eval),100), type='l',col='red')


r = 0.7
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

eval = eigen(cov)$values
plot(seq_along(eval), eval, xlab = '', ylab = "eigenvalues", main = 'r=0.7')
lines(seq_along(eval),rep(mean(eval),100), type='l',col='red')


r = 0.9
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

eval = eigen(cov)$values
plot(seq_along(eval), eval, xlab = '', ylab = "eigenvalues", main = 'r=0.9')
lines(seq_along(eval),rep(mean(eval),100), type='l',col='red')