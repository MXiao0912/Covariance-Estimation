n=30
p=100
r=0.3

corr = matrix(0, ncol=p, nrow=p)
for (i in 1:p){
  for (j in 1:p){
    corr[i,j]=r^abs(i-j)
  }
}
get_err <- function(cov_est, cov_t){
  error = sum((cov_est-cov_t)^2)/(100^2)
  return(error)
}
sd = diag(rlnorm(p, meanlog = 0, sdlog=0.3))
cov = sd %*% corr %*% sd
tr_t = sum(diag(cov))
tr_t2 = sum(diag(t(cov) %*% cov))
tr_tdiag2 = sum(diag(diag(diag(cov)) %*% diag(diag(cov))))

P = tr_tdiag2
R = tr_t2
Q = tr_t^2

rou_b_o = (Q+R-2*P)/((n+1)*R+Q-(n+2)*P)
gamma_b = 1-(1/rou_b_o)*(((2/n)*P-(2/(n*p))*R)/(((n+2)/n)*P-(2/(n*p))*R-(1/p)*Q))


toiter = function(iter){
  
  sample = mvrnorm(n, rep(0,p), cov)
  scov = (t(sample) %*% sample)/nrow(sample)
  f1 = (sum(diag(scov))/p)*diag(p)
  f2 = diag(diag(scov))
  
  erri = data.frame(matrix(nrow=40,ncol=2))
  colnames(erri) = c('gamma','err')
  
  erri['gamma']=seq(-2,2,length.out=40)
  for (i in 1:40){
    gamma_i=erri[i,'gamma']
    scov_b = (1-rou_b_o)*scov+rou_b_o*(gamma_i*f2+(1-gamma_i)*f1)
    erri[i,'err']= get_err(scov_b, cov)
  }
  erri['iter']=iter
  
  tr_s = sum(diag(scov))
  E=(tr_s/p)*diag(p)
  D=diag(diag(scov))
  a1=sum(diag((E-D) %*% (E-scov)))
  a2=sum(diag((E-D) %*% (E-D)))
  b1=sum(diag((scov-cov) %*% (D-E)))
  b2=sum(diag((D-E) %*% (D-E)))
  
  erri['a1']=a1
  erri['a2']=a2
  erri['b1']=b1
  erri['b2']=b2
  
  return(erri)
}
res = map_dfr(1:1000, toiter, .progress=TRUE)
# mean(unique(res$gamma_b))
# res = res %>% group_by(gamma) %>% summarise(meanerr=mean(err))
# ggplot(res, aes(x=gamma,y=meanerr))+geom_line()
