# check what happens to the eigenvalues and eigenvectors when the true covariance
# matrix is the identity matrix
library(Matrix)
library(pracma)

evalue = diag(rep(1,100))
evec = randortho(100)
cov = evec %*% evalue %*% t(evec)

n = 50
sample = mvrnorm(n, rep(0,100), cov)
S = (t(sample) %*% sample)/n

evals = eigen(S)$values
evecs = eigen(S)$vectors
diags = diag(diag(S))
evalsd = eigen(diags)$values
evecsd = eigen(diags)$vectors

tr_s = sum(diag(S))
tr_s2 = sum(diag(t(S) %*% S))
tr_sdiag2 = sum(diag(diag(diag(S)) %*% diag(diag(S))))
rhos = (tr_s2-2*tr_sdiag2+tr_s^2)/((n+1)*tr_s2+tr_s^2-(n+2)*tr_sdiag2)  
phi = (tr_s2-tr_sdiag2)/(tr_s2+tr_s^2-2*tr_sdiag2)
rhooas = min(1/((n+1)*phi),1)

soas = (1-rhooas)*S+rhooas*diag(diag(S))
evals_oas = eigen(soas)$values
evec_oas = eigen(soas)$vectors