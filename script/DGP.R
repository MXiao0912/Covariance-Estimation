p = 100
n_min = 6
n_max = 30
r_min = 0
r_max = 0.9
sd_min = 1
sd_max = 20
small = 1
iter = 1000

get_err <- function(e, t){
  err_e = sum((e-t)^2)
  return(err_e)
}