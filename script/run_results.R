source("script/packages.R")
source("script/DGP.R")
source("script/sim_fct.R")
source("plot.R")

########### baseline ##########################
# varying diagonal size
res <- future_lapply(
    X = seq(sd_min, sd_max, 1),
    FUN = sim_full,
    r_  = (r_min + r_max) / 2,
    n_  = (n_min + n_max) / 2,
    iter = 1000,
    future.seed=TRUE
)
res_error = map_dfr(res, function(x){x[[1]]})
res_coef = map_dfr(res, function(x){x[[2]]})
saveRDS(res_error, 'results/sd.rds')
saveRDS(res_coef, 'results/coef_sd.rds')

# varying sparsity
res <- future_lapply(
  X = seq(r_min, r_max, 0.1),
  FUN = sim_full,
  sd_ = (sd_min + sd_max) / 2,
  n_  = (n_min + n_max) / 2,
  iter = 1000,
  future.seed=TRUE
)


res_error = map_dfr(res, function(x){x[[1]]})
res_coef = map_dfr(res, function(x){x[[2]]})
saveRDS(res_error, 'results/r.rds')
saveRDS(res_coef, 'results/coef_r.rds')

# varying sample size
res <- future_lapply(
  X = seq(n_min, n_max, 1),
  FUN = sim_full,
  sd_ = (sd_min + sd_max) / 2,
  r_  = (r_min + r_max) / 2,
  iter = 1000,
  future.seed=TRUE
)


res_error = map_dfr(res, function(x){x[[1]]})
res_coef = map_dfr(res, function(x){x[[2]]})
saveRDS(res_error, 'results/n.rds')
saveRDS(res_coef, 'results/coef_n.rds')

############## plotting ####################
plot_function("sd")
plot_function("n")
plot_function("r")

############## apply to corr ####################
# varying diagonal size
res <- future_lapply(
  X = seq(sd_min, sd_max, 1),
  FUN = sim_full,
  r_  = (r_min + r_max) / 2,
  n_  = (n_min + n_max) / 2,
  iter = 1000,
  corr_ind=TRUE,
  future.seed=TRUE
)


res_error = map_dfr(res, function(x){x[[1]]})
res_coef = map_dfr(res, function(x){x[[2]]})
saveRDS(res_error, 'results/sd_corr.rds')
saveRDS(res_coef, 'results/coef_sd_corr.rds')

# varying sparsity
res <- future_lapply(
  X = seq(r_min, r_max, 0.1),
  FUN = sim_full,
  sd_ = (sd_min + sd_max) / 2,
  n_  = (n_min + n_max) / 2,
  iter = 1000,
  corr_ind=TRUE,
  future.seed=TRUE
)


res_error = map_dfr(res, function(x){x[[1]]})
res_coef = map_dfr(res, function(x){x[[2]]})
saveRDS(res_error, 'results/r_corr.rds')
saveRDS(res_coef, 'results/coef_r_corr.rds')

# varying sample size
res <- future_lapply(
  X = seq(n_min, n_max, 1),
  FUN = sim_full,
  sd_ = (sd_min + sd_max) / 2,
  r_  = (r_min + r_max) / 2,
  iter = 1000,
  corr_ind=TRUE,
  future.seed=TRUE
)

res_error = map_dfr(res, function(x){x[[1]]})
res_coef = map_dfr(res, function(x){x[[2]]})
saveRDS(res_error, 'results/n_corr.rds')
saveRDS(res_coef, 'results/coef_n_corr.rds')

############## plotting ####################
# doesn't display OD
plot_agg_function("corr")

############## apply to inv ####################
# varying diagonal size
res <- future_lapply(
  X = seq(sd_min, sd_max, 1),
  FUN = sim_full,
  r_  = (r_min + r_max) / 2,
  n_  = (n_min + n_max) / 2,
  iter = 1000,
  inv=TRUE,
  future.seed=TRUE
)

res_error = map_dfr(res, function(x){x[[1]]})
res_coef = map_dfr(res, function(x){x[[2]]})
saveRDS(res_error, 'results/sd_inv.rds')
saveRDS(res_coef, 'results/coef_sd_inv.rds')

# varying sparsity
res <- future_lapply(
  X = seq(r_min, r_max, 0.1),
  FUN = sim_full,
  sd_ = (sd_min + sd_max) / 2,
  n_  = (n_min + n_max) / 2,
  iter = 1000,
  inv=TRUE,
  future.seed=TRUE
)


res_error = map_dfr(res, function(x){x[[1]]})
res_coef = map_dfr(res, function(x){x[[2]]})
saveRDS(res_error, 'results/r_inv.rds')
saveRDS(res_coef, 'results/coef_r_inv.rds')

# varying sample size
res <- future_lapply(
  X = seq(n_min, n_max, 1),
  FUN = sim_full,
  sd_ = (sd_min + sd_max) / 2,
  r_  = (r_min + r_max) / 2,
  iter = 1000,
  corr_ind=TRUE,
  inv=TRUE,
  future.seed=TRUE
)


res_error = map_dfr(res, function(x){x[[1]]})
res_coef = map_dfr(res, function(x){x[[2]]})
saveRDS(res_error, 'results/n_inv.rds')
saveRDS(res_coef, 'results/coef_n_inv.rds')

############## plotting ####################
# doesn't display OD
plot_agg_function("inv")


########### t_distibution ##########################
# varying diagonal size
res <- future_lapply(
  X = seq(sd_min, sd_max, 1),
  FUN = sim_full,
  r_  = (r_min + r_max) / 2,
  n_  = (n_min + n_max) / 2,
  iter = 1000,
  dist="t",
  future.seed=TRUE
)


res_error = map_dfr(res, function(x){x[[1]]})
res_coef = map_dfr(res, function(x){x[[2]]})
saveRDS(res_error, 'results/t/sd.rds')
saveRDS(res_coef, 'results/t/coef_sd.rds')

# varying sparsity
res <- future_lapply(
  X = seq(r_min, r_max, 0.1),
  FUN = sim_full,
  sd_ = (sd_min + sd_max) / 2,
  n_  = (n_min + n_max) / 2,
  iter = 1000,
  dist="t",
  future.seed=TRUE
)


res_error = map_dfr(res, function(x){x[[1]]})
res_coef = map_dfr(res, function(x){x[[2]]})
saveRDS(res_error, 'results/t/r.rds')
saveRDS(res_coef, 'results/t/coef_r.rds')

# varying sample size
res <- future_lapply(
  X = seq(n_min, n_max, 1),
  FUN = sim_full,
  sd_ = (sd_min + sd_max) / 2,
  r_  = (r_min + r_max) / 2,
  iter = 1000,
  corr_ind=TRUE,
  dist="t",
  future.seed=TRUE
)


res_error = map_dfr(res, function(x){x[[1]]})
res_coef = map_dfr(res, function(x){x[[2]]})
saveRDS(res_error, 'results/t/n.rds')
saveRDS(res_coef, 'results/t/coef_n.rds')

plot_function("sd",destination='results/t/')
plot_function("n",destination='results/t/')
plot_function("r",destination='results/t/')


########### uniform distribution ##########################
# varying diagonal size
res <- future_lapply(
  X = seq(sd_min, sd_max, 1),
  FUN = sim_full,
  r_  = (r_min + r_max) / 2,
  n_  = (n_min + n_max) / 2,
  iter = 1000,
  dist="uniform",
  future.seed=TRUE
)

res_error = map_dfr(res, function(x){x[[1]]})
res_coef = map_dfr(res, function(x){x[[2]]})
saveRDS(res_error, 'results/uniform/sd.rds')
saveRDS(res_coef, 'results/uniform/coef_sd.rds')

# varying sparsity
res <- future_lapply(
  X = seq(r_min, r_max, 0.1),
  FUN = sim_full,
  sd_ = (sd_min + sd_max) / 2,
  n_  = (n_min + n_max) / 2,
  iter = 1000,
  dist="uniform",
  future.seed=TRUE
)


res_error = map_dfr(res, function(x){x[[1]]})
res_coef = map_dfr(res, function(x){x[[2]]})
saveRDS(res_error, 'results/uniform/r.rds')
saveRDS(res_coef, 'results/uniform/coef_r.rds')

# varying sample size
res <- future_lapply(
  X = seq(n_min, n_max, 1),
  FUN = sim_full,
  sd_ = (sd_min + sd_max) / 2,
  r_  = (r_min + r_max) / 2,
  iter = 1000,
  corr_ind=TRUE,
  dist="uniform",
  future.seed=TRUE
)


res_error = map_dfr(res, function(x){x[[1]]})
res_coef = map_dfr(res, function(x){x[[2]]})
saveRDS(res_error, 'results/uniform/n.rds')
saveRDS(res_coef, 'results/uniform/coef_n.rds')

plot_function("sd",destination='results/uniform/')
plot_function("n",destination='results/uniform/')
plot_function("r",destination='results/uniform/')

