my_min_trace <- function(models, method = c("wls_var", "ols", "wls_struct", "mint_cov", "mint_shrink", "mint_shrink_oas", "mint_shrink_oasd", 
                                         "mint_shrink_lw"), sparse = NULL){
  if(is.null(sparse)){
    sparse <- requireNamespace("Matrix", quietly = TRUE)
  }
  # print(cov_est)
  # 1/0
  # stop("error message")
  structure(models, class = c("lst_mint_mdl", "lst_mdl", "list"),
            method = match.arg(method), sparse = sparse)
}


my_forecast.lst_mint_mdl <- function(object, key_data, 
                                  new_data = NULL, h = NULL,
                                  point_forecast = list(.mean = mean), ...){
  method <- object%@%"method"
  sparse <- object%@%"sparse"

  if(sparse){
    require_package("Matrix")
    as.matrix <- Matrix::as.matrix
    t <- Matrix::t
    diag <- function(x) if(is.vector(x)) Matrix::Diagonal(x = x) else Matrix::diag(x)
    solve <- Matrix::solve
    cov2cor <- Matrix::cov2cor
  } else {
    cov2cor <- stats::cov2cor
  }
  
  point_method <- point_forecast
  point_forecast <- list()
  # Get forecasts
  fc <- NextMethod()
  if(length(unique(map(fc, interval))) > 1){
    abort("Reconciliation of temporal hierarchies is not yet supported.")
  }
  
  # Compute weights (sample covariance)
  res <- map(object, function(x, ...) residuals(x, ...), type = "response")
  if(length(unique(map_dbl(res, nrow))) > 1){
    # Join residuals by index #199
    res <- unname(as.matrix(reduce(res, full_join, by = index_var(res[[1]]))[,-1]))
  } else {
    res <- matrix(invoke(c, map(res, `[[`, 2)), ncol = length(object))
  }
  
  # Construct S matrix - ??GA: have moved this here as I need it for Structural scaling
  agg_data <- build_key_data_smat(key_data)
  
  n <- nrow(res)
  covm <- crossprod(stats::na.omit(res)) / n
  if(method == "ols"){
    # OLS
    W <- diag(rep(1L, nrow(covm)))
  } else if(method == "wls_var"){
    # WLS variance scaling
    W <- diag(diag(covm))
  } else if (method == "wls_struct"){
    # WLS structural scaling
    W <- diag(vapply(agg_data$agg,length,integer(1L)))
  } else if (method == "mint_cov"){
    # min_trace covariance
    W <- covm
  } else if (method == "mint_shrink"){
    # min_trace shrink
    tar <- diag(apply(res, 2, compose(crossprod, stats::na.omit))/n)
    corm <- cov2cor(covm)
    xs <- scale(res, center = FALSE, scale = sqrt(diag(covm)))
    xs <- xs[stats::complete.cases(xs),]
    v <- (1/(n * (n - 1))) * (crossprod(xs^2) - 1/n * (crossprod(xs))^2)
    diag(v) <- 0
    corapn <- cov2cor(tar)
    d <- (corm - corapn)^2
    lambda <- sum(v)/sum(d)
    lambda <- max(min(lambda, 1), 0)
    W <- lambda * tar + (1 - lambda) * covm
    
  } else if (method == "mint_shrink_oas"){
    p = nrow(covm)
    f1 = (sum(diag(covm))/p)*diag(p)
    f2 = diag(diag(covm))
    tr_s = sum(diag(covm))
    tr_s2 = sum(diag(t(covm) %*% covm))
    tr_sdiag2 = sum(diag(diag(diag(covm)) %*% diag(diag(covm))))
    
    phi = (tr_s2-(tr_s^2)/p)/(tr_s^2+(1-2/p)*tr_s2)
    rou_oas_constant = min(1/((n_+1-2/p)*phi),1)
    W = (1-rou_oas_constant)*scov+rou_oas_constant*f1
    
  } else if (method == "mint_shrink_oasd"){
    
    p = nrow(covm)
    f1 = (sum(diag(covm))/p)*diag(p)
    f2 = diag(diag(covm))
    tr_s = sum(diag(covm))
    tr_s2 = sum(diag(t(covm) %*% covm))
    tr_sdiag2 = sum(diag(diag(diag(covm)) %*% diag(diag(covm))))
    
    phi1 = (tr_s2-tr_sdiag2)/(tr_s2+tr_s^2-2*tr_sdiag2)
    rou_oas_variant = min(1/((n_+1)*phi1),1)
    W = (1-rou_oas_variant)*covm+rou_oas_variant*f2
  } else if (method == "mint_shrink_lw"){
    p = nrow(covm)
    f1 = (sum(diag(covm))/p)*diag(p)
    f2 = diag(diag(covm))
    tr_s = sum(diag(covm))
    tr_s2 = sum(diag(t(covm) %*% covm))
    tr_sdiag2 = sum(diag(diag(diag(covm)) %*% diag(diag(covm))))
    
    num = 0
    for (i in 1:n_){
      fill = sum((matrix(sample[i,],ncol=1) %*% matrix(sample[i,],nrow=1) - scov)^2)
      num = num + fill
    }
    den = (n_^2)*(tr_s2-(tr_s^2)/p)
    rou_lw = num/den
    rou_lw = min(rou_lw,1)
    W = (1-rou_lw)*scov+rou_lw*f1
  } else{
    abort("Unknown reconciliation method")
  }
  
  # Check positive definiteness of weights
  eigenvalues <- eigen(W, only.values = TRUE)[["values"]]
  if (any(eigenvalues < 1e-8)) {
    abort("min_trace needs covariance matrix to be positive definite.", call. = FALSE)
  }
  
  # Reconciliation matrices
  if(sparse){ 
    row_btm <- agg_data$leaf
    row_agg <- seq_len(nrow(key_data))[-row_btm]
    S <- Matrix::sparseMatrix(
      i = rep(seq_along(agg_data$agg), lengths(agg_data$agg)),
      j = vec_c(!!!agg_data$agg),
      x = rep(1, sum(lengths(agg_data$agg))))
    J <- Matrix::sparseMatrix(i = S[row_btm,,drop = FALSE]@i+1, j = row_btm, x = 1L, 
                              dims = rev(dim(S)))
    U <- cbind(
      Matrix::Diagonal(diff(dim(J))),
      -S[row_agg,,drop = FALSE]
    )
    U <- U[, order(c(row_agg, row_btm)), drop = FALSE]
    Ut <- t(U)
    WUt <- W %*% Ut
    P <- J - J %*% WUt %*% solve(U %*% WUt, U)
    # P <- J - J%*%W%*%t(U)%*%solve(U%*%W%*%t(U))%*%U
  }
  else {
    S <- matrix(0L, nrow = length(agg_data$agg), ncol = max(vec_c(!!!agg_data$agg)))
    S[length(agg_data$agg)*(vec_c(!!!agg_data$agg)-1) + rep(seq_along(agg_data$agg), lengths(agg_data$agg))] <- 1L
    R <- t(S)%*%solve(W)
    P <- solve(R%*%S)%*%R
  }
  
  
  reconcile_fbl_list(fc, S, P, W, point_forecast = point_method)
}

# fabletools.min_trace = my_min_trace2
# fabletools.forecast.lst_mint_mdl = my_forecast.lst_mint_mdl