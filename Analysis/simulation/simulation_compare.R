rm(list = ls())
gc()

# use default library paths


options(repos = c(CRAN = "https://cloud.r-project.org"))

cat("Current .libPaths():\n")
print(.libPaths())

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(stringr)
  library(tidyr)
  library(qvalue)
})
# =========================================================
# 0. user settings / arguments
# =========================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript run_three_sample_alpha_final.R <var_alpha> [n_sim]")
}

var_alpha <- as.numeric(args[1])
n_sim <- ifelse(length(args) >= 2, as.integer(args[2]), 300)

if (is.na(var_alpha)) stop("var_alpha must be numeric.")
if (is.na(n_sim) || n_sim <= 0) stop("n_sim must be a positive integer.")

cat("Running simulation with var_alpha =", var_alpha, "and n_sim =", n_sim, "\n")

# =========================================================
# 1. paths
# =========================================================
work_dir <- getwd()
out_dir  <- file.path(work_dir, "three_sample_alpha")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

setwd(work_dir)

# =========================================================
# 2. helper functions
# =========================================================
fast_marginal_pvals <- function(X, y) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- nrow(X)
  
  X_center <- sweep(X, 2, colMeans(X), FUN = "-")
  y_center <- y - mean(y)
  
  Sxx <- colSums(X_center^2)
  Sxx[Sxx <= 0] <- NA_real_
  
  Sxy <- colSums(X_center * y_center)
  
  beta_hat <- Sxy / Sxx
  RSS <- sum(y_center^2) - (Sxy^2) / Sxx
  RSS[RSS < 0] <- 0
  
  sigma2_hat <- RSS / (n - 2)
  se_beta <- sqrt(sigma2_hat / Sxx)
  t_stat <- beta_hat / se_beta
  
  pval <- 2 * pt(-abs(t_stat), df = n - 2)
  pval[is.na(pval)] <- 1
  pval
}

fast_cross_pvals <- function(X, Y) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- nrow(X)
  
  Xc <- sweep(X, 2, colMeans(X), FUN = "-")
  Yc <- sweep(Y, 2, colMeans(Y), FUN = "-")
  
  Sxx <- colSums(Xc^2)
  Syy <- colSums(Yc^2)
  Sxx[Sxx <= 0] <- NA_real_
  Syy[Syy <= 0] <- NA_real_
  
  Sxy <- crossprod(Xc, Yc)
  
  beta_hat <- sweep(Sxy, 1, Sxx, "/")
  
  RSS <- matrix(rep(Syy, each = length(Sxx)), nrow = length(Sxx)) -
    sweep(Sxy^2, 1, Sxx, "/")
  RSS[RSS < 0] <- 0
  
  sigma2 <- RSS / (n - 2)
  se_beta <- sqrt(sweep(sigma2, 1, Sxx, "/"))
  
  t_stat <- beta_hat / se_beta
  pval <- 2 * pt(-abs(t_stat), df = n - 2)
  pval[is.na(pval)] <- 1
  pval
}

nonnullPropEst <- function(x, u, sigma) {
  z  <- (x - u) / sigma
  xi <- c(0:100) / 100
  tmax <- sqrt(log(length(x)))
  tt <- seq(0, tmax, 0.1)
  
  epsest <- NULL
  for (j in seq_along(tt)) {
    t <- tt[j]
    f <- t * xi
    f <- exp(f^2 / 2)
    w <- (1 - abs(xi))
    co <- 0 * xi
    for (i in 1:101) {
      co[i] <- mean(cos(t * xi[i] * z))
    }
    epshat <- 1 - sum(w * f * co) / sum(w)
    epsest <- c(epsest, epshat)
  }
  max(epsest)
}

normalize <- function(v) {
  if (sum(v != 0) > 0) return(v / sqrt(sum(v^2)))
  rep(0, length(v))
}

force_sp <- function(df1) {
  a1 <- dim(df1)[2]
  df2 <- df1
  pos <- NULL
  for (i in 1:a1) {
    df2[pos, i] <- 0
    pos <- c(which(df2[, i] != 0))
  }
  df2
}

soft <- function(a, para) {
  b <- sort(abs(a))
  b <- abs(a) - para
  b <- (b + abs(b)) / 2
  sign(a) * b
}

estimate_uv <- function(W, pen.u, pen.v, init_u, init_v,
                        prec = 1e-03, max.iter = 200, verbose = FALSE) {
  dim.u <- dim(W)[1]
  dim.v <- dim(W)[2]
  if (is.null(init_u)) init_u <- normalize(rnorm(dim.u))
  if (is.null(init_v)) init_v <- normalize(rnorm(dim.v))
  
  u0 <- init_u
  v0 <- init_v
  diffs <- 2
  iter <- 1
  
  while ((diffs > prec) && (iter < max.iter)) {
    u1 <- as.numeric(normalize(soft(W %*% v0, para = pen.u)))
    v1 <- as.numeric(normalize(soft(t(W) %*% u0, para = pen.v)))
    q  <- as.numeric(t(u1) %*% W %*% v1)
    
    if (sum(is.nan(u1)) > 0 || sum(is.nan(v1)) > 0) {
      return(list(u = rep(0, dim.u), v = rep(0, dim.v), q = 0))
    }
    
    diffs <- sum((u0 - u1)^2) / dim.u + sum((v0 - v1)^2) / dim.v
    u0 <- u1
    v0 <- v1
    iter <- iter + 1
  }
  
  list(u = u0, v = v0, q = q)
}

estimate_uv_sp <- function(W, pen.u, pen.v, pos.u = NULL, pos.v = NULL,
                           init_u, init_v, prec = 1e-03, max.iter = 200, verbose = FALSE) {
  dim.u <- dim(W)[1]
  dim.v <- dim(W)[2]
  if (is.null(init_u)) init_u <- normalize(rnorm(dim.u))
  if (is.null(init_v)) init_v <- normalize(rnorm(dim.v))
  
  u0 <- init_u
  v0 <- init_v
  diffs <- 2
  iter <- 1
  
  while ((diffs > prec) && (iter < max.iter)) {
    u1 <- as.numeric(soft(W %*% v0, para = pen.u))
    if (!is.null(pos.u)) u1[pos.u] <- 0
    u1 <- normalize(u1)
    
    v1 <- as.numeric(soft(t(W) %*% u0, para = pen.v))
    if (!is.null(pos.v)) v1[pos.v] <- 0
    v1 <- normalize(v1)
    
    q <- as.numeric(t(u1) %*% W %*% v1)
    diffs <- sum((u0 - u1)^2) / dim.u + sum((v0 - v1)^2) / dim.v
    
    u0 <- u1
    v0 <- v1
    iter <- iter + 1
  }
  
  list(u = u0, v = v0, q = q)
}

intersect_min_12 <- function(W, u10, v10, u20, v20, dec.u = 0, dec.v = 0,
                             prec.grid.u = 100, prec.grid.v = 100,
                             max.iter = 200, method = "both", verbose = FALSE) {
  grid.u <- seq(0, 1, by = 1 / prec.grid.u)
  grid.v <- seq(0, 1, by = 1 / prec.grid.v)
  
  max.iter <- min(prec.grid.v, prec.grid.u) - 1
  iter <- 1
  conv_u <- length(u10)
  conv_v <- length(v10)
  move.u <- 1
  move.v <- 1
  
  if (method == "both") {
    while ((conv_u > dec.u || conv_v > dec.v) && (iter < max.iter)) {
      q <- 0
      if (move.u == 1) p.u <- grid.u[iter]
      if (move.v == 1) p.v <- grid.v[iter]
      
      W_n <- W - q * (u10 %*% t(v10))
      A1 <- estimate_uv_sp(W_n, pen.u = p.u, pen.v = p.v, init_u = u10, init_v = v10, verbose = FALSE)
      u1 <- A1$u; v1 <- A1$v; q1 <- A1$q
      
      W_n <- W - q1 * (u1 %*% t(v1))
      A2 <- estimate_uv_sp(W_n, pen.u = p.u, pen.v = p.v, init_u = u20, init_v = v20, verbose = FALSE)
      u2 <- A2$u; v2 <- A2$v
      
      conv_u <- length(intersect(which(abs(u1) > 0), which(abs(u2) > 0)))
      conv_v <- length(intersect(which(abs(v1) > 0), which(abs(v2) > 0)))
      move.u <- as.numeric(conv_u > dec.u)
      move.v <- as.numeric(conv_v > dec.v)
      iter <- iter + 1
    }
  }
  
  list(comp1 = A1, comp2 = A2)
}

intersect_min_3on <- function(W, u, v, u30, v30, dec.u = 0, dec.v = 0,
                              prec.grid.u = 100, prec.grid.v = 100,
                              method = "both", verbose = FALSE, max.iter = 200) {
  grid.u <- seq(0, 1, by = 1 / prec.grid.u)
  grid.v <- seq(0, 1, by = 1 / prec.grid.v)
  iter <- 1
  conv_u <- length(u30)
  conv_v <- length(v30)
  move.u <- 1
  move.v <- 1
  max.iter <- min(prec.grid.v, prec.grid.u) - 1
  W_n <- W
  
  if (method == "both") {
    while ((conv_u > dec.u || conv_v > dec.v) && (iter < max.iter)) {
      q <- 0
      if (move.u == 1) p.u <- grid.u[iter]
      if (move.v == 1) p.v <- grid.v[iter]
      
      u1 <- u
      v1 <- v
      
      A2 <- estimate_uv(W_n, pen.u = p.u, pen.v = p.v, init_u = u30, init_v = v30, verbose = FALSE)
      u2 <- A2$u
      v2 <- A2$v
      
      conv_u <- length(intersect(which(abs(u1) > 0), which(abs(u2) > 0)))
      conv_v <- length(intersect(which(abs(v1) > 0), which(abs(v2) > 0)))
      move.u <- as.numeric(conv_u > dec.u)
      move.v <- as.numeric(conv_v > dec.v)
      iter <- iter + 1
      W_n <- W - q * (u1 %*% t(v1))
    }
  }
  
  list(comp = A2)
}

archie_work <- function(Sigma_GE, Sigma_GG, Sigma_EE, K = NULL,
                        dec.u = 0, dec.v = 0,
                        prec.grid.u = 100, prec.grid.v = 100,
                        method = "both", verbose = FALSE) {
  W <- Sigma_GG %*% Sigma_GE %*% Sigma_EE
  
  svd.g <- eigen(W %*% t(W))
  if (is.null(K)) {
    d1 <- -diff(svd.g$values)
    K <- which.max(d1)
  }
  
  u10 <- normalize(svd.g$vectors[, 1])
  v10 <- normalize(t(W) %*% svd.g$vectors[, 1])
  u20 <- normalize(svd.g$vectors[, 2])
  v20 <- normalize(t(W) %*% svd.g$vectors[, 2])
  
  obj <- intersect_min_12(W, u10, v10, u20, v20, method = "both",
                          prec.grid.u = 100, prec.grid.v = 100,
                          max.iter = 1000, verbose = verbose)
  
  df.u <- cbind(obj$comp1$u, obj$comp2$u)
  df.v <- cbind(obj$comp1$v, obj$comp2$v)
  qs   <- c(obj$comp1$q, obj$comp2$q)
  
  es <- eigen(W %*% t(W))
  ev <- es$values
  uvecs <- es$vectors
  uvecs <- apply(uvecs, 2, normalize)
  vvecs <- t(W) %*% uvecs
  vvecs <- apply(vvecs, 2, normalize)
  
  W_n <- W - ev[1] * (uvecs[, 1] %*% t(vvecs[, 1])) - ev[2] * (uvecs[, 2] %*% t(vvecs[, 2]))
  qs[1] <- qs[1]^2 / sum(svd.g$values^2)
  W_2 <- W - ev[1] * (uvecs[, 1] %*% t(vvecs[, 1]))
  svd.2 <- eigen(W_2 %*% t(W_2))
  qs[2] <- qs[2]^2 / sum(svd.2$values^2)
  
  if (K > 2) {
    for (i in 3:K) {
      init.u <- normalize(uvecs[, i])
      init.v <- normalize(vvecs[, i])
      obj.1 <- intersect_min_3on(
        W_n, u30 = init.u, v30 = init.v,
        u = df.u[, (i - 1)], v = df.v[, (i - 1)],
        dec.u = 0, dec.v = 0, method = "both",
        prec.grid.u = 100, prec.grid.v = 100, verbose = verbose
      )
      df.u <- cbind(df.u, obj.1$comp$u)
      df.u <- force_sp(df.u)
      df.v <- cbind(df.v, obj.1$comp$v)
      df.v <- force_sp(df.v)
      W_n <- W_n - ev[i] * df.u[, i] %*% t(df.v[, i])
      svd.n <- eigen(W_n %*% t(W_n))
      qs <- c(qs, obj.1$comp$q / sum(svd.n$values^2))
    }
  }
  
  df.u.sel <- df.u
  df.u.sel[df.u.sel != 0] <- 1
  df.v.sel <- df.v
  df.v.sel[df.v.sel != 0] <- 1
  
  list(us = df.u.sel, vs = df.v.sel, qs = qs, K = K)
}

res <- function(p_alpha, p_beta, FDR_list, true_set, num_trans) {
  result <- NULL
  num_set <- nrow(p_alpha)
  select.idx <- vector("list", length(FDR_list))
  
  for (i in 1:num_set) {
    p_a <- p_alpha[i, ]
    p_b <- p_beta
    
    p_a[p_a == 0] <- min(p_a[p_a != 0])
    p_b[p_b == 0] <- min(p_b[p_b != 0])
    p_a[p_a == 1] <- max(p_a[p_a != 1])
    p_b[p_b == 1] <- max(p_b[p_b != 1])
    
    Z_a <- stats::qnorm(p_a, lower.tail = FALSE)
    Z_b <- stats::qnorm(p_b, lower.tail = FALSE)
    
    pi0a <- 1 - nonnullPropEst(Z_a, 0, 1)
    pi0b <- 1 - nonnullPropEst(Z_b, 0, 1)
    
    if (!is.na(pi0a) && !is.na(pi0b) && pi0a >= 0 && pi0b >= 0) {
      if (pi0a > 1) pi0a <- 1
      if (pi0b > 1) pi0b <- 1
      
      p.mat <- cbind(p_a, p_b)
      p3 <- (apply(p.mat, 1, max))^2
      wg1 <- pi0a * (1 - pi0b)
      wg2 <- (1 - pi0a) * pi0b
      wg3 <- pi0a * pi0b
      wg.sum <- wg1 + wg2 + wg3
      wg.std <- c(wg1, wg2, wg3) / wg.sum
      p_dact <- wg.std[1] * p_a + wg.std[2] * p_b + wg.std[3] * p3
      
      for (fdr in seq_along(FDR_list)) {
        select.idx[[fdr]] <- c(
          select.idx[[fdr]],
          which(p.adjust(p_dact, method = "fdr") <= FDR_list[fdr])
        )
      }
    }
  }
  
  for (fdr in seq_along(FDR_list)) {
    target_FDR <- FDR_list[[fdr]]
    sel_unique <- unique(select.idx[[fdr]])
    sel_all <- select.idx[[fdr]]
    
    tpp_dact <- length(intersect(sel_unique, true_set))
    TPP_dact.gene <- tpp_dact / max(1, length(true_set))
    
    Ec <- setdiff(seq_len(num_trans), true_set)
    
    if (length(sel_unique) == 0) {
      FDR_DACT <- 0
      len_dact <- 0
    } else {
      FDR_DACT <- length(intersect(sel_unique, Ec)) / length(sel_unique)
      len_dact <- length(sel_unique)
    }
    
    if (length(sel_all) == 0) {
      FDR_DACT.pair <- 0
      TPP_dact.pair <- 0
    } else {
      FDR_DACT.pair <- sum(sel_all %in% Ec) / length(sel_all)
      TPP_dact.pair <- sum(sel_all %in% true_set) / (length(true_set) * num_set)
    }
    
    result <- rbind(result, c(TPP_dact.gene, FDR_DACT, FDR_DACT.pair, TPP_dact.pair, len_dact, target_FDR))
  }
  
  colnames(result) <- c("TPP_DACT", "FDR_DACT", "FDR_DACT.pair", "TPP_DACT.pair", "len_dact", "target_FDR")
  list(result = result, select.idx = select.idx)
}

res.per <- function(p_alpha, p_beta, FDR_list, true_set, num_trans, idx = 1) {
  result <- NULL
  num_set <- nrow(p_alpha)
  select.idx <- vector("list", length(FDR_list))
  
  i <- idx
  p_a <- p_alpha[i, ]
  p_b <- p_beta
  
  p_a[p_a == 0] <- min(p_a[p_a != 0])
  p_b[p_b == 0] <- min(p_b[p_b != 0])
  p_a[p_a == 1] <- max(p_a[p_a != 1])
  p_b[p_b == 1] <- max(p_b[p_b != 1])
  
  Z_a <- stats::qnorm(p_a, lower.tail = FALSE)
  Z_b <- stats::qnorm(p_b, lower.tail = FALSE)
  
  pi0a <- 1 - nonnullPropEst(Z_a, 0, 1)
  pi0b <- 1 - nonnullPropEst(Z_b, 0, 1)
  
  if (!is.na(pi0a) && !is.na(pi0b) && pi0a >= 0 && pi0b >= 0) {
    if (pi0a > 1) pi0a <- 1
    if (pi0b > 1) pi0b <- 1
    
    p.mat <- cbind(p_a, p_b)
    p3 <- (apply(p.mat, 1, max))^2
    wg1 <- pi0a * (1 - pi0b)
    wg2 <- (1 - pi0a) * pi0b
    wg3 <- pi0a * pi0b
    wg.sum <- wg1 + wg2 + wg3
    wg.std <- c(wg1, wg2, wg3) / wg.sum
    p_dact <- wg.std[1] * p_a + wg.std[2] * p_b + wg.std[3] * p3
    
    for (fdr in seq_along(FDR_list)) {
      select.idx[[fdr]] <- c(
        select.idx[[fdr]],
        which(p.adjust(p_dact, method = "fdr") <= FDR_list[fdr])
      )
    }
  }
  
  for (fdr in seq_along(FDR_list)) {
    target_FDR <- FDR_list[[fdr]]
    sel_unique <- unique(select.idx[[fdr]])
    sel_all <- select.idx[[fdr]]
    
    tpp_dact <- length(intersect(sel_unique, true_set))
    TPP_dact <- tpp_dact / max(1, length(true_set))
    
    Ec <- setdiff(seq_len(num_trans), true_set)
    
    if (length(sel_unique) == 0) {
      FDR_DACT <- 0
      len_dact <- 0
    } else {
      FDR_DACT <- length(intersect(sel_unique, Ec)) / length(sel_unique)
      len_dact <- length(sel_unique)
    }
    
    if (length(sel_all) == 0) {
      FDR_DACT.pair <- 0
      TPP_dact.pair <- 0
    } else {
      FDR_DACT.pair <- sum(sel_all %in% Ec) / length(sel_all)
      TPP_dact.pair <- sum(sel_all %in% true_set) / (length(true_set) * num_set)
    }
    
    result <- rbind(result, c(TPP_dact, FDR_DACT, FDR_DACT.pair, TPP_dact.pair, len_dact, target_FDR))
  }
  
  colnames(result) <- c("TPP_DACT", "FDR_DACT", "FDR_DACT.pair", "TPP_DACT.pair", "len_dact", "target_FDR")
  list(result = result, select.idx = select.idx)
}

# =========================================================
# 3. fixed simulation settings
# =========================================================
FDR_list <- c(0.01, 0.05, 0.1)

result.all <- NULL
result.per <- NULL

# =========================================================
# 4. main loop
# =========================================================
for (iter in 1:n_sim) {
  repeat {
    cat("iter =", iter, "\n")
    
    n1 <- 900
    n2 <- 200000
    n3 <- 30000
    
    p_snp <- 90
    set_num <- 5
    num_trans <- 1000
    num_set <- p_snp / set_num
    
    true_set <- seq_len(10)
    var_beta <- 0.004
    var_snp <- 0.04
    non_null <- 10
    
    alpha <- matrix(0, nrow = num_set, ncol = num_trans)
    beta  <- rep(0, num_trans)
    
    snp_b <- rnorm(p_snp, 0, sqrt(var_snp))
    b <- matrix(0, nrow = p_snp, ncol = num_set)
    for (i in 1:num_set) {
      b[(5 * i - 4):(5 * i), i] <- snp_b[(5 * i - 4):(5 * i)]
    }
    
    alpha[, 1:10] <- matrix(
      rnorm(num_set * 10, 0, sqrt(var_alpha)),
      nrow = num_set, ncol = 10
    )
    
    beta_true <- rnorm(non_null, 0, sqrt(var_beta / non_null))
    beta[1:10] <- beta_true
    
    pSNP <- runif(p_snp, 0.1, 0.4)
    
    # sample 1
    Z1 <- sapply(pSNP, rbinom, n = n1, size = 2)
    cis_gene1 <- Z1 %*% b + matrix(rnorm(num_set * n1, 0, sqrt(1 - 5 * var_snp)), nrow = n1, ncol = num_set)
    e1 <- matrix(rnorm(num_trans * n1, sd = sqrt(1 - var_alpha * 18)), nrow = n1, ncol = num_trans)
    trans_gene1 <- cis_gene1 %*% alpha + e1
    Y1 <- as.numeric(trans_gene1 %*% beta + rnorm(n1, 0, sqrt(1 - var_beta)))
    
    # sample 2
    Z2 <- sapply(pSNP, rbinom, n = n2, size = 2)
    cis_gene2 <- Z2 %*% b + matrix(rnorm(num_set * n2, 0, sqrt(1 - 5 * var_snp)), nrow = n2, ncol = num_set)
    e2 <- matrix(rnorm(num_trans * n2, sd = sqrt(1 - var_alpha * 18)), nrow = n2, ncol = num_trans)
    trans_gene2 <- cis_gene2 %*% alpha + e2
    Y2 <- as.numeric(trans_gene2 %*% beta + rnorm(n2, 0, sqrt(1 - var_beta)))
    
    # sample 3
    Z3 <- sapply(pSNP, rbinom, n = n3, size = 2)
    cis_gene3 <- Z3 %*% b + matrix(rnorm(num_set * n3, 0, sqrt(1 - 5 * var_snp)), nrow = n3, ncol = num_set)
    e3 <- matrix(rnorm(num_trans * n3, sd = sqrt(1 - var_alpha * 18)), nrow = n3, ncol = num_trans)
    trans_gene3 <- cis_gene3 %*% alpha + e3
    Y3 <- as.numeric(trans_gene3 %*% beta + rnorm(n3, 0, sqrt(1 - var_beta)))
    
    # GWAS p-values
    pval.gwas <- sapply(1:p_snp, function(i) summary(lm(Y2 ~ Z2[, i]))$coefficients[2, 4])
    sig.SNP <- which(p.adjust(pval.gwas, method = "bonferroni") <= 0.05)
    
    if (length(sig.SNP) == 0) next
    
    # ARCHIE small sample
    E1_nor <- scale(trans_gene1)
    Z1_nor <- scale(Z1[, sig.SNP, drop = FALSE])
    sigma_EE <- cor(E1_nor)
    sigma_GE <- cov(Z1_nor, E1_nor)
    sigma_GG <- cor(Z1_nor)
    res.alt.s <- try(
      archie_work(as.matrix(sigma_GE), as.matrix(sigma_GG), as.matrix(sigma_EE), verbose = FALSE),
      silent = TRUE
    )
    
    # ARCHIE larger sample
    E3_nor <- scale(trans_gene3)
    Z3_nor <- scale(Z3[, sig.SNP, drop = FALSE])
    sigma_EE <- cor(E3_nor)
    sigma_GE <- cov(Z3_nor, E3_nor)
    sigma_GG <- cor(Z3_nor)
    res.alt.l <- try(
      archie_work(as.matrix(sigma_GE), as.matrix(sigma_GG), as.matrix(sigma_EE), verbose = FALSE),
      silent = TRUE
    )
    
    if (!inherits(res.alt.s, "try-error") && !inherits(res.alt.l, "try-error")) break
  }
  
  # ARCHIE metrics
  selected_gene.s <- res.alt.s$vs[, which.max(res.alt.s$qs)]
  archie_gene.s <- which(selected_gene.s != 0)
  
  TPP_archie.s <- length(intersect(archie_gene.s, true_set)) / max(1, length(true_set))
  Ec <- setdiff(seq_len(num_trans), true_set)
  FDR_archie.s <- ifelse(length(unique(archie_gene.s)) == 0, 0,
                         length(intersect(archie_gene.s, Ec)) / length(unique(archie_gene.s)))
  
  selected_gene.l <- res.alt.l$vs[, which.max(res.alt.l$qs)]
  archie_gene.l <- which(selected_gene.l != 0)
  
  TPP_archie.l <- length(intersect(archie_gene.l, true_set)) / max(1, length(true_set))
  FDR_archie.l <- ifelse(length(unique(archie_gene.l)) == 0, 0,
                         length(intersect(archie_gene.l, Ec)) / length(unique(archie_gene.l)))
  
  len_archie.l <- length(archie_gene.l)
  len_archie.s <- length(archie_gene.s)
  
  # DACT beta-side p-values
  p_beta <- fast_marginal_pvals(trans_gene2, Y2)
  
  # gene1 -> gene2
  p_alpha_gene <- fast_cross_pvals(cis_gene1, trans_gene1)
  result.1 <- res(p_alpha_gene, p_beta, FDR_list, true_set, num_trans)
  per.gene <- res.per(p_alpha_gene, p_beta, FDR_list, true_set, num_trans, idx = 1)
  
  # snp -> gene2
  p_alpha_snp <- fast_cross_pvals(Z3, trans_gene3)
  res.snp <- res(p_alpha_snp, p_beta, FDR_list, true_set, num_trans)
  per.snp <- res.per(p_alpha_snp, p_beta, FDR_list, true_set, num_trans, idx = 1)
  
  result.all <- rbind(
    result.all,
    matrix(
      cbind(
        result.1$result,
        res.snp$result,
        rep(FDR_archie.s, nrow(res.snp$result)),
        rep(TPP_archie.s, nrow(res.snp$result)),
        rep(FDR_archie.l, nrow(res.snp$result)),
        rep(TPP_archie.l, nrow(res.snp$result)),
        rep(len_archie.l, nrow(res.snp$result)),
        rep(len_archie.s, nrow(res.snp$result))
      ),
      ncol = 18,
      nrow = nrow(result.1$result)
    )
  )
  
  result.per <- rbind(
    result.per,
    matrix(
      cbind(per.gene$result, per.snp$result),
      ncol = 12,
      nrow = nrow(per.snp$result)
    )
  )
  
  colnames(result.all) <- c(
    "TPP_DACT.g2g", "FDR_DACT.g2g", "FDR_DACT.g2g.pair", "TPP_DACT.g2g.pair", "len_dact.g2g", "target_FDR.g2g",
    "TPP_DACT.snp2g", "FDR_DACT.snp2g", "FDR_DACT.g2snp.pair", "TPP_DACT.g2snp.pair", "len_dact.snp2g", "target_FDR",
    "FDR_archie.s", "TPP_archie.s", "FDR_archie.l", "TPP_archie.l", "len_archie.l", "len_archie.s"
  )
  
  colnames(result.per) <- c(
    "TPP_per.g2g", "FDR_per.g2g", "FDR_DACT.g2g.pair", "TPP_DACT.g2g.pair", "len_per.g2g", "target_FDR.g2g",
    "TPP_per.snp2g", "FDR_per.snp2g", "FDR_DACT.g2snp.pair", "TPP_DACT.g2snp.pair", "len_per.snp2g", "target_FDR"
  )
  
  cat("\n---------- global gene/snp result ----------\n")
  result.print <- as.data.frame(result.all) %>%
    dplyr::select(
      TPP_DACT.g2g, FDR_DACT.g2g, FDR_DACT.g2g.pair, TPP_DACT.g2g.pair,
      TPP_DACT.snp2g, FDR_DACT.snp2g, FDR_DACT.g2snp.pair, TPP_DACT.g2snp.pair,
      target_FDR
    )
  
  print(result.print %>%
          group_by(target_FDR) %>%
          summarise(across(everything(), mean, na.rm = TRUE)))
  
  cat("\n---------- per gene/snp result ----------\n")
  per.print <- as.data.frame(result.per) %>%
    dplyr::select(
      TPP_per.g2g, FDR_per.g2g, FDR_DACT.g2g.pair, len_per.g2g,
      TPP_per.snp2g, FDR_per.snp2g, len_per.snp2g, target_FDR
    )
  
  print(per.print %>%
          group_by(target_FDR) %>%
          summarise(across(everything(), mean, na.rm = TRUE)))
  
  save(result.all, file = file.path(out_dir, paste0("sim_snp90_alpha_var", var_alpha, "_beta_var",var_beta,"_global.Rdata")))
  save(result.per, file = file.path(out_dir, paste0("sim_snp90_alpha_var", var_alpha,"_beta_var",var_beta, "_per.Rdata")))
  
  gc()
}

cat("Finished var_alpha =", var_alpha," var_beta = ", var_beta, "\n")