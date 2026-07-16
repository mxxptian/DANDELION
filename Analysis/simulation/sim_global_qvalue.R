#!/usr/bin/env Rscript

rm(list = ls())
gc()

args <- commandArgs(trailingOnly = TRUE)

## -------------------------
## parse arguments
## -------------------------
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) stop(paste("Missing value for", flag))
  args[idx + 1]
}

num_gene1 <- as.integer(get_arg("--num_gene1", 18))
var_alpha <- as.numeric(get_arg("--var_alpha", 0.0025))
var_beta  <- as.numeric(get_arg("--var_beta", 0.1))
n_sim     <- as.integer(get_arg("--n_sim", 100))
seed      <- as.integer(get_arg("--seed", 123))
n1        <- as.integer(get_arg("--n1", 900))
n2        <- as.integer(get_arg("--n2", 200000))
num_trans <- as.integer(get_arg("--num_trans", 1000))
non_null  <- as.integer(get_arg("--non_null", 10))

cat("Arguments received:\n")
cat("num_gene1 =", num_gene1, "\n")
cat("var_alpha =", var_alpha, "\n")
cat("var_beta  =", var_beta, "\n")
cat("n_sim     =", n_sim, "\n")
cat("seed      =", seed, "\n")
cat("n1        =", n1, "\n")
cat("n2        =", n2, "\n")
cat("num_trans =", num_trans, "\n")
cat("non_null  =", non_null, "\n")

set.seed(seed)

## =========================
## library path on Midway3
## =========================
user_lib <- "./R_libs"
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

cat("Current .libPaths():\n")
print(.libPaths())

if (!requireNamespace("qvalue", quietly = TRUE, lib.loc = user_lib)) {
  stop(paste0("Package 'qvalue' not found in user_lib: ", user_lib))
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(qvalue, lib.loc = user_lib)
})

if ("package:plyr" %in% search()) {
  detach("package:plyr", unload = TRUE)
}

proj_dir <- getwd()


dir.create(file.path(proj_dir, "summary_results"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(proj_dir, "logs"), showWarnings = FALSE, recursive = TRUE)

## =========================
## helper functions
## =========================
fast_marginal_pvals <- function(X, y) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- nrow(X)
  
  X_center <- sweep(X, 2, colMeans(X), FUN = "-")
  y_center <- y - mean(y)
  
  Sxx <- colSums(X_center^2)
  Sxy <- colSums(sweep(X_center, 1, y_center, FUN = "*"))
  
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

safe_fix_p <- function(p) {
  p <- as.numeric(p)
  
  if (all(is.na(p))) return(rep(1, length(p)))
  
  p[is.na(p)] <- 1
  
  if (any(p == 0)) {
    nz_min <- suppressWarnings(min(p[p > 0], na.rm = TRUE))
    if (is.finite(nz_min)) {
      p[p == 0] <- nz_min
    } else {
      p[p == 0] <- .Machine$double.xmin
    }
  }
  
  if (any(p == 1)) {
    lt1_max <- suppressWarnings(max(p[p < 1], na.rm = TRUE))
    if (is.finite(lt1_max)) {
      p[p == 1] <- lt1_max
    } else {
      p[p == 1] <- 1 - 1e-16
    }
  }
  
  p <- pmin(pmax(p, .Machine$double.xmin), 1 - 1e-16)
  p
}

nonnullPropEst <- function(x,u,sigma){
  # x is a vector
  # u is the mean
  # sigma is the standard deviation
  
  z  = (x - u)/sigma
  xi = c(0:100)/100
  tmax=sqrt(log(length(x)))
  tt=seq(0,tmax,0.1)
  
  epsest=NULL
  
  for (j in 1:length(tt)) {
    
    t=tt[j]
    f  = t*xi
    f  = exp(f^2/2)
    w  = (1 - abs(xi))
    co  = 0*xi
    
    for (i in 1:101) {
      co[i] = mean(cos(t*xi[i]*z));
    }
    epshat = 1 - sum(w*f*co)/sum(w)
    epsest=c(epsest,epshat)
  }
  return(epsest= max(epsest))
}


safe_qvalue <- function(p) {
  p <- pmin(pmax(as.numeric(p), .Machine$double.xmin), 1 - 1e-16)
  
  qobj <- tryCatch(
    qvalue::qvalue(p),
    error = function(e) NULL,
    warning = function(w) {
      invokeRestart("muffleWarning")
    }
  )
  
  if (is.null(qobj)) {
    return(list(
      qvalues = p.adjust(p, method = "BH"),
      method = "BH_fallback"
    ))
  }
  
  list(
    qvalues = qobj$qvalues,
    method = "qvalue"
  )
}

## =========================
## simulation settings
## =========================
n_true_snp  <- 18
n_true_gene <- 10
FDR_list <- c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1)

if (num_gene1 < n_true_snp) {
  stop("num_gene1 is smaller than n_true_snp. Please increase num_gene1 or reduce n_true_snp.")
}
if (num_trans < n_true_gene) {
  stop("num_trans is smaller than n_true_gene. Please increase num_trans or reduce n_true_gene.")
}
if (1 - var_alpha * n_true_snp <= 0) {
  stop("1 - var_alpha * n_true_snp must be > 0 for residual SD.")
}
if (1 - var_beta <= 0) {
  stop("1 - var_beta must be > 0 for Y noise SD.")
}

true_set   <- seq_len(n_true_gene)
true_pairs <- expand.grid(snp = seq_len(n_true_snp), gene = seq_len(n_true_gene))
true_key   <- paste(true_pairs$snp, true_pairs$gene, sep = "_")

result <- matrix(NA, nrow = 0, ncol = 9)
colnames(result) <- c(
  "iter",
  "TPP_gene", "FDR_gene", "len_gene", "num_gene1",
  "TPP_pair", "FDR_pair", "len_pair", "target_FDR"
)

tag <- paste0(
  "num", num_gene1,
  "_alpha", format(var_alpha, scientific = FALSE, trim = TRUE),
  "_beta", format(var_beta, scientific = FALSE, trim = TRUE),
  "_seed", seed
)

## =========================
## main loop
## =========================
for (iter in seq_len(n_sim)) {
  
  cat("\n=============================\n")
  cat("iter =", iter, "/", n_sim, "\n")
  cat("=============================\n")
  
  iter_start <- Sys.time()
  
  num_set <- num_gene1
  
  alpha <- matrix(0, nrow = num_set, ncol = num_trans)
  beta  <- rep(0, num_trans)
  
  alpha[seq_len(n_true_snp), seq_len(n_true_gene)] <- matrix(
    rnorm(n_true_snp * n_true_gene, 0, sqrt(var_alpha)),
    nrow = n_true_snp, ncol = n_true_gene
  )
  
  beta_true <- rnorm(non_null, 0, sqrt(var_beta / non_null))
  beta[seq_len(non_null)] <- beta_true
  
  ## sample 1
  cis_gene1 <- matrix(rnorm(n1 * num_set, 0, 0.8), nrow = n1, ncol = num_set)
  e1 <- matrix(
    rnorm(num_trans * n1, sd = sqrt(1 - var_alpha * n_true_snp)),
    nrow = n1, ncol = num_trans
  )
  trans_gene1 <- cis_gene1 %*% alpha + e1
  Y1 <- as.numeric(trans_gene1 %*% beta + rnorm(n1, 0, sqrt(1 - var_beta)))
  
  ## sample 2
  cis_gene2 <- matrix(rnorm(n2 * num_set, 0, 0.8), nrow = n2, ncol = num_set)
  e2 <- matrix(
    rnorm(num_trans * n2, sd = sqrt(1 - var_alpha * n_true_snp)),
    nrow = n2, ncol = num_trans
  )
  trans_gene2 <- cis_gene2 %*% alpha + e2
  Y2 <- as.numeric(trans_gene2 %*% beta + rnorm(n2, 0, sqrt(1 - var_beta)))
  
  ## fast p-values
  p_beta  <- fast_marginal_pvals(trans_gene2, Y2)
  p_alpha <- fast_cross_pvals(cis_gene1, trans_gene1)
  
  ## build all pair-level DACT p-values
  all_pair_list <- vector("list", length = num_set)
  
  for (i in seq_len(num_set)) {
    p_a <- safe_fix_p(p_alpha[i, ])
    p_b <- safe_fix_p(p_beta)
    
    Z_a <- qnorm(p_a, lower.tail = FALSE)
    Z_b <- qnorm(p_b, lower.tail = FALSE)
    
    pi0a <- 1 - nonnullPropEst(Z_a, 0, 1)
    pi0b <- 1 - nonnullPropEst(Z_b, 0, 1)
    
    p_dact <- rep(1, num_trans)
    
    if (!is.na(pi0a) && !is.na(pi0b) && pi0a >= 0 && pi0b >= 0) {
      pi0a <- min(pi0a, 1)
      pi0b <- min(pi0b, 1)
      
      p3 <- pmax(p_a, p_b)^2
      
      wg1 <- pi0a * (1 - pi0b)
      wg2 <- (1 - pi0a) * pi0b
      wg3 <- pi0a * pi0b
      wg_sum <- wg1 + wg2 + wg3
      
      if (wg_sum > 0) {
        wg_std <- c(wg1, wg2, wg3) / wg_sum
        p_dact <- wg_std[1] * p_a + wg_std[2] * p_b + wg_std[3] * p3
      }
    }
    
    all_pair_list[[i]] <- data.frame(
      snp    = i,
      gene   = seq_len(num_trans),
      p_dact = p_dact
    )
  }
  
  all_pair_df <- dplyr::bind_rows(all_pair_list)
  
  ## global qvalue over all pairs
  qres <- safe_qvalue(all_pair_df$p_dact)
  all_pair_df$qvalue <- qres$qvalues
  q_method <- qres$method
  
  ## evaluate at each target FDR
  for (target_FDR in FDR_list) {
    
    sel_pair <- all_pair_df %>%
      dplyr::filter(qvalue <= target_FDR) %>%
      dplyr::select(snp, gene) %>%
      dplyr::distinct()
    
    sel_gene_union <- sort(unique(sel_pair$gene))
    
    ## gene-level
    if (length(sel_gene_union) == 0) {
      TPP_gene <- 0
      FDR_gene <- 0
      len_gene <- 0
    } else {
      TPP_gene <- length(intersect(sel_gene_union, true_set)) / length(true_set)
      FDR_gene <- length(setdiff(sel_gene_union, true_set)) / length(sel_gene_union)
      len_gene <- length(sel_gene_union)
    }
    
    ## pair-level
    if (nrow(sel_pair) == 0) {
      TPP_pair <- 0
      FDR_pair <- 0
      len_pair <- 0
    } else {
      sel_key <- paste(sel_pair$snp, sel_pair$gene, sep = "_")
      TP <- sum(sel_key %in% true_key)
      FP <- sum(!(sel_key %in% true_key))
      
      TPP_pair <- TP / length(true_key)
      FDR_pair <- FP / nrow(sel_pair)
      len_pair <- nrow(sel_pair)
    }
    
    new_row <- c(
      iter,
      TPP_gene, FDR_gene, len_gene, num_gene1,
      TPP_pair, FDR_pair, len_pair, target_FDR
    )
    
    result <- rbind(result, new_row)
  }
  
  ## save cumulative result
  result_df <- as.data.frame(result)
  
  cumulative_file <- file.path(
    proj_dir, "summary_results",
    paste0("gridsearch_", tag, "_cumulative.Rdata")
  )
  
  save(
    result_df,
    iter,
    q_method,
    file = cumulative_file
  )
  
  ## print running summary
  running_summary <- result_df %>%
    dplyr::group_by(target_FDR) %>%
    dplyr::summarise(
      mean_TPP_gene = mean(TPP_gene),
      mean_FDR_gene = mean(FDR_gene),
      mean_len_gene = mean(len_gene),
      mean_TPP_pair = mean(TPP_pair),
      mean_FDR_pair = mean(FDR_pair),
      mean_len_pair = mean(len_pair),
      .groups = "drop"
    )
  
  print(running_summary)
  cat("q-value method used:", q_method, "\n")
  cat("iter runtime:", round(as.numeric(difftime(Sys.time(), iter_start, units = "mins")), 2), "mins\n")
  
  rm(
    alpha, beta, beta_true,
    cis_gene1, cis_gene2,
    trans_gene1, trans_gene2,
    e1, e2, Y1, Y2,
    p_alpha, p_beta,
    all_pair_df, all_pair_list,
    running_summary
  )
  gc()
}

## =========================
## final summary
## =========================
result_df <- as.data.frame(result)

summary_df <- result_df %>%
  dplyr::group_by(target_FDR) %>%
  dplyr::summarise(
    mean_TPP_gene = mean(TPP_gene),
    mean_FDR_gene = mean(FDR_gene),
    mean_len_gene = mean(len_gene),
    mean_TPP_pair = mean(TPP_pair),
    mean_FDR_pair = mean(FDR_pair),
    mean_len_pair = mean(len_pair),
    .groups = "drop"
  ) %>%
  dplyr::mutate(dist_to_0.1 = abs(mean_FDR_pair - 0.1)) %>%
  dplyr::arrange(dist_to_0.1)

print(summary_df)

save(
  result_df,
  summary_df,
  file = file.path(
    proj_dir, "summary_results",
    paste0("gridsearch_", tag, "_final.Rdata")
  )
)

write.csv(
  result_df,
  file = file.path(
    proj_dir, "summary_results",
    paste0("gridsearch_", tag, "_all_rows.csv")
  ),
  row.names = FALSE
)

write.csv(
  summary_df,
  file = file.path(
    proj_dir, "summary_results",
    paste0("gridsearch_", tag, "_summary.csv")
  ),
  row.names = FALSE
)

cat("\nFinished successfully:\n")
cat("tag =", tag, "\n")
