################################################################################
# Description
#
# The simulation implementation (res and res.per) uses the same DANDELION test statistic as med_gene, 
# including identical Z-score transformation, null proportion estimation, mixture weighting, 
# and p-value combination. 
################################################################################

################################################################################
# Function: res
#
# Description:
#   Apply DANDELION to gene1 -> gene2 -> trait associations across
#   all candidate gene1, and evaluate discovery performance.
#
# Consistency with med_gene():
#   - p_a corresponds to p.trans (gene1 -> gene2)
#   - p_b corresponds to p.wes   (gene2 -> trait)
#   - Same Z-score transform, pi0 estimation, mixture weights,
#     and p_DANDELION construction
#
# Difference:
#   - Here: BH-FDR via p.adjust()
#   - med_gene(): qvalue::qvalue()
################################################################################

res <- function(p_trans, p_wes, FDR_list, select_idx,
                true_gene2_idx, n_gene2) {
  
  result <- NULL
  n_gene1 <- nrow(p_trans)
  
  for (i in seq_len(n_gene1)) {
    
    ## gene1_i -> all gene2
    p_a <- p_trans[i, ]
    p_b <- p_wes
    
    ## numerical stabilization (same spirit as med_gene)
    p_a[p_a == 0] <- min(p_a[p_a != 0])
    p_b[p_b == 0] <- min(p_b[p_b != 0])
    p_a[p_a == 1] <- max(p_a[p_a != 1])
    p_b[p_b == 1] <- max(p_b[p_b != 1])
    
    ## ---- DANDELION core statistic (IDENTICAL to med_gene) ----
    Z_a  <- stats::qnorm(p_a, lower.tail = FALSE)
    Z_b  <- stats::qnorm(p_b, lower.tail = FALSE)
    
    pi0a <- 1 - nonnullPropEst(Z_a, 0, 1)
    pi0b <- 1 - nonnullPropEst(Z_b, 0, 1)
    
    if (!any(is.na(c(pi0a, pi0b))) && pi0a >= 0 && pi0b >= 0) {
      
      pi0a <- min(pi0a, 1)
      pi0b <- min(pi0b, 1)
      
      p3 <- (pmax(p_a, p_b))^2
      
      w <- c(pi0a * (1 - pi0b),
             (1 - pi0a) * pi0b,
             pi0a * pi0b)
      w <- w / sum(w)
      
      p_dandelion <- w[1] * p_a +
        w[2] * p_b +
        w[3] * p3
      
      ## FDR thresholding
      for (f in seq_along(FDR_list)) {
        select_idx[[f]] <- c(
          select_idx[[f]],
          which(p.adjust(p_dandelion, method = "fdr") <= FDR_list[f])
        )
      }
    }
  }
  
  ## -------- performance evaluation --------
  for (f in seq_along(FDR_list)) {
    
    target_FDR <- FDR_list[f]
    
    TPP_DANDELION <-
      length(intersect(select_idx[[f]], true_gene2_idx)) /
      max(1, length(true_gene2_idx))
    
    null_gene2 <- setdiff(seq_len(n_gene2), true_gene2_idx)
    
    FDR_DANDELION <-
      length(intersect(select_idx[[f]], null_gene2)) /
      max(1, length(unique(select_idx[[f]])))
    
    n_selected_gene2 <- length(unique(select_idx[[f]]))
    
    result <- rbind(
      result,
      c(TPP_DANDELION, FDR_DANDELION,
        n_selected_gene2, target_FDR)
    )
  }
  
  colnames(result) <- c("TPP_DANDELION",
                        "FDR_DANDELION",
                        "n_selected_gene2",
                        "target_FDR")
  
  return(list(result = result, select_idx = select_idx))
}


################################################################################
# Function: res.per
#
# Description:
#   Apply DANDELION to a SINGLE gene1 (or SNP),
#   corresponding to per-gene / per-SNP power analysis.
#
# Statistical definition is IDENTICAL to res() and med_gene().
################################################################################

res.per <- function(p_trans, p_wes, FDR_list, select_idx,
                    true_gene2_idx, n_gene2, idx) {
  
  result <- NULL
  
  ## fixed gene1
  p_a <- p_trans[idx, ]
  p_b <- p_wes
  
  p_a[p_a == 0] <- min(p_a[p_a != 0])
  p_b[p_b == 0] <- min(p_b[p_b != 0])
  p_a[p_a == 1] <- max(p_a[p_a != 1])
  p_b[p_b == 1] <- max(p_b[p_b != 1])
  
  ## ---- DANDELION core statistic ----
  Z_a  <- stats::qnorm(p_a, lower.tail = FALSE)
  Z_b  <- stats::qnorm(p_b, lower.tail = FALSE)
  
  pi0a <- 1 - nonnullPropEst(Z_a, 0, 1)
  pi0b <- 1 - nonnullPropEst(Z_b, 0, 1)
  
  if (!any(is.na(c(pi0a, pi0b))) && pi0a >= 0 && pi0b >= 0) {
    
    pi0a <- min(pi0a, 1)
    pi0b <- min(pi0b, 1)
    
    p3 <- (pmax(p_a, p_b))^2
    
    w <- c(pi0a * (1 - pi0b),
           (1 - pi0a) * pi0b,
           pi0a * pi0b)
    w <- w / sum(w)
    
    p_dandelion <- w[1] * p_a +
      w[2] * p_b +
      w[3] * p3
    
    for (f in seq_along(FDR_list)) {
      select_idx[[f]] <- c(
        select_idx[[f]],
        which(p.adjust(p_dandelion, method = "fdr") <= FDR_list[f])
      )
    }
  }
  
  ## evaluation
  for (f in seq_along(FDR_list)) {
    
    target_FDR <- FDR_list[f]
    
    TPP_DANDELION <-
      length(intersect(select_idx[[f]], true_gene2_idx)) /
      max(1, length(true_gene2_idx))
    
    null_gene2 <- setdiff(seq_len(n_gene2), true_gene2_idx)
    
    FDR_DANDELION <-
      length(intersect(select_idx[[f]], null_gene2)) /
      max(1, length(unique(select_idx[[f]])))
    
    n_selected_gene2 <- length(unique(select_idx[[f]]))
    
    result <- rbind(
      result,
      c(TPP_DANDELION, FDR_DANDELION,
        n_selected_gene2, target_FDR)
    )
  }
  
  colnames(result) <- c("TPP_DANDELION",
                        "FDR_DANDELION",
                        "n_selected_gene2",
                        "target_FDR")
  
  return(list(result = result, select_idx = select_idx))
}


############################################################
# Useful function from DANDELION
############################################################


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


library(bigstatsr)
library(tidyverse)

set.seed(123)

############################################################
# Global simulation settings
############################################################

n_sim <- 300

FDR_list <- c(0.01, 0.02, 0.03, 0.04,
              0.05, 0.1, 0.2, 0.3, 0.4)

############################################################
# Storage for results
############################################################

result_global   <- NULL   # gene1 -> gene2
result_per_gene <- NULL   # per-gene1

############################################################
# Start simulation loop
############################################################

for (iter in seq_len(n_sim)) {
  
  cat("\n================ Simulation", iter, "================\n")
  
  ##########################################################
  # 1. Sample size and dimension settings
  ##########################################################
  
  n_gene1_gene2 <- 900
  n_gene2_trait <- 200000
  n_snp_gene1   <- 30000
  
  n_snp           <- 90
  n_snp_per_gene1 <- 5
  n_gene1         <- n_snp / n_snp_per_gene1
  n_gene2         <- 1000
  
  ##########################################################
  # 2. Ground truth definition
  ##########################################################
  
  true_gene2_idx  <- seq_len(10)
  n_nonnull_gene2 <- 110
  
  ##########################################################
  # 3. Effect size parameters
  ##########################################################
  
  var_alpha <- 0.001
  var_beta  <- 0.4
  var_snp   <- 0.04
  
  ##########################################################
  # 4. Generate structural parameters
  ##########################################################
  
  # gene1 -> gene2
  alpha_gene1_gene2 <- matrix(0, n_gene1, n_gene2)
  alpha_gene1_gene2[, true_gene2_idx] <-
    matrix(rnorm(n_gene1 * length(true_gene2_idx),
                 0, sqrt(var_alpha)),
           n_gene1, length(true_gene2_idx))
  
  # gene2 -> trait
  beta_gene2_trait <- rep(0, n_gene2)
  beta_gene2_trait[c(true_gene2_idx,
                     101:(n_nonnull_gene2 - length(true_gene2_idx) + 100))] <-
    rnorm(n_nonnull_gene2, 0,
          sqrt(var_beta / n_nonnull_gene2))
  
  # SNP -> gene1
  beta_snp_gene1 <- rnorm(n_snp, 0, sqrt(var_snp))
  snp_to_gene1   <- matrix(0, n_snp, n_gene1)
  
  for (i in seq_len(n_gene1)) {
    idx <- ((i - 1) * n_snp_per_gene1 + 1):(i * n_snp_per_gene1)
    snp_to_gene1[idx, i] <- beta_snp_gene1[idx]
  }
  
  ##########################################################
  # 5. Generate genotypes
  ##########################################################
  
  maf <- runif(n_snp, 0.1, 0.4)
  
  Z_sample1 <- sapply(maf, rbinom, n = n_gene1_gene2, size = 2)
  Z_sample2 <- sapply(maf, rbinom, n = n_gene2_trait, size = 2)
  Z_sample3 <- sapply(maf, rbinom, n = n_snp_gene1,   size = 2)
  
  ##########################################################
  # 6. Generate gene1 expression
  ##########################################################
  
  G_gene1_sample1 <- Z_sample1 %*% snp_to_gene1 +
    matrix(rnorm(n_gene1_gene2 * n_gene1,
                 0, sqrt(1 - n_snp_per_gene1 * var_snp)),
           n_gene1_gene2, n_gene1)
  
  G_gene1_sample2 <- Z_sample2 %*% snp_to_gene1 +
    matrix(rnorm(n_gene2_trait * n_gene1,
                 0, sqrt(1 - n_snp_per_gene1 * var_snp)),
           n_gene2_trait, n_gene1)
  
  ##########################################################
  # 7. Generate gene2 (mediators)
  ##########################################################
  
  M_gene2_sample1 <- G_gene1_sample1 %*% alpha_gene1_gene2 +
    matrix(rnorm(n_gene1_gene2 * n_gene2,
                 0, sqrt(1 - var_alpha * 18)),
           n_gene1_gene2, n_gene2)
  
  M_gene2_sample2 <- G_gene1_sample2 %*% alpha_gene1_gene2 +
    matrix(rnorm(n_gene2_trait * n_gene2,
                 0, sqrt(1 - var_alpha * 18)),
           n_gene2_trait, n_gene2)
  
  ##########################################################
  # 8. Generate trait
  ##########################################################
  
  Y_sample2 <- M_gene2_sample2 %*% beta_gene2_trait +
    rnorm(n_gene2_trait, 0, sqrt(1 - var_beta))
  
  ##########################################################
  # 9. Estimate p_wes (gene2 -> trait)
  ##########################################################
  
  fit_beta <- big_univLinReg(as_FBM(M_gene2_sample2), Y_sample2)
  p_wes <- 2 * pnorm(-abs(fit_beta[, 3]))
  
  ##########################################################
  # 10. Estimate p_trans (gene1 -> gene2)
  ##########################################################
  
  p_trans <- matrix(NA, n_gene1, n_gene2)
  
  for (j in seq_len(n_gene2)) {
    for (i in seq_len(n_gene1)) {
      fit_alpha <- lm(M_gene2_sample1[, j] ~ G_gene1_sample1[, i])
      p_trans[i, j] <- summary(fit_alpha)$coef[2, 4]
    }
  }
  
  ##########################################################
  # 11. Apply DANDELION (global gene1 -> gene2)
  ##########################################################
  
  select_idx <- vector("list", length(FDR_list))
  
  res_global <- res(
    p_trans          = p_trans,
    p_wes            = p_wes,
    FDR_list         = FDR_list,
    select_idx       = select_idx,
    true_gene2_idx   = true_gene2_idx,
    n_gene2          = n_gene2
  )
  
  result_global <- rbind(result_global, res_global$result)
  
  ##########################################################
  # 12. Apply DANDELION (per-gene1 analysis)
  ##########################################################
  
  per_idx <- 1  # evaluate first gene1 as example
  
  select_idx_per <- vector("list", length(FDR_list))
  
  res_per_gene <- res.per(
    p_trans          = p_trans,
    p_wes            = p_wes,
    FDR_list         = FDR_list,
    select_idx       = select_idx_per,
    true_gene2_idx   = true_gene2_idx,
    n_gene2          = n_gene2,
    idx              = per_idx
  )
  
  result_per_gene <- rbind(result_per_gene, res_per_gene$result)
}

############################################################
# 13. Summarize results
############################################################

result_global_df <- as.data.frame(result_global)
result_per_df    <- as.data.frame(result_per_gene)

cat("\n===== Global DANDELION Performance =====\n")
print(
  result_global_df %>%
    group_by(target_FDR) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
)

cat("\n===== Per-gene DANDELION Performance =====\n")
print(
  result_per_df %>%
    group_by(target_FDR) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
)
