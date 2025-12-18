################################################################################
# Simulation: Evaluate DANDELION performance under a cis(gene1) -> trans(gene2) -> trait
#
# This script simulates:
#   - gene1 (cis genes): num_set = num_gene1
#   - gene2 (trans genes / mediators): num_trans
#   - Trait Y generated from gene2 effects beta
#
# Two sample sizes are used:
#   - n1: smaller sample for estimating gene1 -> gene2  (distal-> proximal) association (p_alpha, i.e. p.trans)
#   - n2: larger sample for estimating gene2 -> trait (proximal-> trait) association (p_beta, i.e. p.wes-like)
#
# The DANDELION step (per gene1) combines:
#   p_a = p_alpha[i, ]  (gene1_i -> each gene2_j)
#   p_b = p_beta        (gene2_j -> trait)
# using the SAME core formula as med_gene():
#   Z_a, Z_b -> pi0a, pi0b -> p3 = (max(p_a, p_b))^2 -> mixture weights -> p_dandelion
#
# Difference vs med_gene():
#   - med_gene() computes qvalues via qvalue::qvalue() and thresholds by target.fdr
#   - here we use BH/FDR adjustment (p.adjust(..., method="fdr")) for simplicity
################################################################################

## ================================
## Required packages
## ================================
library(bigstatsr)   # big_univLinReg, as_FBM
library(stats)

## ================================
## Global parameters
## ================================
n_sim     <- 100
num_gene1 <- 500  # number of gene1 (cis genes / candidate regulators)

FDR_list <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4)

## Result matrix (one row per (iter, target_FDR) appended)
result <- matrix(NA, nrow = 1, ncol = 5)
colnames(result) <- c("TPP_DANDELION", "FDR_DANDELION", "len_dandelion",
                      "num_gene1", "target_FDR")

## ================================
## Simulation loop
## ================================
for (iter in 1:n_sim) {
  
  cat("iter =", iter, "\n")
  
  ## ----------------------------
  ## Sample sizes
  ## ----------------------------
  n1 <- 900       # sample size for gene1->gene2 association (p_alpha)
  n2 <- 200000    # sample size for gene2->trait association (p_beta)
  
  ## ----------------------------
  ## Dimensions
  ## ----------------------------
  set_num   <- 5          # number of trans genes influenced per causal gene1 (block size)
  num_trans <- 1000       # number of gene2
  num_set   <- num_gene1  # number of gene1
  c         <- 1
  
  ## ----------------------------
  ## Truth definition (for evaluation)
  ## IMPORTANT: In your current code, true_set is defined on gene2 index space.
  ## Here true_set = 1:50 are treated as truly disease-relevant gene2.
  ## ----------------------------
  true_set <- 1:50
  
  ## ----------------------------
  ## Effect sizes: alpha (gene1 -> gene2), beta (gene2 -> trait)
  ## ----------------------------
  alpha <- matrix(0, nrow = num_set, ncol = num_trans)  # gene1 -> gene2
  beta  <- rep(0, num_trans)                            # gene2 -> trait
  
  ## Make first 20 gene1 causal: each gene1 affects 5 gene2 in a contiguous block
  for (i in 1:20) {
    alpha[i, (5*i - 4):(5*i)] <- 0.3
  }
  
  ## beta: causal gene2 include (1:50) and (101:160) per your original design
  beta_true <- rnorm(110, 0, sqrt(0.1))
  beta[c(1:50, 101:160)] <- beta_true
  
  ## ============================================================
  ## Generate data: small sample (n1) for p_alpha (gene1 -> gene2)
  ## ============================================================
  cis_gene1 <- matrix(rnorm(n1 * num_set, 0, 0.3), nrow = n1, ncol = num_set)
  e1 <- matrix(rnorm(num_trans * n1, sd = sqrt(0.1)), nrow = n1, ncol = num_trans)
  
  ## gene2 in sample 1
  trans_gene1 <- cis_gene1 %*% alpha + e1
  
  ## trait in sample 1 (not used downstream except optional)
  Y1 <- trans_gene1 %*% beta + rnorm(n1, 0, sqrt(0.1))
  
  ## ============================================================
  ## Generate data: large sample (n2) for p_beta (gene2 -> trait)
  ## ============================================================
  cis_gene2 <- matrix(rnorm(n2 * num_set, 0, 0.3), nrow = n2, ncol = num_set)
  e2 <- matrix(rnorm(num_trans * n2, sd = sqrt(0.1)), nrow = n2, ncol = num_trans)
  
  trans_gene2 <- cis_gene2 %*% alpha + e2
  Y2 <- trans_gene2 %*% beta + rnorm(n2, 0, sqrt(0.1))
  
  ## ============================================================
  ## Step A: gene2 -> trait p-values (p_beta) using large sample
  ## This corresponds to p.wes in med_gene(): p_b = gene2 -> trait p-values
  ## ============================================================
  fit <- big_univLinReg(as_FBM(trans_gene2), Y2)
  p_beta <- 2 * pnorm(-abs(fit[, 3]))  # p-values for gene2 -> trait
  
  ## ============================================================
  ## Step B: gene1 -> gene2 p-values (p_alpha) using small sample
  ## This corresponds to p.trans in med_gene(): p_a = gene1 -> gene2 p-values
  ##
  ## NOTE: p_alpha is (gene1 × gene2). This matches med_gene() where p.trans is
  ##       (gene2 × gene1). Here we store row=gene1, col=gene2.
  ## ============================================================
  p_alpha <- matrix(0, nrow = num_set, ncol = num_trans)
  
  for (j in 1:num_trans) {
    for (i in 1:num_set) {
      fit_alpha <- lm(trans_gene1[, j] ~ cis_gene1[, i])
      p_alpha[i, j] <- summary(fit_alpha)$coef[2, 4]
    }
  }
  
  ## ============================================================
  ## Step C: DANDELION per-gene1 scan across all gene2
  ##
  ## This block implements the SAME core computation as med_gene():
  ##   p_a, p_b -> Z_a, Z_b -> pi0a, pi0b -> p3 -> weights -> p_dandelion
  ##
  ## ============================================================
  select.idx <- vector("list", length(FDR_list))  # accumulated discoveries across gene1
  
  for (i in 1:num_set) {
    
    ## For gene1 = i:
    ## p_a: gene1_i -> each gene2_j (matches med_gene(): p.trans.new[gene.trans, i])
    ## p_b: gene2_j -> trait         (matches med_gene(): p.wes.new[gene.trans])
    p_a <- p_alpha[i, ]
    p_b <- p_beta
    
    ## Numerical stabilization (same intent as med_gene() setting 1 -> 0.99)
    p_a[p_a == 0] <- min(p_a[p_a != 0])
    p_b[p_b == 0] <- min(p_b[p_b != 0])
    p_a[p_a == 1] <- max(p_a[p_a != 1])
    p_b[p_b == 1] <- max(p_b[p_b != 1])
    
    ## Estimate null proportions (same as med_gene())
    Z_a  <- qnorm(p_a, lower.tail = FALSE)
    Z_b  <- qnorm(p_b, lower.tail = FALSE)
    pi0a <- 1 - nonnullPropEst(Z_a, 0, 1)
    pi0b <- 1 - nonnullPropEst(Z_b, 0, 1)
    
    ## Skip if pi0 estimates are invalid
    if (!is.na(pi0a) && !is.na(pi0b) && pi0a >= 0 && pi0b >= 0) {
      
      pi0a <- min(pi0a, 1)
      pi0b <- min(pi0b, 1)
      
      ## p3 term and mixture weights (identical structure to med_gene())
      p3 <- (pmax(p_a, p_b))^2
      
      wg1 <- pi0a * (1 - pi0b)
      wg2 <- (1 - pi0a) * pi0b
      wg3 <- pi0a * pi0b
      wg.sum <- wg1 + wg2 + wg3
      wg.std <- c(wg1, wg2, wg3) / wg.sum
      
      ## DANDELION combined p-values 
      p_dandelion <- wg.std[1] * p_a + wg.std[2] * p_b + wg.std[3] * p3
      
      ## Threshold using BH FDR at multiple target_FDR values
      for (fdr in seq_along(FDR_list)) {
        target_FDR <- FDR_list[fdr]
        select.idx[[fdr]] <- c(select.idx[[fdr]],
                               which(p.adjust(p_dandelion, method = "fdr") <= target_FDR))
      }
    }
  }
  
  ## ============================================================
  ## Step D: Compute TPP and empirical FDR across accumulated discoveries
  ## ============================================================
  for (fdr in seq_along(FDR_list)) {
    
    target_FDR <- FDR_list[fdr]
    
    discovered <- unique(select.idx[[fdr]])
    
    ## TPP: proportion of true gene2 (true_set) recovered
    tpp_dandelion <- length(intersect(discovered, true_set))
    TPP_dandelion <- tpp_dandelion / max(1, length(true_set))
    
    ## Empirical FDR among discoveries (on gene2 index space)
    FDR_DANDELION <- length(setdiff(discovered, true_set)) / max(1, length(discovered))
    
    len_dandelion <- length(discovered)
    
    result <- rbind(result, c(TPP_dandelion, FDR_DANDELION,
                              len_dandelion, num_gene1, target_FDR))
    
    ## Print running mean summary by target_FDR
    result1 <- na.omit(result)
    result1 <- as.data.frame(result1)
    result1$target_FDR <- as.numeric(result1$target_FDR)
    
    print(aggregate(. ~ target_FDR, data = result1, FUN = mean))
  }
  
  save(result, file = paste0("result_FDR", "_num", num_gene1, ".Rdata"))
}
