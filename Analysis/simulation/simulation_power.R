

################################################################################
# Simulation study for multi-layer mediation testing
#
# Goal:
#   Evaluate empirical power of DACT-based mediation tests under a
#   two-layer mediation structure:
#
#       X  →  M1  →  M2  →  Y
#
#   Three mediation paths are tested:
#     (1) X  →  M1 → Y
#     (2) X  →  M2 → Y
#     (3) M1 →  M2 → Y
#
#   Power is evaluated under scenarios where signal strength varies
#   in mediation effects (alpha paths).
################################################################################

## ================================
## Required packages
## ================================
library(stats)
library(bigstatsr)      # for big_univLinReg
library(data.table)

## ================================
## Simulation parameters
## ================================
n  <- 500      # total sample size
n1 <- 300      # subsample size
p  <- 500      # number of mediators
set.seed(123)

## True signal indices
true_idx <- 1:80
true_p   <- length(true_idx)

## ================================
## Data generating process
## ================================
sim_dat <- function(n, p, c, change = c("alpha", "beta")) {
  
  ## Initialize effect vectors
  alpha_1 <- rep(0, p)
  alpha_2 <- rep(0, p)
  beta    <- rep(0, p)
  
  ## Configure signal pattern
  if (change == "alpha") {
    
    alpha_true <- rnorm(100, 0, sqrt(c / 50))
    
    alpha_1[1:100]            <- alpha_true
    alpha_2[c(1:90, 101:110)] <- alpha_true
    
    beta[c(1:80, 111:130)] <- 0.5
  }
  
  if (change == "beta") {
    
    alpha_1[1:100]            <- 0.5
    alpha_2[c(1:90, 101:110)] <- 0.5
    
    beta_true <- rnorm(100, 0, sqrt(c / 5))
    beta[c(1:80, 111:130)] <- beta_true
  }
  
  ## Exposure
  X <- rnorm(n, 0, 0.3)
  
  ## First-layer mediators: X → M1
  M_1 <- matrix(0, n, p)
  for (i in 1:n) {
    e <- rnorm(p, sd = sqrt(0.8))
    M_1[i, ] <- X[i] * alpha_1 + e
  }
  
  ## Second-layer mediators: M1 → M2
  M_2 <- matrix(0, n, p)
  for (i in 1:n) {
    e <- rnorm(p, sd = sqrt(0.8))
    M_2[i, ] <- M_1[i, ] * alpha_2 + e
  }
  
  ## Outcome: M2 → Y
  e_y   <- rnorm(n, 0, sqrt(0.5))
  pheno <- M_2 %*% beta + e_y
  
  return(list(X = X, M_1 = M_1, M_2 = M_2, pheno = pheno))
}






########## Useful Functions for DANDELION ##########


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


## ================================
## Storage for results
## ================================
result <- data.frame(
  power_M1 = numeric(),
  power_M2 = numeric(),
  power_M3 = numeric(),
  c        = numeric()
)

## ================================
## Simulation loop
## ================================
for (c in c(0.25)) {
  for (iter in 95:100) {
    
    cat("Running simulation: c =", c, "iter =", iter, "\n")
    
    ## Generate data (alpha-perturbation scenario)
    dat   <- sim_dat(n, p, c, change = "alpha")
    X     <- dat$X
    M_1   <- dat$M_1
    M_2   <- dat$M_2
    pheno <- dat$pheno
    
    ## Subsample individuals
    sub    <- sample(1:n, n1, replace = FALSE)
    X_sub  <- X[sub]
    M1_sub <- M_1[sub, ]
    M2_sub <- M_2[sub, ]
    
    ## =========================================================
    ## Mediation test 1: X → M1 → Y
    ## =========================================================
    p_alpha.1 <- rep(NA, p)
    for (j in 1:p) {
      fit <- lm(M1_sub[, j] ~ X_sub)
      p_alpha.1[j] <- summary(fit)$coef[2, 4]
    }
    
    fit <- big_univLinReg(as_FBM(M_1), pheno)
    p_beta.1 <- 2 * pnorm(-abs(fit[, 3]))
    
    p_a <- p_alpha.1
    p_b <- p_beta.1
    
    ## Numerical stabilization
    p_a[p_a == 0] <- min(p_a[p_a != 0])
    p_b[p_b == 0] <- min(p_b[p_b != 0])
    p_a[p_a == 1] <- max(p_a[p_a != 1])
    p_b[p_b == 1] <- max(p_b[p_b != 1])
    
    Z_a  <- qnorm(p_a, lower.tail = FALSE)
    Z_b  <- qnorm(p_b, lower.tail = FALSE)
    pi0a <- 1 - nonnullPropEst(Z_a, 0, 1)
    pi0b <- 1 - nonnullPropEst(Z_b, 0, 1)
    
    pi0a <- min(max(pi0a, 0), 1)
    pi0b <- min(max(pi0b, 0), 1)
    
    p3 <- (pmax(p_a, p_b))^2
    wg <- c(pi0a * (1 - pi0b),
            (1 - pi0a) * pi0b,
            pi0a * pi0b)
    wg <- wg / sum(wg)
    
    p_dact <- wg[1] * p_a + wg[2] * p_b + wg[3] * p3
    sel1   <- which(p.adjust(p_dact, "BH") <= 0.05)
    power_M1 <- length(intersect(sel1, true_idx)) / true_p
    
    ## =========================================================
    ## Mediation test 2: X → M2 → Y
    ## =========================================================
    p_alpha.2 <- rep(NA, p)
    for (j in 1:p) {
      fit <- lm(M2_sub[, j] ~ X_sub)
      p_alpha.2[j] <- summary(fit)$coef[2, 4]
    }
    
    fit <- big_univLinReg(as_FBM(M_2), pheno)
    p_beta.2 <- 2 * pnorm(-abs(fit[, 3]))
    
    p_a <- p_alpha.2
    p_b <- p_beta.2
    
    p_a[p_a == 0] <- min(p_a[p_a != 0])
    p_b[p_b == 0] <- min(p_b[p_b != 0])
    p_a[p_a == 1] <- max(p_a[p_a != 1])
    p_b[p_b == 1] <- max(p_b[p_b != 1])
    
    Z_a  <- qnorm(p_a, lower.tail = FALSE)
    Z_b  <- qnorm(p_b, lower.tail = FALSE)
    pi0a <- min(max(1 - nonnullPropEst(Z_a, 0, 1), 0), 1)
    pi0b <- min(max(1 - nonnullPropEst(Z_b, 0, 1), 0), 1)
    
    p3 <- (pmax(p_a, p_b))^2
    wg <- c(pi0a * (1 - pi0b),
            (1 - pi0a) * pi0b,
            pi0a * pi0b)
    wg <- wg / sum(wg)
    
    p_dact <- wg[1] * p_a + wg[2] * p_b + wg[3] * p3
    sel2   <- which(p.adjust(p_dact, "BH") <= 0.05)
    power_M2 <- length(intersect(sel2, true_idx)) / true_p
    
    ## =========================================================
    ## Mediation test 3: M1 → M2 → Y
    ## =========================================================
    p_alpha.3 <- rep(NA, p)
    for (j in 1:p) {
      fit <- lm(M2_sub[, j] ~ M1_sub[, j])
      p_alpha.3[j] <- summary(fit)$coef[2, 4]
    }
    
    p_beta.3 <- p_beta.2
    
    p_a <- p_alpha.3
    p_b <- p_beta.3
    
    p_a[p_a == 0] <- min(p_a[p_a != 0])
    p_b[p_b == 0] <- min(p_b[p_b != 0])
    p_a[p_a == 1] <- max(p_a[p_a != 1])
    p_b[p_b == 1] <- max(p_b[p_b != 1])
    
    Z_a  <- qnorm(p_a, lower.tail = FALSE)
    Z_b  <- qnorm(p_b, lower.tail = FALSE)
    pi0a <- min(max(1 - nonnullPropEst(Z_a, 0, 1), 0), 1)
    pi0b <- min(max(1 - nonnullPropEst(Z_b, 0, 1), 0), 1)
    
    p3 <- (pmax(p_a, p_b))^2
    wg <- c(pi0a * (1 - pi0b),
            (1 - pi0a) * pi0b,
            pi0a * pi0b)
    wg <- wg / sum(wg)
    
    p_dact <- wg[1] * p_a + wg[2] * p_b + wg[3] * p3
    sel3   <- which(p.adjust(p_dact, "BH") <= 0.05)
    power_M3 <- length(intersect(sel3, true_idx)) / true_p
    
    ## =========================================================
    ## Store results
    ## =========================================================
    result <- rbind(
      result,
      data.frame(
        power_M1 = power_M1,
        power_M2 = power_M2,
        power_M3 = power_M3,
        c        = c
      )
    )
    
    save(result, file = "power_alpha_var0005.Rdata")
  }
}
