#!/usr/bin/env Rscript

rm(list = ls())
gc()

## =========================
## 0) arguments
## =========================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1) {
  i <- as.integer(args[1])
} else {
  stop("Please provide index i, e.g. Rscript run_dact_midway.R 2")
}

cat("Running i =", i, "\n")

## =========================
## 1) libraries
## =========================


options(repos = c(CRAN = "https://cloud.r-project.org"))

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(plyr)
  library(tidyr)
  library(tidyverse)
  library(readr)
  library(gprofiler2)
  library(VennDiagram)
  library(qvalue)
})

if ("package:plyr" %in% search()) {
  detach("package:plyr", unload = TRUE)
}

## 如果你原来要 source 帮助函数，这里保留
## source("/project/xuanyao/peixin/review/rare_allel/help_func.R")

cat("Loaded packages successfully.\n")

## =========================
## 2) paths
## =========================
BASE_DIR <- getwd()
UK_dir   <- BASE_DIR

trait_code <- get_arg("--trait_code")
trait_dir  <- file.path(BASE_DIR, trait_code)

trans_eqtl_file <- file.path(
  BASE_DIR,
  "2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz"
)

trans_rdata_file <- file.path(
  BASE_DIR,
  "trans-eQTL_genetype_genename.Rdata"
)

gene_pos_file <- file.path(
  BASE_DIR,
  "gene_position.txt"
)

cat("BASE_DIR =", BASE_DIR, "\n")
cat("trait_dir =", trait_dir, "\n")
cat("trans_eqtl_file =", trans_eqtl_file, "\n")
cat("trans_rdata_file =", trans_rdata_file, "\n")
cat("gene_pos_file =", gene_pos_file, "\n")

## =========================
## 3) helper checks
## =========================
check_file <- function(path) {
  if (!file.exists(path)) {
    stop("File does not exist: ", path)
  } else {
    cat("Found file: ", path, "\n")
  }
}

check_dir <- function(path) {
  if (!dir.exists(path)) {
    stop("Directory does not exist: ", path)
  } else {
    cat("Found directory: ", path, "\n")
  }
}

check_dir(BASE_DIR)
check_dir(UK_dir)
check_dir(trait_dir)

check_file(trans_eqtl_file)
check_file(trans_rdata_file)
check_file(gene_pos_file)

## =========================
## 4) trait file selection
##    only keep:
##    - starts with trait_code_
##    - ends with .Rdata
##    exclude:
##    - _pair.Rdata
##    - _DACT.Rdata
## =========================
trait_files <- list.files(
  trait_dir,
  pattern = paste0("^", trait_code, "_.*\\.Rdata$"),
  full.names = FALSE
)

if (length(trait_files) == 0) {
  stop("No trait-specific .Rdata files found in: ", trait_dir)
}

trait_files <- sort(trait_files)

trait_files_use <- trait_files[
  !grepl("_pair\\.Rdata$|_DACT\\.Rdata$", trait_files)
]

if (length(trait_files_use) == 0) {
  stop("No usable trait .Rdata files left after excluding _pair/_DACT.")
}

allel_list <- sub(
  paste0("^", trait_code, "_(.*)\\.Rdata$"),
  "\\1",
  trait_files_use
)

map_df <- data.frame(
  i = seq_along(trait_files_use),
  file = trait_files_use,
  allele = allel_list,
  stringsAsFactors = FALSE
)

cat("Available trait files:\n")
print(trait_files)

cat("Usable trait files:\n")
print(trait_files_use)

cat("Index-file mapping:\n")
print(map_df)

if (i < 1 || i > length(trait_files_use)) {
  stop("i = ", i, " exceeds valid range 1:", length(trait_files_use))
}

input_trait_file <- file.path(trait_dir, trait_files_use[i])
check_file(input_trait_file)

cat("Selected input file:\n", input_trait_file, "\n")
cat("Selected allele:\n", allel_list[i], "\n")

## =========================
## 5) functions
## =========================
JCCorrect <- function(pval) {
  z <- stats::qnorm(pval, lower.tail = FALSE)
  res <- nullParaEst(z)
  pval.JC <- stats::pnorm(z, mean = res$mu, sd = res$s, lower.tail = FALSE)
  return(pval.JC)
}

nullParaEst <- function(x, gamma = 0.1) {
  n <- length(x)
  t <- c(1:1000) / 200
  
  gan <- n^(-gamma)
  phiplus <- rep(1, 1000)
  phiminus <- rep(1, 1000)
  dphiplus <- rep(1, 1000)
  dphiminus <- rep(1, 1000)
  phi <- rep(1, 1000)
  
  for (ii in 1:1000) {
    s <- t[ii]
    phiplus[ii]   <- mean(cos(s * x))
    phiminus[ii]  <- mean(sin(s * x))
    dphiplus[ii]  <- -mean(x * sin(s * x))
    dphiminus[ii] <- mean(x * cos(s * x))
    phi[ii]       <- sqrt(phiplus[ii]^2 + phiminus[ii]^2)
  }
  
  ind <- min(c(1:1000)[(phi - gan) <= 0])
  tt <- t[ind]
  a  <- phiplus[ind]
  b  <- phiminus[ind]
  da <- dphiplus[ind]
  db <- dphiminus[ind]
  c0 <- phi[ind]
  
  shat <- -(a * da + b * db) / (tt * c0 * c0)
  shat <- sqrt(shat)
  uhat <- -(da * b - db * a) / (c0 * c0)
  
  return(list(mu = uhat, s = shat))
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




med_gene <- function(p.trans, p.wgs, ref.table, trans.p, target.fdr = 0.05, dist = 5e6) {
  
  cat("\n----- Start matching gene 2 -----\n")
  
  common_gene2 <- intersect(rownames(p.trans), names(p.wgs))
  
  ref.table.keep <- ref.table[ref.table$type %in% c("lincRNA", "protein_coding"), ]
  ref.table.keep <- ref.table.keep[!(ref.table.keep$Chromosome %in% c("chrM", "chrX", "chrY")), ]
  ref.table.keep <- ref.table.keep[!duplicated(ref.table.keep$gene_name), ]
  
  common_gene2 <- intersect(common_gene2, ref.table.keep$gene_name)
  
  p.trans.new <- p.trans[rownames(p.trans) %in% common_gene2, , drop = FALSE]
  p.wgs.new   <- p.wgs[names(p.wgs) %in% common_gene2]
  
  p.trans.new <- p.trans.new[match(names(p.wgs.new), rownames(p.trans.new)), , drop = FALSE]
  
  cat("Matched gene2 count:", length(common_gene2), "\n")
  cat("----- End matching gene 2 -----\n")
  
  mat.sig <- matrix(0, nrow = nrow(p.trans.new), ncol = ncol(p.trans.new))
  rownames(mat.sig) <- rownames(p.trans.new)
  colnames(mat.sig) <- colnames(p.trans.new)
  
  mat.p <- matrix(NA, nrow = nrow(p.trans.new), ncol = ncol(p.trans.new))
  rownames(mat.p) <- rownames(p.trans.new)
  colnames(mat.p) <- colnames(p.trans.new)
  
  snp.candidate <- c()
  M <- ncol(p.trans.new)
  
  cat("\n----- Start DACT for ", M, " SNPs -----\n", sep = "")
  
  progress_step <- max(1, floor(M / 100))
  
  for (idx in 1:M) {
    
    cand.snp <- colnames(p.trans.new)[idx]
    
    if (cand.snp %in% trans.p$SNP) {
      
      pos.info <- trans.p[trans.p$SNP == cand.snp, ]
      start.pos <- unique(pos.info$SNPPos)
      end.pos   <- unique(pos.info$SNPPos)
      chr       <- unique(pos.info$SNPChr)
      
      chr <- paste0("chr", chr)
      
      start.critical <- start.pos - dist
      end.critical   <- end.pos + dist
      
      target.genes <- ref.table[ref.table$Chromosome == chr, ]
      
      loc.pos.1 <- which(target.genes$start > start.critical & target.genes$start < end.critical)
      loc.pos.2 <- which(target.genes$end > start.critical & target.genes$end < end.critical)
      loc.pos   <- union(loc.pos.1, loc.pos.2)
      
      gene.cis   <- target.genes$gene_name[loc.pos]
      gene.trans <- setdiff(common_gene2, gene.cis)
      
      if (length(gene.trans) > 0) {
        p_a <- p.trans.new[gene.trans, idx]
        names(p_a) <- gene.trans
        p_a[p_a == 1] <- 0.99
        
        p_b <- p.wgs.new[gene.trans]
        p_b[p_b == 1] <- 0.99
        
        Z_a <- stats::qnorm(p_a, lower.tail = FALSE)
        Z_b <- stats::qnorm(p_b, lower.tail = FALSE)
        
        pi0a <- 1 - nonnullPropEst(Z_a, 0, 1)
        pi0b <- 1 - nonnullPropEst(Z_b, 0, 1)
        
        if (!is.na(pi0a) && !is.na(pi0b) && pi0a >= 0 && pi0b >= 0) {
          
          pi0a <- min(pi0a, 1)
          pi0b <- min(pi0b, 1)
          
          p.mat <- cbind(p_a, p_b)
          p3 <- (apply(p.mat, 1, max))^2
          
          wg1 <- pi0a * (1 - pi0b)
          wg2 <- (1 - pi0a) * pi0b
          wg3 <- pi0a * pi0b
          wg.sum <- wg1 + wg2 + wg3
          wg.std <- c(wg1, wg2, wg3) / wg.sum
          
          p_dact <- wg.std[1] * p_a + wg.std[2] * p_b + wg.std[3] * p3
          
          s <- try(qvalue(p_dact), silent = TRUE)
          if ("try-error" %in% class(s)) {
            s <- qvalue(p_dact, pi0 = 1)
          }
          
          q_dact <- s$qvalues
          names(q_dact) <- gene.trans
          
          vec <- rep(0, nrow(p.trans.new))
          sig_dact <- names(q_dact)[which(q_dact <= target.fdr)]
          vec[rownames(p.trans.new) %in% sig_dact] <- 1
          mat.sig[, idx] <- vec
          
          vec.p <- rep(NA, nrow(p.trans.new))
          vec.p[match(gene.trans, rownames(p.trans.new))] <- p_dact
          mat.p[, idx] <- vec.p
          
          snp.candidate <- c(snp.candidate, cand.snp)
        }
      }
    }
    
    if ((idx %% progress_step) == 0 || idx == M) {
      percent <- round(idx / M * 100, 1)
      cat("Complete ", percent, "% SNPs\n", sep = "")
    }
  }
  
  return(list(
    snp.candidate = snp.candidate,
    mat.sig = mat.sig,
    mat.p = mat.p
  ))
}

calc_pair <- function(mat.sig, snps, p.trans, p.wgs,
                      eta.wgs = 1e-5,
                      eta.trans = 2.664087e-07) {
  
  mat.sig.new <- mat.sig[, snps, drop = FALSE]
  mat.sig.new <- mat.sig.new[rowSums(mat.sig.new) != 0, colSums(mat.sig.new) != 0, drop = FALSE]
  
  cat("\nRemoving self-connected genes\n")
  
  if (nrow(mat.sig.new) > 0 && ncol(mat.sig.new) > 0) {
    for (ii in 1:nrow(mat.sig.new)) {
      for (jj in 1:ncol(mat.sig.new)) {
        if (rownames(mat.sig.new)[ii] == colnames(mat.sig.new)[jj]) {
          mat.sig.new[ii, jj] <- 0
        }
      }
    }
  }
  
  mat.sig.new <- mat.sig.new[rowSums(mat.sig.new) != 0, colSums(mat.sig.new) != 0, drop = FALSE]
  
  pairs_dact <- data.frame(
    gene1 = character(),
    gene2 = character(),
    trans_p = character(),
    stringsAsFactors = FALSE
  )
  
  if (nrow(mat.sig.new) > 0 && ncol(mat.sig.new) > 0) {
    for (ii in 1:nrow(mat.sig.new)) {
      for (jj in 1:ncol(mat.sig.new)) {
        if (mat.sig.new[ii, jj] != 0) {
          pairs_dact[nrow(pairs_dact) + 1, ] <- c(
            colnames(mat.sig.new)[jj],
            rownames(mat.sig.new)[ii],
            p.trans[rownames(mat.sig.new)[ii], colnames(mat.sig.new)[jj]]
          )
        }
      }
    }
  }
  
  if (nrow(pairs_dact) == 0) {
    return(list(
      pairs_dact = pairs_dact,
      num_pair_sig_gene2_wgs = 0,
      num_pair_sig_gene2_gene1_wgs = 0,
      num_pair_sig_gene2_wgs_trans = 0,
      num_pair_sig_gene1_wgs = 0,
      num_pair_sig_gene1_gene2_wgs = 0,
      num_pair_sig_gene1_wgs_trans = 0
    ))
  }
  
  pairs_dact$trans_p <- as.numeric(pairs_dact$trans_p)
  pairs_dact$wgs_gene1 <- p.wgs[match(pairs_dact$gene1, names(p.wgs))]
  pairs_dact$wgs_gene2 <- p.wgs[match(pairs_dact$gene2, names(p.wgs))]
  
  num_pair_sig_gene2_wgs <- sum(na.omit(pairs_dact$wgs_gene2 <= eta.wgs))
  num_pair_sig_gene2_gene1_wgs <- sum(na.omit(pairs_dact$wgs_gene2 <= eta.wgs & pairs_dact$wgs_gene1 <= eta.wgs))
  num_pair_sig_gene2_wgs_trans <- sum(na.omit(pairs_dact$wgs_gene2 <= eta.wgs & pairs_dact$trans_p <= eta.trans))
  
  num_pair_sig_gene1_wgs <- sum(na.omit(pairs_dact$wgs_gene1 <= eta.wgs))
  num_pair_sig_gene1_gene2_wgs <- sum(na.omit(pairs_dact$wgs_gene2 <= eta.wgs & pairs_dact$wgs_gene1 <= eta.wgs))
  num_pair_sig_gene1_wgs_trans <- sum(na.omit(pairs_dact$wgs_gene1 <= eta.wgs & pairs_dact$trans_p <= eta.trans))
  
  return(list(
    pairs_dact = pairs_dact,
    num_pair_sig_gene2_wgs = num_pair_sig_gene2_wgs,
    num_pair_sig_gene2_gene1_wgs = num_pair_sig_gene2_gene1_wgs,
    num_pair_sig_gene2_wgs_trans = num_pair_sig_gene2_wgs_trans,
    num_pair_sig_gene1_wgs = num_pair_sig_gene1_wgs,
    num_pair_sig_gene1_gene2_wgs = num_pair_sig_gene1_gene2_wgs,
    num_pair_sig_gene1_wgs_trans = num_pair_sig_gene1_wgs_trans
  ))
}

## 你自己的 med_gene / calc_pair / 其他函数
## 如果在 help_func.R 里，就 source 进来
## source("/project/xuanyao/peixin/review/rare_allel/help_func.R")

## =========================
## 6) read input
## =========================
cat("\nReading trans-eQTL summary...\n")
trans.p <- fread(trans_eqtl_file, sep = "\t")

out <- strsplit(as.character(trans.p$GeneSymbol), "\\.")
info.out <- do.call(rbind, out)
trans.p$GeneSymbol <- info.out[, 1]

cat("Loading trans-eQTL Rdata...\n")
load(trans_rdata_file)

if (!exists("pval.sub")) {
  stop("Object 'pval.sub' not found in: ", trans_rdata_file)
}

p.trans <- t(pval.sub)

cat("Loading trait-specific data...\n")
load(input_trait_file)

if (!exists("dat_temp")) {
  stop("Object 'dat_temp' not found in: ", input_trait_file)
}

head(dat_temp)

out.pheno <- strsplit(as.character(dat_temp$Name), "\\.")
info.out.pheno <- do.call(rbind, out.pheno)

out.new <- strsplit(as.character(info.out.pheno[, 1]), "\\(")
info.out.new <- do.call(rbind, out.new)

dat_temp$Name <- info.out.new[, 1]
dat_temp_new <- dat_temp[!dat_temp$Name %in% names(which(table(dat_temp$Name) != 1)), ]

cau <- data.frame(pval = dat_temp_new$p_value)
rownames(cau) <- dat_temp_new$Name

p.wgs <- cau[, 1]
p.wgs <- na.omit(p.wgs)
names(p.wgs) <- dat_temp_new$Name

ref.table <- read.delim(gene_pos_file, sep = "\t")

threshold <- 0.05 / nrow(cau)

cat("trait_code =", trait_code, "\n")
cat("allele =", allel_list[i], "\n")
cat("threshold =", threshold, "\n")

## =========================
## 7) run analysis
## =========================
cat("\nRunning med_gene() ...\n")
result <- med_gene(
  p.trans = p.trans,
  p.wgs = p.wgs,
  ref.table = ref.table,
  trans.p = trans.p,
  target.fdr = 0.05
)

cat("\nRunning calc_pair() ...\n")
result.pair <- calc_pair(
  mat.sig = result$mat.sig,
  snps = result$snp.candidate,
  p.trans = p.trans,
  p.wgs = p.wgs,
  eta.wgs = threshold
)

## =========================
## 8) save output
## =========================
out_dact_file <- file.path(
  trait_dir,
  paste0(trait_code, "_", allel_list[i], "_DACT.Rdata")
)

out_pair_file <- file.path(
  trait_dir,
  paste0(trait_code, "_", allel_list[i], "_pair.Rdata")
)

out_thre_file <- file.path(
  trait_dir,
  paste0(trait_code, "_thre.txt")
)

save(result, file = out_dact_file)
save(result.pair, file = out_pair_file)

write.table(
  threshold,
  file = out_thre_file,
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)

cat("\nFinished successfully.\n")
cat("Saved:\n")
cat(out_dact_file, "\n")
cat(out_pair_file, "\n")
cat(out_thre_file, "\n")
