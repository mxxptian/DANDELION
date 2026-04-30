# examples/simulation_compare_template.R
#
# DANDELION simulation template
# ------------------------------------------------------------
# This script demonstrates how to organize a simulation workflow around
# the DANDELION package without including the full manuscript-scale
# data-generation mechanism.
#
# The public version intentionally uses a lightweight toy-data generator
# only to demonstrate the required input format, function calls, and
# output summaries. Replace `simulate_dandelion_inputs()` with your own
# simulation design when running a full study.
#
# Usage:
#   Rscript examples/simulation_template.R 0.0025 10
#
# Arguments:
#   1. var_alpha: user-defined simulation parameter passed to the toy generator
#   2. n_sim:     number of simulation replicates
#
# Output:
#   examples/results/simulation_template_alpha<var_alpha>.rds

suppressPackageStartupMessages({
  library(DANDELION)
})

# -----------------------------
# 1. User arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

var_alpha <- if (length(args) >= 1) as.numeric(args[1]) else 0.0025
n_sim <- if (length(args) >= 2) as.integer(args[2]) else 10

if (is.na(var_alpha) || var_alpha <= 0) {
  stop("var_alpha must be a positive numeric value.")
}
if (is.na(n_sim) || n_sim <= 0) {
  stop("n_sim must be a positive integer.")
}

message("Running DANDELION simulation template")
message("  var_alpha = ", var_alpha)
message("  n_sim     = ", n_sim)

# -----------------------------
# 2. Public toy settings
# -----------------------------
# These settings are intentionally small and are not intended to reproduce
# the manuscript-scale simulation. They only define dimensions for a toy
# example that demonstrates the DANDELION interface.
settings <- list(
  num_gene1 = 5,
  num_snp = 5,
  num_gene2 = 20,
  num_true_gene2 = 3,
  fdr_levels = c(0.01, 0.05, 0.1)
)

# -----------------------------
# 3. Toy input generator
# -----------------------------
# This function is a placeholder for users' own simulation design.
#
# It returns objects in the same format required by DANDELION:
#   p.trans.gene: gene1 -> gene2 trans-association p-value matrix
#   p.trans.snp:  SNP -> gene2 trans-association p-value matrix
#   p.wes:        named gene-level disease association p-values
#   ref.table:    gene annotation table
#   SNP.ref:      SNP annotation table
#   uniq_snp:     SNP-to-cis-gene mapping table
#   true_gene2:   known positive gene2 set for evaluating toy performance
#
# The full manuscript-scale data-generation mechanism is not included in
# this public template during peer review.
simulate_dandelion_inputs <- function(seed, var_alpha, settings) {
  set.seed(seed)

  gene2 <- paste0("G", seq_len(settings$num_gene2))
  gene1 <- paste0("E", seq_len(settings$num_gene1))
  snps <- paste0("rs", seq_len(settings$num_snp))

  true_gene2 <- gene2[seq_len(settings$num_true_gene2)]

  # Generate toy p-values for demonstration only.
  # The first few gene2 values are assigned smaller p-values so that the
  # template produces visible selected signals in small examples.
  p.trans.gene <- matrix(
    stats::runif(settings$num_gene2 * settings$num_gene1, 0.05, 1),
    nrow = settings$num_gene2,
    ncol = settings$num_gene1,
    dimnames = list(gene2, gene1)
  )

  p.trans.snp <- matrix(
    stats::runif(settings$num_gene2 * settings$num_snp, 0.05, 1),
    nrow = settings$num_gene2,
    ncol = settings$num_snp,
    dimnames = list(gene2, snps)
  )

  signal_strength <- max(1e-6, min(0.05, var_alpha * 10))
  p.trans.gene[true_gene2, ] <- stats::runif(
    length(true_gene2) * settings$num_gene1,
    min = 1e-6,
    max = signal_strength
  )
  p.trans.snp[true_gene2, ] <- stats::runif(
    length(true_gene2) * settings$num_snp,
    min = 1e-6,
    max = signal_strength
  )

  p.wes <- stats::runif(settings$num_gene2, 0.05, 1)
  names(p.wes) <- gene2
  p.wes[true_gene2] <- stats::runif(length(true_gene2), 1e-6, signal_strength)

  # Minimal gene annotation table.
  # gene2 genes are placed on chr1; exposure genes are placed on chr2 so that
  # they are treated as distal in this toy example.
  ref.table <- data.frame(
    gene_name = c(gene2, gene1),
    type = "protein_coding",
    Chromosome = c(rep("chr1", length(gene2)), rep("chr2", length(gene1))),
    start = seq_len(length(gene2) + length(gene1)) * 1e7,
    end = seq_len(length(gene2) + length(gene1)) * 1e7 + 1000,
    stringsAsFactors = FALSE
  )

  SNP.ref <- data.frame(
    SNP = snps,
    SNPPos = seq_len(length(snps)) * 2e7,
    SNPChr = "2",
    stringsAsFactors = FALSE
  )

  uniq_snp <- data.frame(
    SNP = snps,
    GeneSymbol = gene1[seq_along(snps)],
    stringsAsFactors = FALSE
  )

  list(
    p.trans.gene = p.trans.gene,
    p.trans.snp = p.trans.snp,
    p.wes = p.wes,
    ref.table = ref.table,
    SNP.ref = SNP.ref,
    uniq_snp = uniq_snp,
    gene1.list = gene1,
    snp.list = snps,
    true_gene2 = true_gene2
  )
}

# -----------------------------
# 4. DANDELION wrappers
# -----------------------------
make_ref_table_keep <- function(ref.table) {
  ref.table[
    ref.table$type %in% c("lincRNA", "protein_coding") &
      !(ref.table$Chromosome %in% c("chrM", "chrX", "chrY")),
    ,
    drop = FALSE
  ]
}

run_gene_based_dandelion <- function(dat, target.fdr) {
  fit <- med_gene(
    p.trans = dat$p.trans.gene,
    p.wes = dat$p.wes,
    ref.table = dat$ref.table,
    gene1.list = dat$gene1.list,
    gene1.type = "Gene",
    target.fdr = target.fdr,
    dist = 5e6,
    verbose = FALSE
  )

  calc_pair.gene(
    mat.sig = fit$mat.sig,
    mat.p = fit$mat.p,
    p.wes = dat$p.wes,
    gene1 = fit$gene1,
    ref.table.keep = make_ref_table_keep(dat$ref.table),
    eta.wgs = 1e-5,
    verbose = FALSE
  )
}

run_snp_based_dandelion <- function(dat, target.fdr) {
  fit <- med_gene(
    p.trans = dat$p.trans.snp,
    p.wes = dat$p.wes,
    ref.table = dat$ref.table,
    gene1.list = dat$snp.list,
    gene1.type = "SNP",
    SNP.ref = dat$SNP.ref,
    target.fdr = target.fdr,
    dist = 5e6,
    verbose = FALSE
  )

  calc_pair.snp(
    mat.sig = fit$mat.sig,
    mat.p = fit$mat.p,
    p.wes = dat$p.wes,
    gene1 = fit$gene1,
    uniq_snp = dat$uniq_snp,
    ref.table.keep = make_ref_table_keep(dat$ref.table),
    SNP.ref = dat$SNP.ref,
    eta.wgs = 1e-5,
    verbose = FALSE
  )
}

# -----------------------------
# 5. Evaluation helper
# -----------------------------
evaluate_dandelion_pairs <- function(pair.res, true_gene2) {
  if (is.null(pair.res$gene.pair) || nrow(pair.res$gene.pair) == 0) {
    selected_gene2 <- character(0)
  } else {
    selected_gene2 <- unique(pair.res$gene.pair$gene2)
  }

  tp <- length(intersect(selected_gene2, true_gene2))
  fp <- length(setdiff(selected_gene2, true_gene2))

  data.frame(
    n_selected = length(selected_gene2),
    TPP = tp / length(true_gene2),
    FDR = ifelse(length(selected_gene2) == 0, 0, fp / length(selected_gene2))
  )
}

run_one_replicate <- function(rep_id, var_alpha, settings) {
  dat <- simulate_dandelion_inputs(
    seed = rep_id,
    var_alpha = var_alpha,
    settings = settings
  )

  out <- data.frame()

  for (fdr in settings$fdr_levels) {
    gene_pair <- run_gene_based_dandelion(dat, target.fdr = fdr)
    gene_eval <- evaluate_dandelion_pairs(gene_pair, dat$true_gene2)
    gene_eval$method <- "gene-based DANDELION"
    gene_eval$target_fdr <- fdr

    snp_pair <- run_snp_based_dandelion(dat, target.fdr = fdr)
    snp_eval <- evaluate_dandelion_pairs(snp_pair, dat$true_gene2)
    snp_eval$method <- "SNP-based DANDELION"
    snp_eval$target_fdr <- fdr

    out <- rbind(out, gene_eval, snp_eval)
  }

  out$replicate <- rep_id
  out$var_alpha <- var_alpha
  out
}

# -----------------------------
# 6. Run template simulation
# -----------------------------
sim_results <- do.call(
  rbind,
  lapply(seq_len(n_sim), function(i) {
    message("Simulation replicate ", i, " / ", n_sim)
    run_one_replicate(
      rep_id = i,
      var_alpha = var_alpha,
      settings = settings
    )
  })
)

summary_results <- stats::aggregate(
  cbind(TPP, FDR, n_selected) ~ method + target_fdr + var_alpha,
  data = sim_results,
  FUN = mean
)

print(summary_results)

# -----------------------------
# 7. Save output
# -----------------------------
out_dir <- file.path("examples", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(
  out_dir,
  paste0("simulation_template_alpha", var_alpha, ".rds")
)

saveRDS(
  list(
    settings = settings,
    replicate_results = sim_results,
    summary_results = summary_results
  ),
  file = out_file
)

message("Saved results to: ", out_file)
