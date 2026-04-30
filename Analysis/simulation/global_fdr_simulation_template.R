# examples/global_fdr_simulation_template.R
#
# DANDELION global-FDR simulation template
# ------------------------------------------------------------
# This public template demonstrates how to structure a simulation
# workflow for evaluating DANDELION at both the gene level and the
# exposure-gene pair level.
#
# The full manuscript-scale data-generation mechanism, HPC paths,
# parameter grid, and internal helper scripts are intentionally not
# included in this public example during peer review. Replace
# `simulate_dandelion_inputs()` with your own data-generation design
# when running a full simulation study.
#
# Usage:
#   Rscript examples/global_fdr_simulation_template.R \
#     --num_gene1 100 --var_alpha 0.004 --var_beta 0.1 --n_sim 10 --seed 123
#
# Output:
#   examples/results/global_fdr_template_summary.csv
#   examples/results/global_fdr_template_results.rds

suppressPackageStartupMessages({
  library(DANDELION)
})

# -----------------------------
# 1. Command-line arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) stop(paste("Missing value for", flag))
  args[idx + 1]
}

num_gene1 <- as.integer(get_arg("--num_gene1", 20))
var_alpha <- as.numeric(get_arg("--var_alpha", 0.004))
var_beta <- as.numeric(get_arg("--var_beta", 0.1))
n_sim <- as.integer(get_arg("--n_sim", 10))
seed <- as.integer(get_arg("--seed", 123))
num_gene2 <- as.integer(get_arg("--num_gene2", 100))
num_true_gene1 <- as.integer(get_arg("--num_true_gene1", 5))
num_true_gene2 <- as.integer(get_arg("--num_true_gene2", 5))

if (is.na(num_gene1) || num_gene1 <= 0) stop("num_gene1 must be a positive integer.")
if (is.na(num_gene2) || num_gene2 <= 0) stop("num_gene2 must be a positive integer.")
if (is.na(var_alpha) || var_alpha <= 0) stop("var_alpha must be a positive numeric value.")
if (is.na(var_beta) || var_beta <= 0) stop("var_beta must be a positive numeric value.")
if (is.na(n_sim) || n_sim <= 0) stop("n_sim must be a positive integer.")
if (is.na(seed)) stop("seed must be an integer.")

num_true_gene1 <- min(num_true_gene1, num_gene1)
num_true_gene2 <- min(num_true_gene2, num_gene2)

message("Running DANDELION global-FDR simulation template")
message("  num_gene1      = ", num_gene1)
message("  num_gene2      = ", num_gene2)
message("  var_alpha      = ", var_alpha)
message("  var_beta       = ", var_beta)
message("  n_sim          = ", n_sim)
message("  seed           = ", seed)
message("  num_true_gene1 = ", num_true_gene1)
message("  num_true_gene2 = ", num_true_gene2)

# -----------------------------
# 2. Public toy settings
# -----------------------------
settings <- list(
  num_gene1 = num_gene1,
  num_gene2 = num_gene2,
  num_true_gene1 = num_true_gene1,
  num_true_gene2 = num_true_gene2,
  var_alpha = var_alpha,
  var_beta = var_beta,
  fdr_levels = c(0.01, 0.05, 0.1, 0.2)
)

# -----------------------------
# 3. Toy input generator
# -----------------------------
# This function is a public placeholder. It generates toy p-value matrices
# with known positive gene1-gene2 pairs only to demonstrate how to evaluate
# DANDELION results.
#
# The manuscript-scale simulation used a separate data-generation mechanism
# involving latent effect matrices, independent samples, and large-scale
# p-value generation. That full mechanism is not included in this template.
simulate_dandelion_inputs <- function(seed, settings) {
  set.seed(seed)

  gene1 <- paste0("E", seq_len(settings$num_gene1))
  gene2 <- paste0("G", seq_len(settings$num_gene2))

  true_gene1 <- gene1[seq_len(settings$num_true_gene1)]
  true_gene2 <- gene2[seq_len(settings$num_true_gene2)]

  # Start from null-like p-values.
  p.trans <- matrix(
    stats::runif(settings$num_gene2 * settings$num_gene1, 0.05, 1),
    nrow = settings$num_gene2,
    ncol = settings$num_gene1,
    dimnames = list(gene2, gene1)
  )

  p.wes <- stats::runif(settings$num_gene2, 0.05, 1)
  names(p.wes) <- gene2

  # Add toy signals to the known true set. This is not the manuscript
  # data-generation mechanism; it is only a compact public demonstration.
  signal_alpha <- max(1e-6, min(0.05, settings$var_alpha * 10))
  signal_beta <- max(1e-6, min(0.05, settings$var_beta / 2))

  p.trans[true_gene2, true_gene1] <- stats::runif(
    length(true_gene2) * length(true_gene1),
    min = 1e-6,
    max = signal_alpha
  )
  p.wes[true_gene2] <- stats::runif(
    length(true_gene2),
    min = 1e-6,
    max = signal_beta
  )

  # Minimal annotation table. Put candidate gene2 genes on chr1 and exposure
  # genes on chr2 so that they are treated as distal in this toy example.
  ref.table <- data.frame(
    gene_name = c(gene2, gene1),
    type = "protein_coding",
    Chromosome = c(rep("chr1", length(gene2)), rep("chr2", length(gene1))),
    start = seq_len(length(gene2) + length(gene1)) * 1e7,
    end = seq_len(length(gene2) + length(gene1)) * 1e7 + 1000,
    stringsAsFactors = FALSE
  )

  true_pairs <- expand.grid(
    gene1 = true_gene1,
    gene2 = true_gene2,
    stringsAsFactors = FALSE
  )

  list(
    p.trans = p.trans,
    p.wes = p.wes,
    ref.table = ref.table,
    gene1.list = gene1,
    true_gene1 = true_gene1,
    true_gene2 = true_gene2,
    true_pairs = true_pairs
  )
}

make_ref_table_keep <- function(ref.table) {
  ref.table[
    ref.table$type %in% c("lincRNA", "protein_coding") &
      !(ref.table$Chromosome %in% c("chrM", "chrX", "chrY")),
    ,
    drop = FALSE
  ]
}

# -----------------------------
# 4. Run DANDELION at a target FDR
# -----------------------------
run_dandelion_at_fdr <- function(dat, target.fdr) {
  fit <- med_gene(
    p.trans = dat$p.trans,
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

# -----------------------------
# 5. Gene-level and pair-level metrics
# -----------------------------
evaluate_global_fdr <- function(pair.res, true_gene2, true_pairs) {
  if (is.null(pair.res$gene.pair) || nrow(pair.res$gene.pair) == 0) {
    selected_gene2 <- character(0)
    selected_pairs <- data.frame(gene1 = character(), gene2 = character())
  } else {
    selected_gene2 <- unique(pair.res$gene.pair$gene2)
    selected_pairs <- unique(pair.res$gene.pair[, c("gene1", "gene2"), drop = FALSE])
  }

  # Gene-level metrics: selected gene2 union across all exposures.
  if (length(selected_gene2) == 0) {
    TPP_gene <- 0
    FDR_gene <- 0
    len_gene <- 0
  } else {
    TPP_gene <- length(intersect(selected_gene2, true_gene2)) / length(true_gene2)
    FDR_gene <- length(setdiff(selected_gene2, true_gene2)) / length(selected_gene2)
    len_gene <- length(selected_gene2)
  }

  # Pair-level metrics: selected gene1-gene2 pairs.
  true_key <- paste(true_pairs$gene1, true_pairs$gene2, sep = "_")

  if (nrow(selected_pairs) == 0) {
    TPP_pair <- 0
    FDR_pair <- 0
    len_pair <- 0
  } else {
    selected_key <- paste(selected_pairs$gene1, selected_pairs$gene2, sep = "_")
    TP <- sum(selected_key %in% true_key)
    FP <- sum(!(selected_key %in% true_key))

    TPP_pair <- TP / length(true_key)
    FDR_pair <- FP / nrow(selected_pairs)
    len_pair <- nrow(selected_pairs)
  }

  data.frame(
    TPP_gene = TPP_gene,
    FDR_gene = FDR_gene,
    len_gene = len_gene,
    TPP_pair = TPP_pair,
    FDR_pair = FDR_pair,
    len_pair = len_pair
  )
}

run_one_replicate <- function(rep_id, settings) {
  dat <- simulate_dandelion_inputs(
    seed = rep_id,
    settings = settings
  )

  out <- data.frame()

  for (fdr in settings$fdr_levels) {
    pair.res <- run_dandelion_at_fdr(dat, target.fdr = fdr)
    eval <- evaluate_global_fdr(
      pair.res = pair.res,
      true_gene2 = dat$true_gene2,
      true_pairs = dat$true_pairs
    )

    eval$target_FDR <- fdr
    eval$replicate <- rep_id
    eval$num_gene1 <- settings$num_gene1
    eval$var_alpha <- settings$var_alpha
    eval$var_beta <- settings$var_beta

    out <- rbind(out, eval)
  }

  out
}

# -----------------------------
# 6. Main loop
# -----------------------------
set.seed(seed)

result <- do.call(
  rbind,
  lapply(seq_len(n_sim), function(i) {
    message("Simulation replicate ", i, " / ", n_sim)
    run_one_replicate(i, settings)
  })
)

summary_df <- stats::aggregate(
  cbind(TPP_gene, FDR_gene, len_gene, TPP_pair, FDR_pair, len_pair) ~
    target_FDR + num_gene1 + var_alpha + var_beta,
  data = result,
  FUN = mean
)

summary_df$dist_to_0.1 <- abs(summary_df$FDR_pair - 0.1)
summary_df <- summary_df[order(summary_df$dist_to_0.1), ]

print(summary_df)

# -----------------------------
# 7. Save output
# -----------------------------
out_dir <- file.path("examples", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(
  list(
    settings = settings,
    replicate_results = result,
    summary_results = summary_df
  ),
  file = file.path(out_dir, "global_fdr_template_results.rds")
)

utils::write.csv(
  summary_df,
  file = file.path(out_dir, "global_fdr_template_summary.csv"),
  row.names = FALSE
)

message("Saved results to: ", out_dir)
