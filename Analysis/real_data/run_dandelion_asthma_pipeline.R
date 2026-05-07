#!/usr/bin/env Rscript

# Analysis/real_data/run_dandelion_asthma_pipeline.R
#
# SNP-based DANDELION pipeline template for asthma / trait-specific analyses
# -------------------------------------------------------------------------
# This script demonstrates how to run SNP-based DANDELION using the exported
# functions from the DANDELION R package:
#   - med_gene()
#   - calc_pair.snp()
#
# This public GitHub version removes local HPC paths, local library settings,
# package installation code, and manually redefined DANDELION functions.
# Users should install DANDELION before running this script.
#
# Usage example:
#   Rscript examples/run_dandelion_asthma_pipeline.R \
#     --trait_code GCST90085447 \
#     --trait_dir path/to/trait_folder \
#     --trans_eqtl_file path/to/trans-eQTLs.txt.gz \
#     --trans_rdata_file path/to/trans-eQTL_genetype_genename.Rdata \
#     --gene_pos_file path/to/gene_position.txt \
#     --trait_index 1 \
#     --out_dir path/to/output \
#     --target_fdr 0.05 \
#     --dist 5000000 \
#     --n_cores 1
#
# Expected input objects:
#   1. trans_rdata_file should contain an object named `pval.sub`.
#      After transposition, `p.trans <- t(pval.sub)` should have:
#        rows    = candidate disease-proximal genes (gene2)
#        columns = SNP exposures
#
#   2. trait-specific RData file should contain an object named `dat_temp`.
#      `dat_temp` should contain columns:
#        Name
#        p_value
#
#   3. trans_eqtl_file should contain SNP annotation columns:
#        SNP
#        SNPPos
#        SNPChr
#      If a column `GeneSymbol` exists, it will be used as the SNP-to-cis-gene
#      mapping for `uniq_snp`.
#
#   4. gene_pos_file should contain:
#        gene_name
#        type
#        Chromosome
#        start
#        end
#
# Output:
#   <trait_code>_<allele>_DACT.Rdata
#   <trait_code>_<allele>_pair.Rdata
#   <trait_code>_thre.txt
#   <trait_code>_<allele>_summary.csv

suppressPackageStartupMessages({
  library(DANDELION)
  library(data.table)
})

# -----------------------------
# 1. Command-line arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) stop("Missing value for ", flag, call. = FALSE)
  args[idx + 1]
}

trait_code <- get_arg("--trait_code", "GCST90085447")
trait_dir <- get_arg("--trait_dir", NULL)
trans_eqtl_file <- get_arg("--trans_eqtl_file", NULL)
trans_rdata_file <- get_arg("--trans_rdata_file", NULL)
gene_pos_file <- get_arg("--gene_pos_file", NULL)
out_dir <- get_arg("--out_dir", trait_dir)

trait_index <- as.integer(get_arg("--trait_index", 1))
target_fdr <- as.numeric(get_arg("--target_fdr", 0.05))
dist <- as.numeric(get_arg("--dist", 5e6))
n_cores <- as.integer(get_arg("--n_cores", 1))
dry_run <- as.logical(get_arg("--dry_run", FALSE))

if (is.null(trait_dir)) stop("Please provide --trait_dir.", call. = FALSE)
if (is.null(trans_eqtl_file)) stop("Please provide --trans_eqtl_file.", call. = FALSE)
if (is.null(trans_rdata_file)) stop("Please provide --trans_rdata_file.", call. = FALSE)
if (is.null(gene_pos_file)) stop("Please provide --gene_pos_file.", call. = FALSE)
if (is.na(trait_index) || trait_index < 1) stop("--trait_index must be a positive integer.", call. = FALSE)
if (is.na(target_fdr) || target_fdr <= 0 || target_fdr >= 1) {
  stop("--target_fdr must be between 0 and 1.", call. = FALSE)
}
if (is.na(dist) || dist <= 0) stop("--dist must be positive.", call. = FALSE)
if (is.na(n_cores) || n_cores < 1) stop("--n_cores must be a positive integer.", call. = FALSE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("DANDELION asthma / trait-specific SNP pipeline")
message("  trait_code       = ", trait_code)
message("  trait_dir        = ", trait_dir)
message("  trans_eqtl_file  = ", trans_eqtl_file)
message("  trans_rdata_file = ", trans_rdata_file)
message("  gene_pos_file    = ", gene_pos_file)
message("  out_dir          = ", out_dir)
message("  trait_index      = ", trait_index)
message("  target_fdr       = ", target_fdr)
message("  dist             = ", dist)
message("  n_cores          = ", n_cores)
message("  dry_run          = ", dry_run)

# -----------------------------
# 2. Input checks
# -----------------------------
check_file <- function(path, label = path) {
  if (!file.exists(path)) stop("File does not exist: ", label, " = ", path, call. = FALSE)
  invisible(TRUE)
}

check_dir <- function(path, label = path) {
  if (!dir.exists(path)) stop("Directory does not exist: ", label, " = ", path, call. = FALSE)
  invisible(TRUE)
}

check_dir(trait_dir, "trait_dir")
check_file(trans_eqtl_file, "trans_eqtl_file")
check_file(trans_rdata_file, "trans_rdata_file")
check_file(gene_pos_file, "gene_pos_file")

# -----------------------------
# 3. Helper functions
# -----------------------------
clean_gene_name <- function(x) {
  x <- as.character(x)
  x <- sub("\\..*$", "", x)
  x <- sub("\\(.*$", "", x)
  trimws(x)
}

get_trait_files <- function(trait_dir, trait_code) {
  files <- list.files(
    trait_dir,
    pattern = paste0("^", trait_code, "_.*\\.[Rr]data$|^", trait_code, "_.*\\.[Rr]Data$"),
    full.names = FALSE
  )

  files <- files[
    !grepl("_pair\\.[Rr]data$|_pair\\.[Rr]Data$|_DACT\\.[Rr]data$|_DACT\\.[Rr]Data$", files)
  ]

  sort(files)
}

load_object_from_rdata <- function(path, object_name) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!(object_name %in% ls(env))) {
    stop("Object '", object_name, "' not found in: ", path, call. = FALSE)
  }
  get(object_name, envir = env)
}

prepare_trait_pvalues <- function(dat_temp) {
  required_cols <- c("Name", "p_value")
  if (!all(required_cols %in% colnames(dat_temp))) {
    stop("dat_temp must contain columns: Name and p_value.", call. = FALSE)
  }

  dat_temp$Name <- clean_gene_name(dat_temp$Name)

  duplicated_names <- names(which(table(dat_temp$Name) != 1))
  dat_keep <- dat_temp[!(dat_temp$Name %in% duplicated_names), , drop = FALSE]

  p.wes <- dat_keep$p_value
  names(p.wes) <- dat_keep$Name
  p.wes <- p.wes[!is.na(p.wes)]

  p.wes
}

prepare_ref_table <- function(gene_pos_file) {
  ref.table <- utils::read.delim(gene_pos_file, sep = "\t", stringsAsFactors = FALSE)
  required_cols <- c("gene_name", "type", "Chromosome", "start", "end")
  if (!all(required_cols %in% colnames(ref.table))) {
    stop(
      "gene_pos_file must contain columns: ",
      paste(required_cols, collapse = ", "),
      call. = FALSE
    )
  }
  ref.table$gene_name <- clean_gene_name(ref.table$gene_name)
  ref.table
}

prepare_snp_ref <- function(trans.p) {
  required_cols <- c("SNP", "SNPPos", "SNPChr")
  if (!all(required_cols %in% colnames(trans.p))) {
    stop(
      "trans_eqtl_file must contain SNP annotation columns: ",
      paste(required_cols, collapse = ", "),
      call. = FALSE
    )
  }

  snp.ref <- unique(trans.p[, ..required_cols])
  as.data.frame(snp.ref)
}

prepare_uniq_snp <- function(trans.p) {
  if (!all(c("SNP", "GeneSymbol") %in% colnames(trans.p))) {
    warning(
      "trans_eqtl_file does not contain both SNP and GeneSymbol. ",
      "Using SNP as a placeholder GeneSymbol mapping. For real analyses, ",
      "provide a biologically meaningful SNP-to-cis-gene mapping."
    )
    return(data.frame(
      SNP = unique(trans.p$SNP),
      GeneSymbol = unique(trans.p$SNP),
      stringsAsFactors = FALSE
    ))
  }

  out <- unique(trans.p[, .(SNP, GeneSymbol)])
  out$GeneSymbol <- clean_gene_name(out$GeneSymbol)
  as.data.frame(out)
}

make_summary_table <- function(result, result.pair, trait_code, allele, target_fdr, threshold) {
  data.frame(
    trait_code = trait_code,
    allele = allele,
    target_fdr = target_fdr,
    wes_threshold = threshold,
    n_exposures_with_valid_results = length(result$gene1),
    n_significant_pairs = ifelse(is.null(result.pair$pairs_dact), 0, nrow(result.pair$pairs_dact)),
    n_locus_level_pairs = ifelse(is.null(result.pair$gene.pair), 0, nrow(result.pair$gene.pair)),
    n_wes_significant_gene2 = length(result.pair$sig_gene2),
    n_non_wes_significant_gene2 = length(result.pair$non_sig.gene2),
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# 4. Read inputs
# -----------------------------
message("Reading trans-eQTL summary table ...")
trans.p <- data.table::fread(trans_eqtl_file, sep = "\t")

if ("GeneSymbol" %in% colnames(trans.p)) {
  trans.p$GeneSymbol <- clean_gene_name(trans.p$GeneSymbol)
}

SNP.ref <- prepare_snp_ref(trans.p)
uniq_snp <- prepare_uniq_snp(trans.p)

message("Loading trans p-value matrix from RData ...")
pval.sub <- load_object_from_rdata(trans_rdata_file, "pval.sub")
p.trans <- t(pval.sub)
p.trans <- as.matrix(p.trans)

message("Reading gene annotation table ...")
ref.table <- prepare_ref_table(gene_pos_file)
ref.table.keep <- ref.table[
  ref.table$type %in% c("lincRNA", "protein_coding") &
    !(ref.table$Chromosome %in% c("chrM", "chrX", "chrY")),
  ,
  drop = FALSE
]
ref.table.keep <- ref.table.keep[!duplicated(ref.table.keep$gene_name), , drop = FALSE]

trait_files <- get_trait_files(trait_dir, trait_code)
if (length(trait_files) == 0) {
  stop("No usable trait RData files found in: ", trait_dir, call. = FALSE)
}
if (trait_index > length(trait_files)) {
  stop(
    "trait_index = ", trait_index,
    " exceeds number of usable trait files = ", length(trait_files),
    call. = FALSE
  )
}

allele_list <- sub(
  paste0("^", trait_code, "_(.*)\\.[Rr][Dd]ata$"),
  "\\1",
  trait_files
)

trait_file <- file.path(trait_dir, trait_files[trait_index])
allele <- allele_list[trait_index]

message("Selected trait file: ", trait_file)
message("Selected allele/mask: ", allele)

message("Loading trait-specific data ...")
dat_temp <- load_object_from_rdata(trait_file, "dat_temp")
p.wes <- prepare_trait_pvalues(dat_temp)

threshold <- 0.05 / length(p.wes)

message("Prepared inputs:")
message("  dim(p.trans)       = ", paste(dim(p.trans), collapse = " x "))
message("  length(p.wes)      = ", length(p.wes))
message("  nrow(ref.table)    = ", nrow(ref.table))
message("  nrow(SNP.ref)      = ", nrow(SNP.ref))
message("  nrow(uniq_snp)     = ", nrow(uniq_snp))
message("  WES threshold      = ", threshold)

if (dry_run) {
  message("Dry run completed. Exiting before running DANDELION.")
  quit(save = "no", status = 0)
}

# -----------------------------
# 5. Run SNP-based DANDELION
# -----------------------------
message("Running med_gene() using gene1.type = 'SNP' ...")
result <- med_gene(
  p.trans = p.trans,
  p.wes = p.wes,
  ref.table = ref.table,
  gene1.list = colnames(p.trans),
  target.fdr = target_fdr,
  dist = dist,
  gene1.type = "SNP",
  SNP.ref = SNP.ref,
  n.cores = n_cores,
  verbose = TRUE
)

message("Organizing significant SNP-gene pairs with calc_pair.snp() ...")
result.pair <- calc_pair.snp(
  mat.sig = result$mat.sig,
  mat.p = result$mat.p,
  p.wes = p.wes,
  gene1 = result$gene1,
  uniq_snp = uniq_snp,
  ref.table.keep = ref.table.keep,
  eta.wgs = threshold,
  SNP.ref = SNP.ref,
  verbose = TRUE
)

summary_df <- make_summary_table(
  result = result,
  result.pair = result.pair,
  trait_code = trait_code,
  allele = allele,
  target_fdr = target_fdr,
  threshold = threshold
)

print(summary_df)

# -----------------------------
# 6. Save outputs
# -----------------------------
out_dact_file <- file.path(
  out_dir,
  paste0(trait_code, "_", allele, "_DACT.Rdata")
)

out_pair_file <- file.path(
  out_dir,
  paste0(trait_code, "_", allele, "_pair.Rdata")
)

out_threshold_file <- file.path(
  out_dir,
  paste0(trait_code, "_thre.txt")
)

out_summary_file <- file.path(
  out_dir,
  paste0(trait_code, "_", allele, "_summary.csv")
)

save(result, file = out_dact_file)
save(result.pair, file = out_pair_file)

write.table(
  threshold,
  file = out_threshold_file,
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)

utils::write.csv(summary_df, file = out_summary_file, row.names = FALSE)

message("Finished successfully.")
message("Saved files:")
message("  ", out_dact_file)
message("  ", out_pair_file)
message("  ", out_threshold_file)
message("  ", out_summary_file)
