#' Apply DANDELION to Identify Candidate Disease Proximal Genes
#'
#' @description
#' Applies the DANDELION procedure to integrate trans-regulatory association
#' p-values and gene-level trait association p-values. The exposure side can be
#' either distal genes or SNPs. For clarity, exposures are referred to as
#' `gene1`, and candidate disease proximal genes are referred to as `gene2`.
#'
#' @param p.trans A numeric matrix of trans-regulatory association p-values.
#'   Rows are candidate disease proximal genes (`gene2`), and columns are
#'   exposures (`gene1`), either distal genes or SNPs. Row names and column names
#'   must be provided. (e.g., from trans-eQTL studies like eQTLGen)
#' @param p.wes A named numeric vector of gene-level trait association p-values (e.g., from WES burden tests or GWAS-based gene-level tests).
#'    Note: Ensure that names (Gene Symbols) match the row names of `p.trans`
#' @param target.fdr False Discovery Rate threshold. DANDELION uses this to decide which gene pairs are statistically significant. Default is 0.1.
#' @param ref.table A data frame containing gene annotation and genomic positions.
#'   Required columns are `gene_name`, `type`, `Chromosome`, `start`, and `end`.
#'   `Chromosome` should use the format `chr1`, `chr2`, etc.
#' @param gene1.list A character vector of candidate exposures to analyze. These
#'   values must be a subset of `colnames(p.trans)`.
#' @param dist  The cis-window size (in base pairs). Genes within this distance from the exposure will be excluded to focus on true distal (trans) effects.
#'   Default is 5e6.
#' @param gene1.type Exposure type. Must be either `"SNP"` or `"Gene"`.
#' @param SNP.ref Optional SNP annotation data frame used when `gene1.type = "SNP"`.
#'   Required columns are `SNP`, `SNPPos`, and `SNPChr`. `SNPChr` can be `1` or
#'   `chr1`; both formats are accepted.
#' @param n.cores Number of cores for parallel execution. On non-Windows systems,
#'   values larger than 1 use `parallel::mclapply()`. On Windows, the function
#'   falls back to single-core execution.
#' @param verbose Logical. If `TRUE`, progress messages are printed. Default is `FALSE`.
#'
#' @return A list with three elements: `gene1`, the analyzed exposures with at
#'   least one valid DANDELION result; `mat.sig`, a matrix encoding significant
#'   gene1-gene2 pairs; and `mat.p`, a matrix of DANDELION p-values.
#'
#' @examples
#' set.seed(1)
#' p.trans <- matrix(runif(60, 0.001, 0.9), nrow = 12, ncol = 5)
#' rownames(p.trans) <- paste0("G", 1:12)
#' colnames(p.trans) <- paste0("E", 1:5)
#' p.wes <- runif(12, 0.001, 0.9)
#' names(p.wes) <- paste0("G", 1:12)
#' ref.table <- data.frame(
#'   gene_name = c(paste0("G", 1:12), paste0("E", 1:5)),
#'   type = "protein_coding",
#'   Chromosome = "chr1",
#'   start = seq_len(17) * 1e7,
#'   end = seq_len(17) * 1e7 + 1000
#' )
#' res <- med_gene(
#'   p.trans = p.trans,
#'   p.wes = p.wes,
#'   ref.table = ref.table,
#'   gene1.list = colnames(p.trans),
#'   gene1.type = "Gene"
#' )
#'
#' @export
med_gene <- function(p.trans,
                     p.wes,
                     ref.table,
                     gene1.list,
                     target.fdr = 0.1,
                     dist = 5e6,
                     gene1.type = c("SNP", "Gene"),
                     SNP.ref = NULL,
                     n.cores = 1,
                     verbose = FALSE) {
  gene1.type <- match.arg(gene1.type)
  gene1.type.lc <- tolower(gene1.type)

  check_med_gene_inputs(
    p.trans = p.trans,
    p.wes = p.wes,
    ref.table = ref.table,
    gene1.list = gene1.list,
    gene1.type = gene1.type,
    SNP.ref = SNP.ref
  )
  
  # Restrict analysis to requested exposures and retain only common gene2 values
  # shared by the trans-association matrix, burden-test p-values, and annotation table.
  gene1.list <- intersect(gene1.list, colnames(p.trans))
  if (length(gene1.list) == 0) {
    stop("No values in gene1.list are present in colnames(p.trans).")
  }

  p.trans <- p.trans[, gene1.list, drop = FALSE]

  dandelion_message("Start matching gene2.", verbose = verbose)

  # Filter the gene annotation table to autosomal lincRNA and protein-coding genes.
  # This table is used both for gene matching and for cis-window exclusion.
  ref.table.keep <- ref.table[ref.table$type %in% c("lincRNA", "protein_coding"), , drop = FALSE]
  ref.table.keep <- ref.table.keep[!(ref.table.keep$Chromosome %in% c("chrM", "chrX", "chrY")), , drop = FALSE]
  ref.table.keep <- ref.table.keep[!duplicated(ref.table.keep$gene_name), , drop = FALSE]

  common_gene2 <- intersect(rownames(p.trans), names(p.wes))
  common_gene2 <- intersect(common_gene2, ref.table.keep$gene_name)

  if (length(common_gene2) == 0) {
    stop("No common gene2 values among rownames(p.trans), names(p.wes), and ref.table$gene_name.")
  }

  p.trans.new <- p.trans[common_gene2, , drop = FALSE]
  p.wes.new <- p.wes[common_gene2]
  p.trans.new <- p.trans.new[match(names(p.wes.new), rownames(p.trans.new)), , drop = FALSE]

  dandelion_message("Finished matching gene2.", verbose = verbose)

  # Initialize output matrices:
  # mat.sig stores whether each exposure-gene2 pair is significant;
  # mat.p stores the corresponding DANDELION p-value.
  mat.sig <- matrix(0L, nrow = nrow(p.trans.new), ncol = ncol(p.trans.new))
  rownames(mat.sig) <- rownames(p.trans.new)
  colnames(mat.sig) <- colnames(p.trans.new)

  mat.p <- matrix(NA_real_, nrow = nrow(p.trans.new), ncol = ncol(p.trans.new))
  rownames(mat.p) <- rownames(p.trans.new)
  colnames(mat.p) <- colnames(p.trans.new)

  M <- ncol(p.trans.new)

  # For gene-based DANDELION, use the genomic coordinates of the exposure gene
  # to remove nearby cis genes before testing trans-mediated paths.
  run_single_gene <- function(i) {
    gene1 <- colnames(p.trans.new)[i]

    if (!(gene1 %in% ref.table.keep$gene_name)) {
      return(NULL)
    }

    pos.info <- ref.table.keep[ref.table.keep$gene_name == gene1, , drop = FALSE]
    pos.info <- pos.info[1, , drop = FALSE]

    gene.cis <- find_cis_genes_for_region(
      chr = pos.info$Chromosome,
      start.pos = pos.info$start,
      end.pos = pos.info$end,
      ref.table = ref.table.keep,
      dist = dist
    )

    run_dandelion_for_exposure(
      exposure.id = gene1,
      exposure.index = i,
      gene.trans = setdiff(common_gene2, gene.cis),
      p.trans.new = p.trans.new,
      p.wes.new = p.wes.new,
      target.fdr = target.fdr
    )
  }

  # For SNP-based DANDELION, use the SNP position to remove genes within the
  # cis-window before testing trans-mediated paths.
  run_single_snp <- function(i) {
    cand.snp <- colnames(p.trans.new)[i]

    if (is.null(SNP.ref)) {
      return(NULL)
    }

    snp.info <- SNP.ref[SNP.ref$SNP == cand.snp, , drop = FALSE]
    if (nrow(snp.info) == 0) {
      return(NULL)
    }

    snp.info <- snp.info[1, , drop = FALSE]
    snp.chr <- normalize_chr(snp.info$SNPChr)
    snp.pos <- as.numeric(snp.info$SNPPos)

    if (is.na(snp.pos)) {
      return(NULL)
    }

    gene.cis <- find_cis_genes_for_region(
      chr = snp.chr,
      start.pos = snp.pos,
      end.pos = snp.pos,
      ref.table = ref.table.keep,
      dist = dist
    )

    run_dandelion_for_exposure(
      exposure.id = cand.snp,
      exposure.index = i,
      gene.trans = setdiff(common_gene2, gene.cis),
      p.trans.new = p.trans.new,
      p.wes.new = p.wes.new,
      target.fdr = target.fdr
    )
  }

  if (gene1.type.lc == "gene") {
    dandelion_message(sprintf("Start DANDELION for %s genes.", M), verbose = verbose)
    worker <- run_single_gene
  } else {
    dandelion_message(sprintf("Start DANDELION for %s SNPs.", M), verbose = verbose)
    worker <- run_single_snp
  }

  # Run one exposure at a time. On Unix-like systems, users can enable
  # multi-core execution through n.cores.
  index.seq <- seq_len(M)
  if (.Platform$OS.type != "windows" && n.cores > 1) {
    results <- parallel::mclapply(index.seq, worker, mc.cores = n.cores)
  } else {
    results <- lapply(index.seq, worker)
  }

  gene1.candidate <- character(0)
  for (res in results) {
    if (!is.null(res)) {
      gene1.candidate <- c(gene1.candidate, res$id)
      mat.sig[, res$id] <- res$sig.vec
      mat.p[, res$id] <- res$p.vec
    }
  }

  if (gene1.type.lc == "gene") {
    dandelion_message(
      sprintf("Number of gene1 detected: %s.", length(gene1.candidate)),
      verbose = verbose
    )
  } else {
    dandelion_message(
      sprintf("Number of SNPs detected: %s.", length(gene1.candidate)),
      verbose = verbose
    )
  }

  out <- list(
    gene1 = gene1.candidate,
    mat.sig = mat.sig,
    mat.p = mat.p,
    gene1.type = gene1.type,
    target.fdr = target.fdr
  )
  class(out) <- "dandelion_result"
  out
}


#' Organize DANDELION SNP-Gene Results
#'
#' @description
#' Cleans significant DANDELION pairs when the exposure side contains SNPs.
#' SNPs are mapped to cis genes using `uniq_snp` first. If a SNP is not present
#' in `uniq_snp`, the function optionally maps it to the nearest or overlapping
#' gene using `SNP.ref` and `ref.table.keep`.
#'
#' @param mat.sig Matrix encoding significant pairs, returned by `med_gene()`.
#' @param mat.p Matrix of DANDELION p-values, returned by `med_gene()`.
#' @param p.wes Named numeric vector of gene-level trait association p-values.
#' @param gene1 Candidate SNPs returned by `med_gene()`.
#' @param uniq_snp Data frame containing known SNP-to-cis-gene mapping. Required
#'   columns are `SNP` and `GeneSymbol`.
#' @param ref.table.keep Gene annotation data frame after filtering. Required
#'   columns are `gene_name`, `Chromosome`, `start`, and `end`.
#' @param eta.wgs Significance threshold for WES genes. Default is 1e-5.
#' @param SNP.ref Optional SNP annotation data frame with columns `SNP`, `SNPPos`,
#'   and `SNPChr`. Used only for SNPs not mapped by `uniq_snp`.
#' @param verbose Logical. If `TRUE`, progress messages are printed. Default is `FALSE`.
#'
#' @return A list containing `pairs_dact`, `gene.pair`, `sig_gene2`, and
#'   `non_sig.gene2`.
#'
#' @export
calc_pair.snp <- function(mat.sig,
                          mat.p,
                          p.wes,
                          gene1,
                          uniq_snp,
                          ref.table.keep,
                          eta.wgs = 1e-5,
                          SNP.ref = NULL,
                          verbose = FALSE) {
  check_pair_inputs(mat.sig, mat.p, p.wes, gene1, ref.table.keep)

  if (!all(c("SNP", "GeneSymbol") %in% colnames(uniq_snp))) {
    stop("uniq_snp must contain columns: SNP and GeneSymbol.")
  }
  
  # Build the significant SNP-gene2 pair table from mat.sig and mat.p.
  pairs_dact <- build_pairs_from_mats(mat.sig, mat.p, gene1)
  if (nrow(pairs_dact) == 0) {
    return(list(
      pairs_dact = pairs_dact,
      gene.pair = pairs_dact,
      sig_gene2 = character(0),
      non_sig.gene2 = character(0)
    ))
  }

  pairs_dact$wgs_gene2 <- p.wes[match(pairs_dact$gene2, names(p.wes))]

  # Map SNP exposures to cis genes using the provided SNP-to-gene table.
  snp.map <- uniq_snp[match(pairs_dact$gene1, uniq_snp$SNP), , drop = FALSE]
  pairs_dact$cis_gene1 <- snp.map$GeneSymbol

  # For SNPs not mapped by uniq_snp, optionally infer the closest or overlapping gene.
  not.mapped <- unique(pairs_dact$gene1[is.na(pairs_dact$cis_gene1) | pairs_dact$cis_gene1 == ""])
  if (length(not.mapped) > 0 && !is.null(SNP.ref)) {
    inferred.map <- infer_snp_to_gene(not.mapped, SNP.ref, ref.table.keep)
    if (nrow(inferred.map) > 0) {
      idx <- match(pairs_dact$gene1, inferred.map$SNP)
      fill <- is.na(pairs_dact$cis_gene1) | pairs_dact$cis_gene1 == ""
      pairs_dact$cis_gene1[fill & !is.na(idx)] <- inferred.map$GeneSymbol[idx[fill & !is.na(idx)]]
    }
  }

  colnames(pairs_dact) <- c("rsid", "gene2", "DANDELION_p", "wgs_gene2", "gene1")

  dandelion_message("Combining gene loci.", verbose = verbose)
  gene.pair <- pairs_dact[!is.na(pairs_dact$gene1) & pairs_dact$gene1 != "", , drop = FALSE]
  gene.pair <- gene.pair[!duplicated(gene.pair[, c("rsid", "gene2")]), , drop = FALSE]

  gene.pair$gene1 <- remove_gene_version(gene.pair$gene1)
  
  # Merge nearby exposure genes into locus-level regions for easier interpretation.
  gene.pair <- combine_gene_loci(gene.pair, ref.table.keep, gene1.col = "gene1")

  gene2.p <- gene.pair$wgs_gene2
  sig_gene2 <- unique(gene.pair$gene2[!is.na(gene2.p) & gene2.p < eta.wgs])
  non_sig.gene2 <- unique(gene.pair$gene2[!is.na(gene2.p) & gene2.p >= eta.wgs])

  out <- list(
    pairs_dact = pairs_dact,
    gene.pair = gene.pair,
    sig_gene2 = sig_gene2,
    non_sig.gene2 = non_sig.gene2
  )
  class(out) <- "dandelion_pairs"
  out
}


#' Organize DANDELION Gene-Gene Results
#'
#' @description
#' Cleans significant DANDELION pairs when the exposure side contains genes and
#' merges nearby exposure genes into genomic loci.
#'
#' @param mat.sig Matrix encoding significant pairs, returned by `med_gene()`.
#' @param mat.p Matrix of DANDELION p-values, returned by `med_gene()`.
#' @param p.wes Named numeric vector of gene-level trait association p-values.
#' @param gene1 Candidate exposure genes returned by `med_gene()`.
#' @param ref.table.keep Gene annotation data frame after filtering. Required
#'   columns are `gene_name`, `Chromosome`, `start`, and `end`.
#' @param eta.wgs Significance threshold for WES genes. Default is 1e-5.
#' @param verbose Logical. If `TRUE`, progress messages are printed. Default is `FALSE`.
#'
#' @return A list containing `pairs_dact`, `gene.pair`, `sig_gene2`, and
#'   `non_sig.gene2`.
#'
#' @export
calc_pair.gene <- function(mat.sig,
                           mat.p,
                           p.wes,
                           gene1,
                           ref.table.keep,
                           eta.wgs = 1e-5,
                           verbose = FALSE) {
  check_pair_inputs(mat.sig, mat.p, p.wes, gene1, ref.table.keep)

  mat.sig.new <- mat.sig[, gene1, drop = FALSE]
  mat.sig.new <- mat.sig.new[rowSums(mat.sig.new) != 0, colSums(mat.sig.new) != 0, drop = FALSE]

  dandelion_message("Removing self-connected genes.", verbose = verbose)
  # Remove self-connections where the exposure gene and target gene are identical.
  common.self <- intersect(rownames(mat.sig.new), colnames(mat.sig.new))
  if (length(common.self) > 0) {
    for (g in common.self) {
      mat.sig.new[g, g] <- 0L
    }
  }
  mat.sig.new <- mat.sig.new[rowSums(mat.sig.new) != 0, colSums(mat.sig.new) != 0, drop = FALSE]

  pairs_dact <- build_pairs_from_mats(mat.sig.new, mat.p, colnames(mat.sig.new))
  if (nrow(pairs_dact) == 0) {
    return(list(
      pairs_dact = pairs_dact,
      gene.pair = pairs_dact,
      sig_gene2 = character(0),
      non_sig.gene2 = character(0)
    ))
  }

  pairs_dact$wgs_gene2 <- p.wes[match(pairs_dact$gene2, names(p.wes))]

  dandelion_message("Combining gene loci.", verbose = verbose)
  gene.pair <- pairs_dact[!is.na(pairs_dact$gene1) & pairs_dact$gene1 != "", , drop = FALSE]
  gene.pair <- gene.pair[!duplicated(gene.pair[, c("gene1", "gene2")]), , drop = FALSE]
  gene.pair$gene1 <- remove_gene_version(gene.pair$gene1)
  gene.pair <- combine_gene_loci(gene.pair, ref.table.keep, gene1.col = "gene1")

  gene2.p <- gene.pair$wgs_gene2
  sig_gene2 <- unique(gene.pair$gene2[!is.na(gene2.p) & gene2.p < eta.wgs])
  non_sig.gene2 <- unique(gene.pair$gene2[!is.na(gene2.p) & gene2.p >= eta.wgs])

  out <- list(
    pairs_dact = pairs_dact,
    gene.pair = gene.pair,
    sig_gene2 = sig_gene2,
    non_sig.gene2 = non_sig.gene2
  )
  class(out) <- "dandelion_pairs"
  out
}


#' Print a DANDELION result object
#'
#' @param x A `dandelion_result` object returned by `med_gene()`.
#' @param ... Additional arguments, currently unused.
#'
#' @return Invisibly returns `x`.
#'
#' @method print dandelion_result
#' @export
print.dandelion_result <- function(x, ...) {
  cat("DANDELION result\n")
  cat("  Exposure type: ", x$gene1.type, "\n", sep = "")
  cat("  Number of analyzed exposures with valid results: ", length(x$gene1), "\n", sep = "")
  cat("  Significant pairs: ", sum(x$mat.sig != 0, na.rm = TRUE), "\n", sep = "")
  invisible(x)
}

#' Print DANDELION pair results
#'
#' @param x A `dandelion_pairs` object returned by `calc_pair.gene()` or
#'   `calc_pair.snp()`.
#' @param ... Additional arguments, currently unused.
#'
#' @return Invisibly returns `x`.
#'
#' @method print dandelion_pairs
#' @export
print.dandelion_pairs <- function(x, ...) {
  cat("DANDELION pair results\n")
  cat("  Number of significant pairs: ", nrow(x$pairs_dact), "\n", sep = "")
  cat("  Number of locus-level pairs: ", nrow(x$gene.pair), "\n", sep = "")
  cat("  Number of WES-significant gene2: ", length(x$sig_gene2), "\n", sep = "")
  cat("  Number of non-WES-significant gene2: ", length(x$non_sig.gene2), "\n", sep = "")
  invisible(x)
}

#' Generate DANDELION Network Figure
#'
#' @description
#' Generates a network plot from DANDELION gene-pair results.
#'
#' @param gene.pair Data frame of identified gene pairs, returned by
#'   `calc_pair.gene()` or `calc_pair.snp()`. Must contain columns `region` and
#'   `gene2`.
#' @param p.wes Named numeric vector of gene-level trait association p-values.
#' @param eta.wgs Significance threshold for WES genes. Default is 1e-5.
#' @param pic_dir Directory to save the figure.
#'
#' @return Invisibly returns the path to the saved PDF file.
#'
#' @export
gen_fig <- function(gene.pair, p.wes, eta.wgs = 1e-5, pic_dir) {
  if (!all(c("region", "gene2") %in% colnames(gene.pair))) {
    stop("gene.pair must contain columns: region and gene2.")
  }
  if (missing(pic_dir) || is.null(pic_dir)) {
    stop("pic_dir must be provided.")
  }

  if (!dir.exists(pic_dir)) {
    dir.create(pic_dir, recursive = TRUE, showWarnings = FALSE)
  }

  WES.imp <- names(p.wes)[p.wes <= eta.wgs]

  result.pair <- unique(gene.pair[, c("region", "gene2")])
  result.pair <- result.pair[!is.na(result.pair$region) & !is.na(result.pair$gene2), , drop = FALSE]

  if (nrow(result.pair) == 0) {
    warning("No valid pairs available for plotting.")
    return(invisible(NULL))
  }

  node.g1 <- unique(c(result.pair$region, result.pair$gene2))
  net.g1 <- igraph::graph_from_data_frame(
    d = result.pair[, c("region", "gene2")],
    vertices = data.frame(name = node.g1),
    directed = TRUE
  )

  igraph::V(net.g1)$type <- "Not"
  igraph::V(net.g1)$type[igraph::V(net.g1)$name %in% WES.imp] <- "Important"
  igraph::E(net.g1)$type <- "Not"

  node_sizes <- log(igraph::degree(net.g1) + 1) + 1

  out.file <- file.path(pic_dir, "Network.pdf")
  grDevices::pdf(file = out.file, width = 16, height = 12)
  on.exit(grDevices::dev.off(), add = TRUE)

  plot(
    net.g1,
    edge.arrow.size = 0.1,
    edge.curved = 0,
    vertex.size = node_sizes,
    vertex.color = "orange",
    vertex.frame.color = "#555555",
    edge.color = c("dark red", "slategrey")[(igraph::E(net.g1)$type == "Not") + 1],
    vertex.label.color = c("#bf484a", "#50649d")[(igraph::V(net.g1)$type == "Not") + 1],
    vertex.label.dist = 0.1,
    vertex.label.cex = 0.3
  )

  invisible(out.file)
}


# Internal helpers ---------------------------------------------------------
# Validate input objects and required columns before running DANDELION.
# This keeps downstream error messages clearer for users.
check_med_gene_inputs <- function(p.trans, p.wes, ref.table, gene1.list, gene1.type, SNP.ref) {
  if (!is.matrix(p.trans)) {
    stop("p.trans must be a matrix.")
  }
  if (is.null(rownames(p.trans)) || is.null(colnames(p.trans))) {
    stop("p.trans must have both row names and column names.")
  }
  if (!is.numeric(p.wes) || is.null(names(p.wes))) {
    stop("p.wes must be a named numeric vector.")
  }
  
  common_check <- intersect(rownames(p.trans), names(p.wes))
  if (length(common_check) == 0) {
    stop(
      "Data mismatch: no overlapping gene names were found between ",
      "rownames(p.trans) and names(p.wes). Please ensure both objects use ",
      "the same gene identifiers, such as gene symbols or Ensembl IDs."
    )
  }
  
  req.ref <- c("gene_name", "type", "Chromosome", "start", "end")
  if (!all(req.ref %in% colnames(ref.table))) {
    stop("ref.table must contain columns: gene_name, type, Chromosome, start, end.")
  }
  if (!is.character(gene1.list)) {
    stop("gene1.list must be a character vector.")
  }
  if (gene1.type == "SNP") {
    if (is.null(SNP.ref)) {
      stop("SNP.ref must be provided when gene1.type = 'SNP'.")
    }
    req.snp <- c("SNP", "SNPPos", "SNPChr")
    if (!all(req.snp %in% colnames(SNP.ref))) {
      stop("SNP.ref must contain columns: SNP, SNPPos, SNPChr.")
    }
  }
}

check_pair_inputs <- function(mat.sig, mat.p, p.wes, gene1, ref.table.keep) {
  if (!is.matrix(mat.sig) || !is.matrix(mat.p)) {
    stop("mat.sig and mat.p must be matrices.")
  }
  if (!identical(dim(mat.sig), dim(mat.p))) {
    stop("mat.sig and mat.p must have the same dimensions.")
  }
  if (!identical(rownames(mat.sig), rownames(mat.p)) || !identical(colnames(mat.sig), colnames(mat.p))) {
    stop("mat.sig and mat.p must have the same row and column names.")
  }
  if (!is.numeric(p.wes) || is.null(names(p.wes))) {
    stop("p.wes must be a named numeric vector.")
  }
  if (length(gene1) == 0) {
    stop("gene1 is empty.")
  }
  if (!all(gene1 %in% colnames(mat.sig))) {
    stop("All values in gene1 must be present in colnames(mat.sig).")
  }
  req.ref <- c("gene_name", "Chromosome", "start", "end")
  if (!all(req.ref %in% colnames(ref.table.keep))) {
    stop("ref.table.keep must contain columns: gene_name, Chromosome, start, end.")
  }
}

normalize_chr <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^chr", x), x, paste0("chr", x))
}

remove_gene_version <- function(x) {
  sub("\\..*$", "", as.character(x))
}

# Internal helper for optional progress messages.
dandelion_message <- function(..., verbose = FALSE) {
  if (isTRUE(verbose)) {
    message(...)
  }
}

clamp_p <- function(p) {
  p <- as.numeric(p)
  p[p <= 0] <- .Machine$double.xmin
  p[p >= 1] <- 1 - 1e-15
  p
}

find_cis_genes_for_region <- function(chr, start.pos, end.pos, ref.table, dist) {
  start.critical <- as.numeric(start.pos) - dist
  end.critical <- as.numeric(end.pos) + dist

  target.genes <- ref.table[ref.table$Chromosome == chr, , drop = FALSE]
  if (nrow(target.genes) == 0) {
    return(character(0))
  }

  loc.pos <- which(
    (target.genes$start > start.critical & target.genes$start < end.critical) |
      (target.genes$end > start.critical & target.genes$end < end.critical) |
      (target.genes$start <= start.critical & target.genes$end >= end.critical)
  )

  unique(target.genes$gene_name[loc.pos])
}

run_dandelion_for_exposure <- function(exposure.id,
                                       exposure.index,
                                       gene.trans,
                                       p.trans.new,
                                       p.wes.new,
                                       target.fdr) {
  if (length(gene.trans) == 0) {
    return(NULL)
  }

  p_a <- p.trans.new[gene.trans, exposure.index]
  p_b <- p.wes.new[gene.trans]

  # Remove genes with missing trans-association or burden-test p-values.
  valid <- !is.na(p_a) & !is.na(p_b)
  p_a <- p_a[valid]
  p_b <- p_b[valid]
  gene.trans <- gene.trans[valid]

  if (length(gene.trans) < 2) {
    return(NULL)
  }

  p_a <- clamp_p(p_a)
  p_b <- clamp_p(p_b)
  names(p_a) <- gene.trans
  names(p_b) <- gene.trans

  Z_a <- stats::qnorm(p_a, lower.tail = FALSE)
  Z_b <- stats::qnorm(p_b, lower.tail = FALSE)
  
  # Step: Calculate mixture weights for the DACT-style model.
  # These weights (wg1, wg2, wg3) are automatically estimated from the genome-wide 
  # distribution of p-values. This ensures that the model adaptively accounts for 
  # the signal strength in both trans-regulation and trait-association data.
  pi0a <- 1 - nonnullPropEst(Z_a, 0, 1) # Estimated proportion of null trans-eQTLs
  pi0b <- 1 - nonnullPropEst(Z_b, 0, 1) # Estimated proportion of null trait-associations

  if (any(is.na(c(pi0a, pi0b))) || pi0a < 0 || pi0b < 0) {
    return(NULL)
  }

  pi0a <- min(pi0a, 1)
  pi0b <- min(pi0b, 1)

  # Combine the component p-values using DACT-style mixture weights.
  # p3 corresponds to the joint null component based on the MaxP statistic.
  p3 <- (pmax(p_a, p_b))^2
  wg1 <- pi0a * (1 - pi0b)
  wg2 <- (1 - pi0a) * pi0b
  wg3 <- pi0a * pi0b
  wg.sum <- wg1 + wg2 + wg3

  if (is.na(wg.sum) || wg.sum <= 0) {
    return(NULL)
  }

  wg.std <- c(wg1, wg2, wg3) / wg.sum
  p_dact <- wg.std[1] * p_a + wg.std[2] * p_b + wg.std[3] * p3
  p_dact <- clamp_p(p_dact)
  names(p_dact) <- gene.trans


  # Convert DANDELION p-values to q-values when qvalue is available;
  # otherwise use Benjamini-Hochberg adjusted p-values.
  q_dact <- safe_qvalues(p_dact)
  names(q_dact) <- gene.trans

  # Store significant indicators and DANDELION p-values in vectors aligned
  # to the rows of p.trans.new.
  vec.sig <- rep(0L, nrow(p.trans.new))
  names(vec.sig) <- rownames(p.trans.new)
  sig_dact <- names(q_dact)[q_dact <= target.fdr]
  vec.sig[sig_dact] <- 1L

  vec.p <- rep(NA_real_, nrow(p.trans.new))
  names(vec.p) <- rownames(p.trans.new)
  vec.p[gene.trans] <- p_dact

  list(
    id = exposure.id,
    sig.vec = unname(vec.sig),
    p.vec = unname(vec.p)
  )
}

build_pairs_from_mats <- function(mat.sig, mat.p, gene1) {
  mat.sig.new <- mat.sig[, gene1, drop = FALSE]
  if (nrow(mat.sig.new) == 0 || ncol(mat.sig.new) == 0) {
    return(data.frame(
      gene1 = character(),
      gene2 = character(),
      DANDELION_p = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  mat.sig.new <- mat.sig.new[rowSums(mat.sig.new) != 0, colSums(mat.sig.new) != 0, drop = FALSE]
  if (nrow(mat.sig.new) == 0 || ncol(mat.sig.new) == 0) {
    return(data.frame(
      gene1 = character(),
      gene2 = character(),
      DANDELION_p = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  loc <- which(mat.sig.new != 0, arr.ind = TRUE)
  data.frame(
    gene1 = colnames(mat.sig.new)[loc[, "col"]],
    gene2 = rownames(mat.sig.new)[loc[, "row"]],
    DANDELION_p = as.numeric(mat.p[cbind(rownames(mat.sig.new)[loc[, "row"]], colnames(mat.sig.new)[loc[, "col"]])]),
    stringsAsFactors = FALSE
  )
}

infer_snp_to_gene <- function(snps, SNP.ref, ref.table.keep) {
  req.snp <- c("SNP", "SNPPos", "SNPChr")
  if (!all(req.snp %in% colnames(SNP.ref))) {
    stop("SNP.ref must contain columns: SNP, SNPPos, SNPChr.")
  }

  out <- data.frame(SNP = character(), GeneSymbol = character(), stringsAsFactors = FALSE)
  for (snp in snps) {
    snp.info <- SNP.ref[SNP.ref$SNP == snp, , drop = FALSE]
    if (nrow(snp.info) == 0) next
    snp.info <- snp.info[1, , drop = FALSE]

    snp.chr <- normalize_chr(snp.info$SNPChr)
    snp.pos <- as.numeric(snp.info$SNPPos)
    if (is.na(snp.pos)) next

    gene.pos <- ref.table.keep[ref.table.keep$Chromosome == snp.chr, , drop = FALSE]
    if (nrow(gene.pos) == 0) next

    inside <- which(gene.pos$start <= snp.pos & gene.pos$end >= snp.pos)
    if (length(inside) >= 1) {
      mapped.gene <- gene.pos$gene_name[inside[1]]
    } else {
      dist.to.gene <- pmin(abs(gene.pos$start - snp.pos), abs(gene.pos$end - snp.pos))
      mapped.gene <- gene.pos$gene_name[which.min(dist.to.gene)]
    }

    out <- rbind(out, data.frame(SNP = snp, GeneSymbol = mapped.gene, stringsAsFactors = FALSE))
  }
  out
}

# Internal helper:
# The function merges identified exposure genes into "genomic loci" based on 
# physical proximity (default 500kb window). This simplifies the results by 
# grouping nearby genes that likely represent the same regulatory signal, 
# making the output more interpretable for biological follow-up.

combine_gene_loci <- function(gene.pair, ref.table.keep, gene1.col = "gene1", window.size = 5e5) {
  if (nrow(gene.pair) == 0) {
    gene.pair$region <- character(0)
    return(gene.pair)
  }

  gene.pair$region <- gene.pair[[gene1.col]]
  gene1_list <- unique(gene.pair[[gene1.col]])
  ref.subset <- ref.table.keep[ref.table.keep$gene_name %in% gene1_list, , drop = FALSE]

  if (nrow(ref.subset) == 0) {
    return(gene.pair)
  }

  for (chr in unique(ref.subset$Chromosome)) {
    ref.chr <- ref.subset[ref.subset$Chromosome == chr, , drop = FALSE]
    ref.chr <- ref.chr[order(ref.chr$start), , drop = FALSE]

    if (nrow(ref.chr) == 0) next

    region.id <- 1L
    current.start <- ref.chr$start[1]
    current.end <- ref.chr$end[1]
    region.list <- vector("list", 0)
    current.genes <- ref.chr$gene_name[1]

    if (nrow(ref.chr) >= 2) {
      for (i in 2:nrow(ref.chr)) {
        next.start <- ref.chr$start[i]
        next.end <- ref.chr$end[i]
        next.gene <- ref.chr$gene_name[i]

        if (next.start <= current.end + window.size) {
          current.end <- max(current.end, next.end)
          current.genes <- c(current.genes, next.gene)
        } else {
          region.list[[region.id]] <- current.genes
          region.id <- region.id + 1L
          current.start <- next.start
          current.end <- next.end
          current.genes <- next.gene
        }
      }
    }
    region.list[[region.id]] <- current.genes

    for (genes.in.region in region.list) {
      region.name <- paste(unique(genes.in.region), collapse = ", ")
      idx <- gene.pair[[gene1.col]] %in% genes.in.region
      gene.pair$region[idx] <- region.name
    }
  }

  gene.pair
}

#' Estimate the Non-null Proportion
#'
#' @description Internal helper used by DANDELION.
#'
#' @param x Numeric vector.
#' @param u Null mean.
#' @param sigma Null standard deviation.
#'
#' @return Estimated non-null proportion.
#'
#' @keywords internal
nonnullPropEst <- function(x, u, sigma) {
  x <- x[is.finite(x)]
  if (length(x) < 2) {
    return(NA_real_)
  }

  z <- (x - u) / sigma
  xi <- seq(0, 1, by = 0.01)
  tmax <- sqrt(log(length(x)))
  tt <- seq(0, tmax, by = 0.1)

  epsest <- numeric(length(tt))
  for (j in seq_along(tt)) {
    t <- tt[j]
    f <- exp((t * xi)^2 / 2)
    w <- 1 - abs(xi)
    co <- numeric(length(xi))
    for (i in seq_along(xi)) {
      co[i] <- mean(cos(t * xi[i] * z))
    }
    epsest[j] <- 1 - sum(w * f * co) / sum(w)
  }

  max(epsest, na.rm = TRUE)
}

safe_qvalues <- function(p) {
  p <- clamp_p(p)
  
  if (length(p) < 10 || length(unique(p)) < 4) {
    return(stats::p.adjust(p, method = "BH"))
  }
  
  if (requireNamespace("qvalue", quietly = TRUE)) {
    return(tryCatch(
      qvalue::qvalue(p)$qvalues,
      error = function(e) stats::p.adjust(p, method = "BH"),
      warning = function(w) stats::p.adjust(p, method = "BH")
    ))
  }
  
  stats::p.adjust(p, method = "BH")
}