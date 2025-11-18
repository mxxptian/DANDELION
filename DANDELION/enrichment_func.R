run_enrichment_all <- function(
    UK_dir,
    pheno_file,
    mendelian_file,
    gene_position_file,
    n_perm = 1000,
    idx_start = 1,
    idx_end = NULL
){
  library(dplyr)
  
  # --- load phenotype list ---
  pheno_list <- read.csv(pheno_file)
  file_list <- list.files(paste0(UK_dir, "whole_blood"))
  info.out <- do.call(rbind, strsplit(as.character(file_list), "_"))
  
  if(is.null(idx_end)) idx_end <- nrow(info.out)
  
  enrich_data <- NULL
  
  for(p in idx_start:idx_end){
    cat("Running phenotype:", info.out[p,1], "\n")
    
    temp <- run_enrichment_single(
      p = p,
      UK_dir = UK_dir,
      pheno_list = pheno_list,
      mendelian_file = mendelian_file,
      gene_position_file = gene_position_file,
      info.out = info.out,
      n_perm = n_perm
    )
    
    enrich_data <- rbind(enrich_data, temp)
  }
  
  return(enrich_data)
}


run_enrichment_single <- function(
    p,
    UK_dir,
    pheno_list,
    mendelian_file,
    gene_position_file,
    info.out,
    n_perm = 1000
){
  library(dplyr)
  
  trait_name <- info.out[p,1]
  trait_dir  <- paste0(UK_dir, trait_name, "/")
  
  # --- Load TWAS file ---------------------------
  load(paste0(trait_dir, trait_name, ".Rdata"))
  
  # Process gene names
  info.out.pheno <- do.call(rbind, strsplit(as.character(dat_M3.1$Name),'\\.'))
  info.out.new   <- do.call(rbind, strsplit(as.character(info.out.pheno[,1]), "\\("))
  dat_M3.1$Name  <- info.out.new[,1]
  dat_M3.1_new   <- dat_M3.1[!dat_M3.1$Name %in% names(which(table(dat_M3.1$Name) != 1)), ]
  
  # p-values
  cau <- data.frame(pval = dat_M3.1_new$p_value)
  rownames(cau) <- dat_M3.1_new$Name
  p.wgs <- cau$pval
  names(p.wgs) <- rownames(cau)
  
  threshold <- 0.05 / nrow(cau)
  
  # --- Load DACT pair file -----------------------
  load(paste0(trait_dir, trait_name, "_pair.Rdata"))
  
  gene2 <- unique(result.pair$pairs_dact$gene2)
  sig_gene2    <- unique(result.pair$pairs_dact$gene2[result.pair$pairs_dact$wgs_gene2 < threshold])
  non_sig_gene2<- unique(result.pair$pairs_dact$gene2[result.pair$pairs_dact$wgs_gene2 >= threshold])
  
  # --- Mendelian gene list -----------------------
  mendelian <- read.csv(mendelian_file)
  mend.gene <- unique(mendelian$Gene_Symbol_HGNC)
  
  # --- Gene position table -----------------------
  ref.table <- read.delim(gene_position_file, sep="\t")
  ref.table$length <- ref.table$end - ref.table$start
  ref.table <- ref.table[!duplicated(ref.table$gene_name),]
  
  dist.twas <- data.frame(
    gene   = ref.table$gene_name,
    length = ref.table$length,
    start  = ref.table$length - 100000,
    end    = ref.table$length + 100000
  )
  
  cau$gene <- rownames(cau)
  cau_len <- right_join(dist.twas, cau, by="gene")
  
  # =====================================================
  # Helper function: permutation overlap
  # =====================================================
  perm_overlap <- function(target_genes){
    overlap_A <- intersect(target_genes, mend.gene)
    B_list <- numeric(n_perm)
    
    if(length(overlap_A)==0) return(rep(NA, n_perm))
    
    exclude_genes <- cau_len[!cau_len$gene %in% target_genes, ]
    
    for(i in 1:n_perm){
      choose_B <- c()
      for(g in target_genes){
        pos <- which(cau_len$gene == g)
        start <- cau_len$start[pos]
        end   <- cau_len$end[pos]
        
        idx <- which(exclude_genes$length > start &
                       exclude_genes$length < end)
        
        sim.genes <- exclude_genes$gene[idx]
        if(length(sim.genes) == 0){
          choose_B <- c(choose_B, NA)
        }else{
          choose_B <- c(choose_B, sample(sim.genes, 1))
        }
      }
      
      B_list[i] <- sum(choose_B %in% mend.gene, na.rm=TRUE)
    }
    
    return(list(
      overlapA = length(overlap_A),
      permB    = B_list
    ))
  }
  
  # =====================================================
  # Run for different categories
  # =====================================================
  res_gene2      <- perm_overlap(gene2)
  res_sig_gene2  <- perm_overlap(sig_gene2)
  res_nonsig     <- perm_overlap(non_sig_gene2)
  
  imp <- names(p.wgs)[which(p.wgs <= threshold)]
  burden_gene2 <- imp[!imp %in% gene2]
  res_burden <- perm_overlap(burden_gene2)
  
  # =====================================================
  # Construct return table
  # =====================================================
  out <- data.frame(
    p_overlap           = mean(res_gene2$permB >= res_gene2$overlapA, na.rm=TRUE),
    p_sig_overlap       = mean(res_sig_gene2$permB >= res_sig_gene2$overlapA, na.rm=TRUE),
    p_nonsig_overlap    = mean(res_nonsig$permB >= res_nonsig$overlapA, na.rm=TRUE),
    p_burden_overlap    = mean(res_burden$permB >= res_burden$overlapA, na.rm=TRUE),
    
    avg_overlap         = mean(res_gene2$permB, na.rm=TRUE),
    avg_sig_overlap     = mean(res_sig_gene2$permB, na.rm=TRUE),
    avg_nonsig_overlap  = mean(res_nonsig$permB, na.rm=TRUE),
    avg_burden_overlap  = mean(res_burden$permB, na.rm=TRUE),
    
    num_overlap_A       = res_gene2$overlapA,
    num_sig_overlap_A   = res_sig_gene2$overlapA,
    num_nonsig_overlap_A= res_nonsig$overlapA,
    num_burden_overlap_A= res_burden$overlapA,
    
    trait = pheno_list[pheno_list$Study_ID == trait_name, "Abbr"]
  )
  
  return(out)
}


# Run an example
# result <- run_enrichment_all(
# UK_dir = "/Users/px/Desktop/Uchicago/Mediation/UK Biobank/",
# pheno_file = "/Users/px/Desktop/Uchicago/Mediation/UK Biobank/phenotypes.csv",
# mendelian_file = "/Users/px/Desktop/Uchicago/Mediation/UK Biobank/Mendelian genes of disease disorders.csv",
# gene_position_file = "/Users/px/Desktop/Uchicago/Mediation/Testosterone/gene_position.txt",
# n_perm = 1000,
#  idx_start = 1,
#  idx_end = 3
# )

# save(result, file = "enrichment_replicate1000.Rdata")
