#' @title Apply DANDELION to identify candidate disease proximal genes.
#' @description Apply DANDELION using p-values from WES and eQTLGen data to generate DANDELION's p-values and identified pairs. For code clarity, disease distal genes are abbreviated as 'gene1', and disease proximal genes as 'gene2' within the function.
#' @param p.trans A matrix of p-values for distal-> proximal (gene1 -> gene2) associations.
#'        The matrix has dimensions M (disease proximal genes, gene2) x N (disease distal genes, gene1).
#'        Each entry represents the p-value for a potential trans effect from a gene1 to a gene2.
#'        Both row names (gene2) and column names (gene1) should be gene symbols.
#' @param p.wes A named numeric vector containing p-values for associations between disease proximal genes (gene2) and the trait.
#'        The names should correspond to core gene identifiers (gene2), and must match the row names of p.trans.
#' @param gene1.list A character vector of candidate disease distal gene names (gene1),
#'        which must be a subset of the column names of p.trans.
#' @param ref.table reference data for position information.
#         should contain five columns:
#         1. gene_name: gene identifier, which can be a gene symbol, Ensembl ID,
#         2. type: gene type
#         3. Chromosome: 'chr1'
#         4. start: start position of gene
#         5. end: end position of gene
#' @param target.fdr the pre-defined FDR level in the BH step, the default value is 0.1.
#' @param dist the distant between the gene and trans-genes, the default value is 5e6
#' @param gene1.type specify the exposure type. The expousre can be SNP, or Gene.
#' @param SNP.ref reference data for position information if the exposure is 'SNP'
#'        should contain three columns:
#'        1. SNP: rsid
#'        2. SNPPos: position of SNP
#'        3. SNPChr: '1'
#'
#' @return A list of candidate genes 1, identified pair matrix, and DANDELION's p-values
#' @export



med_gene <- function(p.trans, p.wes, ref.table, gene1.list, target.fdr=0.1, dist=5e6, gene1.type = c('SNP','Gene'), SNP.ref = NULL){


  # normalize type
  gene1.type <- match.arg(gene1.type)
  gene1.type.lc <- tolower(gene1.type)
  ##### Step 0: sample the candidate gene. #####

  p.trans = p.trans[,gene1.list]

  ##### Step 1: match the common gene. #####

  cat('\n-----Start matching gene 2-----\n')

  ## match gene 2 of p.trans and p.wes
  common_gene2 <- intersect(rownames(p.trans), names(p.wes))

  ## keep only lincRNA and protein coding genes
  ref.table.keep <- ref.table[ref.table$type %in% c('lincRNA', 'protein_coding'),]

  ## remove genes on chr M, X, and Y
  ref.table.keep <- ref.table.keep[!(ref.table.keep$Chromosome %in% c('chrM', 'chrX', 'chrY')),]


  ref.table.keep = ref.table.keep[!duplicated(ref.table.keep$gene_name),]


  ## match common_gene2 to reference table
  common_gene2 <- intersect(common_gene2, ref.table.keep$gene_name)

  ## keep common gene2 for p.trans
  p.trans.new <- p.trans[rownames(p.trans)%in%common_gene2,]

  ## keep common gene2 for p.wes
  p.wes.new <- p.wes[names(p.wes)%in%common_gene2]

  ## reorder p.trans.new and p.wes.new
  p.trans.new <- p.trans.new[match(names(p.wes.new), rownames(p.trans.new)),]

  cat('\n-----End matching gene 2-----\n')

  ##### Step 2: remove the cis genes. #####

  mat.sig <- matrix(0, nrow = nrow(p.trans.new), ncol = ncol(p.trans.new))
  rownames(mat.sig) <- rownames(p.trans.new)
  colnames(mat.sig) <- colnames(p.trans.new)
  # matrix that encodes whether the gene1-gene2 pair are significant

  mat.p <- matrix(NA, nrow = nrow(p.trans.new), ncol = ncol(p.trans.new))
  rownames(mat.p) <- rownames(p.trans.new)
  colnames(mat.p) <- colnames(p.trans.new)
  # matrix that stores p-values of DACT results

  gene1.candidate <- c()  # list of gene 1 included in analysis

  M <- ncol(p.trans.new) # number of candidate gene 1

  ##### Step 3: worker for Gene #####
  run_single_gene <- function(i) {
    gene1 <- colnames(p.trans.new)[i]
    # only if gene1 exists in ref.table.keep
    if (gene1 %in% ref.table.keep$gene_name){

      pos.info  <- ref.table.keep[ref.table.keep$gene_name == gene1, ]
      start.pos <- pos.info$start
      end.pos   <- pos.info$end
      chr       <- pos.info$Chromosome

      # define cis-window
      start.critical <- start.pos - dist
      end.critical   <- end.pos   + dist

      # genes on same chr
      target.genes <- ref.table[ref.table$Chromosome == chr, ]
      loc.pos.1 <- which(target.genes$start > start.critical & target.genes$start < end.critical)
      loc.pos.2 <- which(target.genes$end   > start.critical & target.genes$end   < end.critical)
      loc.pos    <- union(loc.pos.1, loc.pos.2)
      gene.cis   <- target.genes$gene_name[loc.pos]

      gene.trans <- setdiff(common_gene2, gene.cis)
      if (length(gene.trans) == 0) return(NULL)

      # DACT step
      p_a <- p.trans.new[gene.trans, i]; names(p_a) <- gene.trans
      p_b <- p.wes.new[gene.trans];     p_b[p_b == 1] <- 0.99

      Z_a  <- stats::qnorm(p_a, lower.tail = FALSE)
      Z_b  <- stats::qnorm(p_b, lower.tail = FALSE)
      pi0a <- 1 - nonnullPropEst(Z_a, 0, 1)
      pi0b <- 1 - nonnullPropEst(Z_b, 0, 1)

      if (any(is.na(c(pi0a, pi0b))) || pi0a < 0 || pi0b < 0) return(NULL)
      pi0a <- min(pi0a, 1); pi0b <- min(pi0b, 1)

      p3    <- (pmax(p_a, p_b))^2
      wg1   <- pi0a * (1 - pi0b)
      wg2   <- (1 - pi0a) * pi0b
      wg3   <- pi0a *  pi0b
      wg.sum<- wg1 + wg2 + wg3
      wg.std<- c(wg1, wg2, wg3) / wg.sum

      p_dact <- wg.std[1]*p_a + wg.std[2]*p_b + wg.std[3]*p3

      s <- try(qvalue::qvalue(p_dact), silent = TRUE)
      if (inherits(s, "try-error")) s <- qvalue::qvalue(p_dact, pi0 = 1)
      q_dact <- s$qvalues; names(q_dact) <- gene.trans

      vec.sig <- rep(0L, nrow(p.trans.new))
      sig_dact <- names(q_dact)[q_dact <= target.fdr]
      vec.sig[rownames(p.trans.new) %in% sig_dact] <- 1L

      vec.p <- rep(NA_real_, nrow(p.trans.new))
      vec.p[match(gene.trans, rownames(p.trans.new))] <- p_dact

      list(id = gene1, sig.vec = vec.sig, p.vec = vec.p)
    }
  }
  ##### Step 4: worker for SNP #####
  run_single_snp <- function(i) {
    cand.snp <- colnames(p.trans.new)[i]

    if (is.null(SNP.ref)) return(NULL)
    req_cols <- c("SNP", "SNPPos", "SNPChr")
    if (!all(req_cols %in% colnames(SNP.ref))) {
      stop("SNP.ref must contain columns: SNP, SNPPos, SNPChr")
    }
    if(cand.snp %in% SNP.ref$SNP){

      pos.info  <- SNP.ref[SNP.ref$SNP == cand.snp, ]
      start.pos <- unique(pos.info$SNPPos)
      end.pos   <- unique(pos.info$SNPPos)
      chr       <- paste0("chr", unique(pos.info$SNPChr))

      start.critical <- start.pos - dist
      end.critical   <- end.pos   + dist

      target.genes <- ref.table[ref.table$Chromosome == chr, ]
      loc.pos.1 <- which(target.genes$start > start.critical & target.genes$start < end.critical)
      loc.pos.2 <- which(target.genes$end   > start.critical & target.genes$end   < end.critical)
      loc.pos    <- union(loc.pos.1, loc.pos.2)
      gene.cis   <- target.genes$gene_name[loc.pos]

      gene.trans <- setdiff(common_gene2, gene.cis)
      if (length(gene.trans) == 0) return(NULL)

      p_a <- p.trans.new[gene.trans, i]; names(p_a) <- gene.trans
      p_b <- p.wes.new[gene.trans]
      p_a[p_a == 1] <- 0.99
      p_b[p_b == 1] <- 0.99

      Z_a  <- stats::qnorm(p_a, lower.tail = FALSE)
      Z_b  <- stats::qnorm(p_b, lower.tail = FALSE)
      pi0a <- 1 - nonnullPropEst(Z_a, 0, 1)
      pi0b <- 1 - nonnullPropEst(Z_b, 0, 1)

      if (any(is.na(c(pi0a, pi0b))) || pi0a < 0 || pi0b < 0) return(NULL)
      pi0a <- min(pi0a, 1); pi0b <- min(pi0b, 1)

      p3    <- (pmax(p_a, p_b))^2
      wg1   <- pi0a * (1 - pi0b)
      wg2   <- (1 - pi0a) * pi0b
      wg3   <- pi0a *  pi0b
      wg.sum<- wg1 + wg2 + wg3
      wg.std<- c(wg1, wg2, wg3) / wg.sum

      p_dact <- wg.std[1]*p_a + wg.std[2]*p_b + wg.std[3]*p3

      s <- try(qvalue::qvalue(p_dact), silent = TRUE)
      if (inherits(s, "try-error")) s <- qvalue::qvalue(p_dact, pi0 = 1)
      q_dact <- s$qvalues; names(q_dact) <- gene.trans

      vec.sig <- rep(0L, nrow(p.trans.new))
      sig_dact <- names(q_dact)[q_dact <= target.fdr]
      vec.sig[rownames(p.trans.new) %in% sig_dact] <- 1L

      vec.p <- rep(NA_real_, nrow(p.trans.new))
      vec.p[match(gene.trans, rownames(p.trans.new))] <- p_dact

      list(id = cand.snp, sig.vec = vec.sig, p.vec = vec.p)
    }
  }

  ##### Step 5: run sequentially #####
  if (gene1.type == "gene") {
    cat(paste0('\n-----Start DACT for ', M, ' genes (single-core) -----\n'))
    results <- lapply(seq_len(M), run_single_gene)
  } else {
    cat(paste0('\n-----Start DACT for ', M, ' SNPs (single-core) -----\n'))
    results <- lapply(seq_len(M), run_single_snp)
  }

  # collect
  for (res in results) {
    if (!is.null(res)) {
      gene1.candidate <- c(gene1.candidate, res$id)
      mat.sig[, res$id] <- res$sig.vec
      mat.p[,   res$id] <- res$p.vec
    }
  }


  if (gene1.type.lc == "gene") {
    cat(paste0('\nNumber of Gene1 detected: ', length(gene1.candidate), '\n'))
  } else {
    cat(paste0('\nNumber of SNPs detected: ', length(gene1.candidate), '\n'))
  }



  return(list(gene1=gene1.candidate, mat.sig=mat.sig, mat.p=mat.p))
}


#' @title Cleanse the pairs and disease proximal genes identified by DANDELION and map the disease proximal genes to WES significant genes (FOR genes 1 is SNPs)
#' @description Organize the output results of med_gene() and merge nearby genes into a single locus.
#' @param mat.sig matrix encoding significant pair, output by med_gene().
#' @param mat.p matrix encoding DANDELION P-values, output by med_gene()
#' @param gene1 candidate SNPs names used for analysis, output by med_gene()
#' @param uniq_snp matrix should contains columns about the snp (named as 'SNP') and the corresponding significant cis gene name (named as 'GeneSymbol')
#' @param p.wes a vector of length K. p value for gene 2 -> trait. Both p.trans and p.wes should contain the gene name.
#' @param ref.table.keep reference data for position information for 'lincRNA' and 'protein_coding' which are not located in 'chrX', 'chrY' and 'chrM'.
#' @param eta.wgs the significant threshold for WES gene
#'
#' @return a list includes identified pairs, combined pairs (loci based), identified disease proximal genes (also WES significant), and identified disease proximal genes (NOT WES significant)
#' @export



calc_pair.snp <- function(mat.sig, mat.p, p.wes, gene1, uniq_snp, ref.table.keep, eta.wgs=1e-5, GRCh = '37'){

  # mat.sig: matrix encoding significant pair, output by med_gene()
  # mat.p: matrix encoding DANDELION P-values, output by med_gene()
  # gene1: gene 1 names used for analysis, output by med_gene()
  # uniq_snp: matrix should contains columns about the snp (named as 'SNP') and the corresponding significant cis gene name (named as 'GeneSymbol')


  ensembl <- biomaRt::useEnsembl("snp", dataset = "hsapiens_snp", GRCh = GRCh) #keep consistent with the build version of gene position file

  mat.sig.new <- mat.sig[,gene1] # only consider gene 1 included in analysis


  mat.sig.new <- mat.sig.new[rowSums(mat.sig.new)!=0,colSums(mat.sig.new)!=0]

  pairs_dact <- data.frame(gene1=character(), gene2=character(), DENDELION_p=character(), stringsAsFactors = FALSE) # significant pairs by DACT
  for (i in 1:nrow(mat.sig.new)) {
    for (j in 1:ncol(mat.sig.new)) {
      if(mat.sig.new[i,j]!=0){
        pairs_dact[nrow(pairs_dact)+1,] <- c(colnames(mat.sig.new)[j],
                                             rownames(mat.sig.new)[i],
                                             mat.p[rownames(mat.sig.new)[i],colnames(mat.sig.new)[j]])
      }
    }
  }

  pairs_dact$wgs_gene2 <- p.wes[match(pairs_dact$gene2, names(p.wes))]

  ####### check the genome build of gene position file

  gene_test = ref.table.keep$gene_name[1]

  grch_gene = get_gene(gene_test)

  gene_sta = grch_gene$gene_start

  if(gene_sta==ref.table.keep[ref.table.keep$gene_name==gene_test,'start']){
    cat('\nGene positions are from GRCH 38!\n')
    gene_grch = 38
    if(GRCh!=gene_grch){

      cat('\nWarning: The SNP and gene positions are from different genome builds!!! Please correct it!!!\n')
    }else{
      cat('\nThe SNP and gene positions are from the same genome builds\n')
    }
  }else{
    cat('\nGene positions are from GRCH 37!\n')
    gene_grch = 37
    if(GRCh!=gene_grch){

      cat('\nWarning: The SNP and gene positions are from different genome builds!!! Please correct it!!!\n')
    }else{
      cat('\nThe SNP and gene positions are from the same genome builds.\n')
    }
  }

  ########################## match snp with gene ##############################

  snp.1 = data.frame(SNP = pairs_dact[,1])


  snp.merge = uniq_snp[match(snp.1$SNP, uniq_snp$SNP),]

  snp2cis = snp.merge$GeneSymbol

  pairs_dact$cis_gene1 = snp2cis



  ####### if there is no signficiant cis gene for the snps


  if(length(non_cis)==0){
    cat('\nAll the candidate genes 1 are significant cis genes!\n')
  }else{

    #get genomic position

    cat('\nThere are some insignificant cis genes 1\n!')
    SNPs <- getBM(attributes=c("refsnp_id",
                               "chr_name",
                               "chrom_start",
                               "chrom_end"),
                  filters ="snp_filter", values =non_cis, mart = ensembl, uniqueRows=TRUE)



    SNP.uniq = SNPs[!is.na(as.numeric(SNPs$chr_name)),]




    for (i in 1:nrow(SNP.uniq)) {
      rs.map = SNP.uniq[i,]

      gene.pos = ref.table.keep[ref.table.keep$Chromosome==paste0('chr',rs.map$chr_name),]

      loc.gene1 = which(gene.pos$start<rs.map$chrom_start&gene.pos$end>rs.map$chrom_start)

      # print(length(loc.gene1))

      if(length(loc.gene1)==1){
        idx = which( pairs_dact$gene1==rs.map$refsnp_id)

        pairs_dact$cis_gene1[idx] = gene.pos[loc.gene1,'gene_name']
      }else{
        gene.pos$start.diff = abs(gene.pos$start-rs.map$chrom_start)
        gene.pos$end.diff = abs(gene.pos$end-rs.map$chrom_start)
        if(min(gene.pos$start.diff)< min(gene.pos$end.diff)){

          idx = which( pairs_dact$gene1==rs.map$refsnp_id)

          idx.gene = which(gene.pos$start.diff==min(gene.pos$start.diff))

          pairs_dact$cis_gene1[idx] = gene.pos[idx.gene,'gene_name']
        }else{
          idx = which( pairs_dact$gene1==rs.map$refsnp_id)

          idx.gene = which(gene.pos$end.diff==min(gene.pos$end.diff))

          pairs_dact$cis_gene1[idx] = gene.pos[idx.gene,'gene_name']
        }
      }

    }
  }


  colnames( pairs_dact) = c('rsid', "gene2",  "DENDELION_p", 'wgs_gene2',  "gene1")



  ####################### combine into one loci ###################################

  # pair = pairs_dact[,c('gene1','gene2')]



  #remove gene1 is NA (there is no corresponding significant cis-gene mapped to the candidate snp)
  cat(paste0('\nCombineing gene loci\n'))
  gene.pair = pairs_dact[which(is.na(pairs_dact$gene1)!=TRUE),]

  gene.pair = gene.pair[!duplicated(gene.pair[,c(1,2)]),]


  gene1.out = strsplit(as.character(gene.pair$gene1), '\\.')

  gene1.info.out = do.call(rbind, gene1.out)
  gene.pair$gene1 = gene1.info.out[,1]


  gene.pair$region = gene.pair$gene1



  gene1_list = gene.pair$gene1

  ref.subset = ref.table.keep[ref.table.keep$gene_name%in%gene1_list,]

  dup.chr = names(table(ref.subset$Chromosome))[which(table(ref.subset$Chromosome)>=2)]

  ref.subset.1 = ref.subset[ref.subset$Chromosome%in%dup.chr,]

  if(nrow(ref.subset.1)!=0){

    for (chr in 1:length(unique(ref.subset.1$Chromosome))) {
      ref.subset.2 = ref.subset.1[ref.subset.1$Chromosome%in%unique(ref.subset.1$Chromosome)[chr],]

      ref.subset.2$region = NA
      ref.subset.2 = ref.subset.2[order(ref.subset.2$start),]

      window = 1

      g1 = 1

      while (g1<= nrow(ref.subset.2)) {

        start = ref.subset.2$start[g1]-500000
        end = ref.subset.2$end[g1]+500000

        loc.pos.1 <- which(ref.subset.2$start>start & ref.subset.2$end<end)

        ref.subset.2$region[loc.pos.1] = window

        g.end = max(loc.pos.1)

        if(g.end == g1){
          g1 = g1+1
          window = window + 1
        }else{
          g1 = g.end
        }


      }

      region_list = unique(ref.subset.2$region)

      for (r1 in 1:length(region_list)) {
        gene.list <- ref.subset.2$gene_name[which(ref.subset.2$region==r1)]

        gene.sub = gene.pair[gene.pair$gene1%in%gene.list,]



        gene.sub = gene.sub%>%distinct()
        gene.sub = gene.sub[order(gene.sub$gene1),]

        # gene.pair$comb_trans_p[(gene.pair$gene1%in%gene.sub$gene1)] = paste(gene.sub$trans_p, collapse = ", ")

        gene.pair$region[(gene.pair$gene1%in%gene.sub$gene1)] = paste(unique(gene.sub$gene1), collapse = ", ")
      }

    }
  }


  gene2 = gene.pair[,"wgs_gene2"]

  imp.var = names(p.wes)[p.wes<=threshold]


  sig_gene2 = gene.pair[which(gene2<threshold),"gene2"]
  non_sig.gene2 = gene.pair[which(gene2>=threshold),"gene2"]

  sig_gene2 = unique(sig_gene2)
  non_sig.gene2 = unique(non_sig.gene2)



  return(list(pairs_dact=pairs_dact,
              gene.pair = gene.pair,
              sig_gene2 = sig_gene2,
              non_sig.gene2 = non_sig.gene2))
}


#' @title Cleanse the pairs and disease proximal genes identified by DENDELION and map the disease proximal genes to WES significant genes (FOR genes 1 is gene.)
#' @description Organize the output results of med_gene() and merge nearby genes into a single locus.
#' @param mat.sig matrix encoding significant pair, output by med_gene().
#' @param mat.p matrix encoding DENDELION P-values, output by med_gene()
#' @param gene1 gene 1 names used for analysis, output by med_gene()
#' @param p.wes a vector of length K. p value for gene 2 -> trait. Both p.trans and p.wes should contain the gene name.
#' @param ref.table.keep reference data for position information for 'lincRNA' and 'protein_coding' which are not located in 'chrX', 'chrY' and 'chrM'.
#' @param eta.wgs the significant threshold for WES gene
#'
#' @return a list includes identified pairs, combined pairs (loci based), identified disease proximal genes (also WES significant), and identified disease proximal genes (NOT WES significant)
#' @export



calc_pair.gene <- function(mat.sig, mat.p, p.wes, gene1, ref.table.keep, eta.wgs=1e-5){

  # mat.sig: matrix encoding significant pair, output by med_gene()
  # mat.p: matrix encoding DENDELION P-values, output by med_gene()
  # gene1: gene 1 names used for analysis, output by med_gene()
  # uniq_snp: matrix should contains columns about the snp (named as 'SNP') and the corresponding significant cis gene name (named as 'GeneSymbol')




  mat.sig.new <- mat.sig[,gene1] # only consider gene 1 included in analysis

  mat.sig.new <- mat.sig.new[rowSums(mat.sig.new)!=0,colSums(mat.sig.new)!=0]

  cat(paste0('\nRemoving self-connected genes\n'))

  for (i in 1:nrow(mat.sig.new)) {
    for (j in 1:ncol(mat.sig.new)) {
      if(rownames(mat.sig.new)[i]==colnames(mat.sig.new)[j]){
        mat.sig.new[i,j] = 0
      }
    }
  }

  mat.sig.new <- mat.sig.new[rowSums(mat.sig.new)!=0,colSums(mat.sig.new)!=0]

  pairs_dact <- data.frame(gene1=character(), gene2=character(), DENDELION_p=character(), stringsAsFactors = FALSE) # significant pairs by DACT
  for (i in 1:nrow(mat.sig.new)) {
    for (j in 1:ncol(mat.sig.new)) {
      if(mat.sig.new[i,j]!=0){
        pairs_dact[nrow(pairs_dact)+1,] <- c(colnames(mat.sig.new)[j],
                                             rownames(mat.sig.new)[i],
                                             mat.p[rownames(mat.sig.new)[i],colnames(mat.sig.new)[j]])
      }
    }
  }

  pairs_dact$wgs_gene2 <- p.wes[match(pairs_dact$gene2, names(p.wes))]



  ####################### combine into one loci ###################################

  # pair = pairs_dact[,c('gene1','gene2')]



  #remove gene1 is NA (there is no corresponding significant cis-gene mapped to the candidate snp)
  cat(paste0('\nCombineing gene loci\n'))
  gene.pair = pairs_dact[which(is.na(pairs_dact$gene1)!=TRUE),]

  gene.pair = gene.pair[!duplicated(gene.pair[,c(1,2)]),]


  gene1.out = strsplit(as.character(gene.pair$gene1), '\\.')

  gene1.info.out = do.call(rbind, gene1.out)
  gene.pair$gene1 = gene1.info.out[,1]


  gene.pair$region = gene.pair$gene1



  gene1_list = gene.pair$gene1

  ref.subset = ref.table.keep[ref.table.keep$gene_name%in%gene1_list,]

  dup.chr = names(table(ref.subset$Chromosome))[which(table(ref.subset$Chromosome)>=2)]

  ref.subset.1 = ref.subset[ref.subset$Chromosome%in%dup.chr,]

  if(nrow(ref.subset.1)!=0){

    for (chr in 1:length(unique(ref.subset.1$Chromosome))) {
      ref.subset.2 = ref.subset.1[ref.subset.1$Chromosome%in%unique(ref.subset.1$Chromosome)[chr],]

      ref.subset.2$region = NA
      ref.subset.2 = ref.subset.2[order(ref.subset.2$start),]

      window = 1

      g1 = 1

      while (g1<= nrow(ref.subset.2)) {

        start = ref.subset.2$start[g1]-500000
        end = ref.subset.2$end[g1]+500000

        loc.pos.1 <- which(ref.subset.2$start>start & ref.subset.2$end<end)

        ref.subset.2$region[loc.pos.1] = window

        g.end = max(loc.pos.1)

        if(g.end == g1){
          g1 = g1+1
          window = window + 1
        }else{
          g1 = g.end
        }


      }

      region_list = unique(ref.subset.2$region)

      for (r1 in 1:length(region_list)) {
        gene.list <- ref.subset.2$gene_name[which(ref.subset.2$region==r1)]

        gene.sub = gene.pair[gene.pair$gene1%in%gene.list,]



        gene.sub = gene.sub%>%distinct()
        gene.sub = gene.sub[order(gene.sub$gene1),]

        # gene.pair$comb_trans_p[(gene.pair$gene1%in%gene.sub$gene1)] = paste(gene.sub$trans_p, collapse = ", ")

        gene.pair$region[(gene.pair$gene1%in%gene.sub$gene1)] = paste(unique(gene.sub$gene1), collapse = ", ")
      }

    }
  }


  gene2 = gene.pair[,"wgs_gene2"]

  imp.var = names(p.wes)[p.wes<=threshold]


  sig_gene2 = gene.pair[which(gene2<threshold),"gene2"]
  non_sig.gene2 = gene.pair[which(gene2>=threshold),"gene2"]

  sig_gene2 = unique(sig_gene2)
  non_sig.gene2 = unique(non_sig.gene2)



  return(list(pairs_dact=pairs_dact,
              gene.pair = gene.pair,
              sig_gene2 = sig_gene2,
              non_sig.gene2 = non_sig.gene2))
}




#' @title Generate figures of DENDELION's result
#' @description Generate analysis plots related to DENDELION results
#' @param gene.pair data frame of identified gene pair, output by calc_pair().
#' @param sig_gene2 gene 2 identified by DENDELION which are also WES significant genes, output by calc_pair().
#' @param non_sig.gene2 gene 2 identified by DENDELION which are NOT WES significant genes, output by calc_pair()
#' @param conservation.score data.frame including the conservation scores which should contain columns: 1. 'pLI.score', 2. 'EDS', 3. 'RVIS', 4. 'GeneSymbol'
#' @param pic_dir the directory to save figures
#' @param p.wes a vector of length K. p value for gene 2 -> trait. Both p.trans and p.wes should contain the gene name.
#' @param eta.wgs the significant threshold for WES gene
#'
#' @return a list includes identified pairs, combined pairs (loci based), identified disease proximal genes (also WES significant), and identified disease proximal genes (NOT WES significant)
#' @export

gen_fig = function(gene.pair, p.wes, eta.wgs=1e-5,  pic_dir){





  imp = which(p.wes<=eta.wgs) #total number of WES significant genes 383

  WES.imp = names(p.wes)[imp] #383


  data = data.frame(gene = c(WES.imp, na.omit(unique(gene.pair[,"gene2"]))),
                    category = c(rep('WES_Genes', length(WES.imp)), rep('DENDELION_Genes', length(na.omit(unique(gene.pair[,"gene2"]))))))



  ############################# Network visualization ####################

  result.pair = gene.pair[,c('region', 'gene2')]


  result.pair = result.pair%>%distinct()


  sig.pair.dact = c()

  for (i in 1:nrow(result.pair)) {

    sig.pair.dact <- c(sig.pair.dact,
                       paste0(result.pair[i,'region'],'-',result.pair[i,'gene2']))

  }


  multiple_gene1 = gene.pair[duplicated(gene.pair[,c(1)]),]

  out.gene1 <- strsplit(as.character(multiple_gene1$gene1),', ')

  info.out.gene1 = as.vector(do.call(rbind, out.gene1))

  node.g1 = unique(c(result.pair[,'region'], result.pair[, 'gene2']))

  colored_gene = node.g1[node.g1%in%WES.imp]


  net.g1 <- graph_from_data_frame(d=result.pair[,c('region','gene2')], vertices=node.g1, directed=T)
  V(net.g1)$type = rep('Not', length(node.g1))
  V(net.g1)$type[which(node.g1%in%WES.imp)] = rep('Important', length(which(node.g1%in%WES.imp)))


  E(net.g1)$type = rep('Not', length(E(net.g1)))



  comps = components(net.g1)

  comp_size = sizes(comps)

  node_sizes = log(degree(net.g1))+1


  pdf(file = paste0(pic_dir, '/',"Network.pdf"),
      width = 16, height = 12)

  plot(net.g1, edge.arrow.size= 0.1, edge.curved=0,
       # vertex.shape="none",
       vertex.size=node_sizes,
       vertex.color="orange", vertex.frame.color="#555555",
       edge.color=c("dark red", "slategrey")[(E(net.g1)$type=="Not")+1],
       # vertex.label = labels,
       vertex.label.color = c( "#bf484a", "#50649d")[1+(V(net.g1)$type=="Not")],
       # layout = layout,
       vertex.label.dist = 0.1,
       vertex.label.cex = 0.3)

  dev.off()

}

########## Useful Functions for DENDELION ##########

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
