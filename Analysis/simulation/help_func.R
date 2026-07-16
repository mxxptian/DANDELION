#library(qvalue)
library(data.table)
library(stringr)
library(tidyr)




med_gene <- function(p.trans, p.wgs, ref.table, target.fdr=0.1, dist=5e6){
  
  # p.trans: p value for trans genes. M*N dim matrix. 
  # M rows for gene 2. N columns for gene 1.
  # Each entry is the p-value of gene 1 -> gene 2.
  
  # p.wgs: vector of length K. p value for gene 2 -> trait.
  # Both p.trans and p.wgs should contain the gene name.
  
  # ref.table: reference data for position information.
  # should contain five columns:
  # 1. gene_name: gene name
  # 2. type: gene type
  # 3. Chromosome
  # 4. start: start position of gene
  # 5. end: end position of gene
  
  ##### Step 1: match the common gene. #####
  
  cat('\n-----Start matching gene 2-----\n')
  
  ## match gene 2 of p.trans and p.wgs
  common_gene2 <- intersect(rownames(p.trans), names(p.wgs))
  
  ## keep only lincRNA and protein coding genes
  ref.table.keep <- ref.table[ref.table$type %in% c('lincRNA', 'protein_coding'),]
  
  ## remove genes on chr M, X, and Y
  ref.table.keep <- ref.table.keep[!(ref.table.keep$Chromosome %in% c('chrM', 'chrX', 'chrY')),]
  
  
  ref.table.keep = ref.table.keep[!duplicated(ref.table.keep$gene_name),]
  
  
  ## match common_gene2 to reference table
  common_gene2 <- intersect(common_gene2, ref.table.keep$gene_name)
  
  ## keep common gene2 for p.trans
  p.trans.new <- p.trans[rownames(p.trans)%in%common_gene2,]
  
  ## keep common gene2 for p.wgs
  p.wgs.new <- p.wgs[names(p.wgs)%in%common_gene2]
  
  ## reorder p.trans.new and p.wgs.new
  p.trans.new <- p.trans.new[match(names(p.wgs.new), rownames(p.trans.new)),]
  
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
  
  cat(paste0('\n-----Start DACT for ',M,' genes-----\n'))
  
  for (i in 1:ncol(p.trans.new)) {
    
    gene1 = colnames(p.trans.new)[i]  # gene 1 name
    
    if(gene1 %in% ref.table.keep$gene_name){ # if gene 1 in lincRNA or protein coding
      
      pos.info <- ref.table.keep[ref.table.keep$gene_name==gene1,]  # position information for gene 1
      start.pos <- pos.info$start
      end.pos <- pos.info$end
      chr <- pos.info$Chromosome
      
      start.critical <- start.pos - dist  # start position of critical region
      end.critical <- end.pos + dist  # end position of critical region
      
      target.genes <- ref.table[ref.table$Chromosome==chr,]  # genes on the same chromosome
      
      loc.pos.1 <- which(target.genes$start>start.critical & target.genes$start<end.critical)  # genes whose start in critical region
      loc.pos.2 <- which(target.genes$end>start.critical & target.genes$end<end.critical)  # genes whose end in critical region
      loc.pos <- union(loc.pos.1, loc.pos.2)  # cis gene location
      
      gene.cis <- target.genes$gene_name[loc.pos] # cis gene name
      
      gene.trans <- setdiff(common_gene2, gene.cis) # trans gene name
      
      
      ##### Step 3: DACT step #####
      
      p_a = p.trans.new[gene.trans, i]; names(p_a) <- gene.trans
      p_b = p.wgs.new[gene.trans]; p_b[which(p_b==1)] <- 0.99 
      
      Z_a = stats::qnorm(p_a,lower.tail = F)
      Z_b = stats::qnorm(p_b,lower.tail = F)
      
      pi0a = 1-nonnullPropEst(Z_a, 0, 1)
      pi0b = 1-nonnullPropEst(Z_b, 0, 1)
      
      if(is.na(pi0a)==FALSE & is.na(pi0b)==FALSE & pi0a>=0 & pi0b>=0){
        
        if(pi0a > 1){pi0a = 1}
        if(pi0b > 1){pi0b = 1}
        
        p.mat = cbind(p_a,p_b)
        p3 = (apply(p.mat,1,max))^2
        wg1 = pi0a*(1-pi0b)
        wg2 = (1-pi0a)*pi0b
        wg3 = pi0a*pi0b
        wg.sum = wg1 + wg2 + wg3
        wg.std = c(wg1,wg2,wg3)/wg.sum
        p_dact = wg.std[1]*p_a + wg.std[2]*p_b + wg.std[3]*p3
        
        s = try(qvalue(p_dact))
        
        if('try-error'%in% class(s)){
          s = qvalue(p_dact, pi0=1)
        }
        
        q_dact = s$qvalues
        
        names(q_dact) = gene.trans 
        
        vec = rep(0, length(rownames(p.trans.new)))  # vector of 0&1. 1 stands for significant gene 2 by DACT.
        
        sig_dact = names(q_dact)[which(q_dact<=target.fdr)]
        
        vec[rownames(p.trans.new)%in%sig_dact]=1
        
        mat.sig[,i] <- vec
        
        vec.p = rep(NA, length(rownames(p.trans.new)))
        
        # vec.p[rownames(p.trans.new)%in%gene.trans] = p_dact
        vec.p[match(gene.trans, rownames(p.trans.new))] = p_dact
        
        mat.p[,i] <- vec.p
        
        gene1.candidate = c(gene1.candidate, gene1)
      }
    }
    
    if((i%%floor(M/100))==0){
      percent = i%/%floor(M/100)
      cat(paste0('\nComplete ',percent,'% genes!\n'))
      cat(paste0('Number of Gene 1 detected: ', length(gene1.candidate), '\n'))
    }
    
  }
  
  return(list(gene1=gene1.candidate, mat.sig=mat.sig, mat.p=mat.p))
}


# results <- med_gene(p.trans, p.wgs, ref.table = ref.table, target.fdr = 0.2)

calc_pair <- function(mat.sig, gene1, p.trans, p.wgs, eta.wgs=1e-5, eta.trans=2.664087e-07){
  
  # mat.sig: matrix encoding significant pair, output by med_gene()
  # gene1: gene 1 names used for analysis, output by med_gene()
  # p.trans: p value for trans genes. M*N dim matrix. 
  # eta.wgs: threshold for wgs significant genes.
  # eta.trans: threshold for trans significant genes.
  
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
  
  pairs_dact <- data.frame(gene1=character(), gene2=character(), trans_p=character()) # significant pairs detected by COTA
  for (i in 1:nrow(mat.sig.new)) {
    for (j in 1:ncol(mat.sig.new)) {
      if(mat.sig.new[i,j]!=0){
        pairs_dact[nrow(pairs_dact)+1,] <- c(colnames(mat.sig.new)[j], 
                                             rownames(mat.sig.new)[i],
                                             p.trans[rownames(mat.sig.new)[i],colnames(mat.sig.new)[j]])
      }
    }
  }
  
  pairs_dact$trans_p <- as.numeric(pairs_dact$trans_p)
  pairs_dact$wgs_gene1 <- p.wgs[match(pairs_dact$gene1, names(p.wgs))]
  pairs_dact$wgs_gene2 <- p.wgs[match(pairs_dact$gene2, names(p.wgs))]
  
  ##### Consider gene 2 #####
  
  # (a)
  num_pair_sig_gene2_wgs <- sum(na.omit(pairs_dact$wgs_gene2<=eta.wgs))
  
  # (b)
  num_pair_sig_gene2_gene1_wgs <- sum(na.omit(pairs_dact$wgs_gene2<=eta.wgs & pairs_dact$wgs_gene1<=eta.wgs))
  
  # (c)
  num_pair_sig_gene2_wgs_trans <-sum(na.omit(pairs_dact$wgs_gene2<=eta.wgs & pairs_dact$trans_p<=eta.trans))
  
  
  ##### Consider gene 1 #####
  
  # (a)
  num_pair_sig_gene1_wgs <- sum(na.omit(pairs_dact$wgs_gene1<=eta.wgs))
  
  # (b)
  num_pair_sig_gene1_gene2_wgs <- sum(na.omit(pairs_dact$wgs_gene2<=eta.wgs & pairs_dact$wgs_gene1<=eta.wgs))
  
  # (c)
  num_pair_sig_gene1_wgs_trans <- sum(na.omit(pairs_dact$wgs_gene1<=eta.wgs & pairs_dact$trans_p<=eta.trans))
  
  
  return(list(pairs_dact=pairs_dact, 
              num_pair_sig_gene2_wgs=num_pair_sig_gene2_wgs,
              num_pair_sig_gene2_gene1_wgs=num_pair_sig_gene2_gene1_wgs,
              num_pair_sig_gene2_wgs_trans=num_pair_sig_gene2_wgs_trans,
              num_pair_sig_gene1_wgs=num_pair_sig_gene1_wgs,
              num_pair_sig_gene1_gene2_wgs=num_pair_sig_gene1_gene2_wgs,
              num_pair_sig_gene1_wgs_trans=num_pair_sig_gene1_wgs_trans))
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
