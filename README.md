# DANDELION



### **DANDELION: Identification of Candidate Disease Proximal Genes**

**Description**

`DANDELION` is a statistical framework designed to identify candidate *disease-proximal genes* by integrating whole-exome sequencing (WES) and eQTL data. It quantifies the relationship between distal and proximal genes through trans-association p-values, combines them with gene–trait association evidence, and outputs a ranked set of significant gene pairs that may mediate disease risk.

For clarity, disease *distal genes* (putative regulatory genes) are denoted as **gene1**, while *proximal genes* (potential mediators) are denoted as **gene2** within the function.

**Details**

DANDELION takes as input:

* A matrix of **p-values** (`p.trans`) representing distal-to-proximal (gene1 → gene2) associations, where rows correspond to proximal genes and columns to distal genes.
* A named numeric vector of **WES-based p-values** (`p.wes`) for proximal gene–trait associations.
* A vector of **candidate distal genes** (`gene1.list`), corresponding to the subset of distal genes under investigation.

The method further incorporates genomic reference information from:

* `ref.table`, containing gene annotation (gene symbol, type, chromosome, start, end), or
* `SNP.ref`, if the exposure type (`gene1.type`) is specified as `"SNP"`.

DANDELION filters out genes beyond a user-defined cis-distance (default: 5 Mb) and applies a Benjamini–Hochberg FDR correction (default FDR = 0.1) to identify statistically significant trans associations.

The function returns:

1. A list of significant **distal (gene1)** candidates.
2. A data frame of **identified trans gene pairs**.
3. A matrix of **DANDELION-adjusted p-values** summarizing the joint evidence across datasets.

**Output**

The results highlight potential mediator genes that bridge genetic variation and disease phenotypes, providing a robust, interpretable framework for integrating genomic and transcriptomic association signals in complex trait studies.




# Installation


You can install the development version of
`DANDELION` from Github via the `devtools` package. I suppose using
the `remotes` package would work as well.

Before installation of COTA, you are also requested the below packages:
``` r
install.packages(c('qvalue', 'data.table', 'stringr', 'tidyr', 'AnnotationDbi', 'org.Hs.eg.db', 'ggplot2', 'igraph', 'VennDiagram', 'biomaRt', 'plyr', 'dplyr'), dependencies=TRUE)

```

``` r
devtools::install_github("mxxptian/DANDELION")
```

# Example

You can access the example code and data through: https://www.dropbox.com/scl/fo/b1r27fxqlq84ywory8qmo/ABRuFJpz1FIAy9vluMhLAMU?rlkey=n3pzgkdd0fqz9waamidlyfcfy&dl=0
```r


library(qvalue)
library(data.table)
library(stringr)
library(plyr)
library(tidyr)
library(tidyverse)
library(readr)
library(gprofiler2)
library(VennDiagram)
library(igraph)
library(DANDELION)  
library(biomaRt)
library(rsnps)

################################# Gene->Gene->Trait CASE ######################################
# Here we use UKBB summary statistics and GBAT result as an example for gene-gene case

# load UKBB summary statistics
dat = read.table(gzfile('GCST90082964_buildGRCh38.tsv.gz'), header = TRUE)

dat_M3.1 = dat[dat$effect_allele=='M3.1',]

head(dat_M3.1)

out.pheno <- strsplit(as.character(dat_M3.1$Name),'\\.') 
info.out.pheno = do.call(rbind, out.pheno)

out.new <- strsplit(as.character(info.out.pheno[,1]),'\\(')
info.out.new = do.call(rbind, out.new)

dat_M3.1$Name = info.out.new[,1]

dat_M3.1_new = dat_M3.1[!dat_M3.1$Name%in%names(which(table(dat_M3.1$Name)!=1)),]

head(dat_M3.1_new)

cau = data.frame(pval = dat_M3.1_new$p_value)
rownames(cau) = dat_M3.1_new$Name

p.wgs <- cau[,1]
p.wgs <- na.omit(p.wgs)
names(p.wgs) = dat_M3.1_new$Name

# load GBAT data
p.trans <- read.delim('/GBAT/trans_pval/cor_chr1_all.txt', sep = '\t')

for (i in 2:22) {
  temp = read.delim(paste0('/GBAT/trans_pval/cor_chr',i,'_all.txt'), sep = '\t')
  p.trans = cbind(p.trans,temp)
}
rm(temp)

threshold = 0.05/nrow(cau)
eta.wgs = threshold

# load reference table including gene position and gene type. Note!!! this position is GRCH 37
ref.table = read.delim('gene_position.txt', sep = '\t')

## only keep gene type in (lincRNA, protein_coding)
ref.table.keep <- ref.table[ref.table$type %in% c('lincRNA', 'protein_coding'),]

## remove genes on chr M, X, and Y
ref.table.keep <- ref.table.keep[!(ref.table.keep$Chromosome %in% c('chrM', 'chrX', 'chrY')),]

ref.table.keep <-ref.table.keep[!duplicated(ref.table.keep$gene_name),]

## standardize gene name
out.ref <- strsplit(as.character(ref.table.keep$gene_name),'\\.')
out.ref.new = do.call(rbind, out.ref)
ref.table.keep$gene_name = out.ref.new[,1]

# Here, we consider all gene1 included in UKBB. You may specify a vector of candidate genes 1
gene1.list = dput(colnames(p.trans))  

# STEP 1: we applied DANDELION with UKBB and GBAT data with target fdr level at 0.1
result = med_gene(p.trans, p.wgs, ref.table, gene1.list, target.fdr=0.1, dist=5e6, gene1.type = 'Gene')

# STEP 2: we obtained the pair identified by DANDELION
result.pair = calc_pair.gene(result$mat.sig, result$mat.p, p.wgs, result$gene1, ref.table.keep, eta.wgs=threshold)

# STEP 3: we plotted the pair identified by DANDELION
setwd('/Users/px/Desktop/test/example/')

## load conservation score table
conv.table = read_excel('conservation_scores.xlsx')

## set saving directory of the figures
pic_dir = '/Users/px/Desktop/test/example/'

gen_fig(result.pair$gene.pair, p.wgs, eta.wgs=threshold, pic_dir) 

################################# SNP->Gene->Trait CASE ######################################

# Here we use UKBB summary statistics and GBAT result as an example for snp-gene case

# You can access this file by (https://www.dropbox.com/scl/fo/b1r27fxqlq84ywory8qmo/ABRuFJpz1FIAy9vluMhLAMU?rlkey=n3pzgkdd0fqz9waamidlyfcfy&dl=0).
# Note!!! The file size is very large.
trans.p = fread('2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt', sep = '\t')

# This matrix includes the p-values from snp -> core genes
load('trans-eQTL_genetype_genename.Rdata')
p.trans = t(pval.sub)

# This file includes the significant cis-gene
load('snp2sig_cisgene.Rdata')

load('GCST90082964.Rdata')
head(dat_M3.1)

out.pheno <- strsplit(as.character(dat_M3.1$Name),'\\.')
info.out.pheno = do.call(rbind, out.pheno)

out.new <- strsplit(as.character(info.out.pheno[,1]),'\\(')
info.out.new = do.call(rbind, out.new)

dat_M3.1$Name = info.out.new[,1]

dat_M3.1_new = dat_M3.1[!dat_M3.1$Name%in%names(which(table(dat_M3.1$Name)!=1)),]

head(dat_M3.1_new)

cau = data.frame(pval = dat_M3.1_new$p_value)
rownames(cau) = dat_M3.1_new$Name

p.wgs <- cau[,1]
p.wgs <- na.omit(p.wgs)
names(p.wgs) = dat_M3.1_new$Name

# Be care about the build version of the position
ref.table = read.delim('gene_position.txt', sep = '\t')

threshold = 0.05/nrow(cau)
eta.wgs = threshold

## only keep gene type in (lincRNA, protein_coding)
ref.table.keep <- ref.table[ref.table$type %in% c('lincRNA', 'protein_coding'),]

## remove genes on chr M, X, and Y
ref.table.keep <- ref.table.keep[!(ref.table.keep$Chromosome %in% c('chrM', 'chrX', 'chrY')),]

ref.table.keep <-ref.table.keep[!duplicated(ref.table.keep$gene_name),]

## standardize gene name
out.ref <- strsplit(as.character(ref.table.keep$gene_name),'\\.')
out.ref.new = do.call(rbind, out.ref)
ref.table.keep$gene_name = out.ref.new[,1]

# Here, we consider all gene1 included in UKBB. You may specify a vector of candidate genes 1
gene1.list = dput(colnames(p.trans)) 

# STEP 1: we applied DANDELION with UKBB and eQTLGen data with target fdr level at 0.1
result = med_gene(p.trans, p.wgs, ref.table, gene1.list, target.fdr=0.1, dist=5e6, gene1.type = 'SNP', SNP.ref = trans.p)

# STEP 2: we obtained the pair identified by DANDELION
result.pair = calc_pair.snp(mat.sig = result$mat.sig, mat.p = result$mat.p, p.wgs, gene1 = result$gene1, 
                           uniq_snp = uniq_snp, ref.table.keep = ref.table.keep, eta.wgs=threshold, GRCh = '37')

# STEP 3: we plotted the pair identified by DANDELION
setwd('/Users/px/Desktop/test/example/')

## load conservation score table
conv.table = read_excel('conservation_scores.xlsx')

## set saving directory of the figures
pic_dir = '/Users/px/Desktop/test/example/'

gen_fig(result.pair$gene.pair, p.wgs, eta.wgs=threshold, pic_dir)  



```
