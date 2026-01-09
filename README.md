# DANDELION



### **DANDELION: Identification of Candidate Disease Proximal Genes**
`DANDELION` is the software implementation accompanying the study “Leveraging trans-gene regulation to prioritize disease-proximal genes” and provides an end-to-end implementation of the proposed mediation-based gene prioritization framework.

## Description

`DANDELION` identifies disease-proximal genes (DPGs) that mediate the effects of disease-associated loci on disease risk. It models DPGs as mediators and trans regulatory targets of disease associated loci, which is named as disease distal genes (Figure below). Using a causal mediation framework, it integrate gene effects on disease obtained from burden tests of whole-exome sequencing (WES) and trans regulatory signals in disease-relevant tissues and cell types.
In the R function and examples below, disease distal genes (trans regulators) are denoted as gene1, while disease proximal genes ( mediators and trans regulation targets) are denoted as gene2.

For clarity, disease *distal genes* (putative regulatory genes) are denoted as **gene1**, while *proximal genes* (potential mediators) are denoted as **gene2** within the function.

**Details**

`DANDELION` takes as input:

* A matrix of **trans regulation p-values** (`p.trans`) representing distal-to-proximal (gene1 → gene2) associations, where rows correspond to proximal genes and columns to distal genes. The trans signals should be from tissues that are relevant to the disease of interest. We provide the trans p value matrices in whole blood datasets(eQTLGEN and DGN) in this package. Users interested in diseases relevant to whole blood can use them directly to nominate candidate DPGs.
* A named numeric vector of **WES-based p-values** (`p.wes`) for proximal gene–trait associations.
* A vector of **candidate distal genes** (`gene1.list`), corresponding to the set of distal genes, typically genes showing association signals with the disease.


`DANDELION` filters out genes beyond a user-defined cis-distance (default: 5 Mb) and applies a Benjamini–Hochberg FDR correction (default FDR = 0.1) to identify statistically significant trans associations. 
`DANDELION` uses a user-defined trans-regulation distance (default: >5 Mb). For a disease and a disease-associated SNP, all genome-wide genes in trans can serve as potential mediators.  `DANDELION` uses a causal mediation method DACT to decompose all paths of no mediation into three NULL cases and robustly identifies mediation paths with high statistical power. It uses Benjamini–Hochberg FDR correction (default FDR = 0.1) to identify statistically significant mediation paths.

The function returns:

1. A list of significant **distal (gene1)** candidates.
2. A data frame of **identified trans gene pairs**.
3. A matrix of **DANDELION-adjusted p-values** summarizing the joint evidence across datasets.

**Output**

The results highlight potential mediator genes that bridge genetic variation and disease phenotypes, providing a robust, interpretable framework for integrating genomic and transcriptomic association signals in complex trait studies.
The mediators on the significant mediation path are the DPGs. `DANDELION` can generate trans regulatory networks of DPGs based on the significant mediation paths (gene1->gene2).

## Analysis

The `Analysis/` directory contains all analyses performed in this project.  
It is organized into **simulation** and **real data**, each serving a distinct methodological purpose.

---

###  Analysis/simulation/

This folder contains simulation studies designed to **evaluate the statistical performance of the DANDELION framework** under controlled data-generating mechanisms.

The simulations focus on assessing:

- **Statistical power (simulation_power.R)**  
  The ability of DANDELION to correctly identify disease-proximal genes (gene2) that truly mediate the effects of gene1 on the trait.

- **False discovery rate (FDR) control (simulation_fdr.R)**  
  The accuracy of DANDELION in maintaining target FDR levels across a wide range of signal sparsity and effect size configurations.


- **Compare with competing method (Archie) (simulation_compare.R)**  
  Performance is evaluated at both:
  - the **global level** (aggregating signals across all gene1), and
  - the **per-gene / per-SNP level**, ensuring stable behavior across analytical resolutions.

All simulation pipelines use the **same DANDELION test statistic** as implemented in the real data analysis (`med_gene()`), ensuring full methodological consistency.

---

### Analysis/real_data 

This folder contains all analyses based on **real human genetic data**, serving as the primary application of the DANDELION framework.

In particular, it includes:

- **Blood trait analyses (blood trait analysis.R)**  
  Application of DANDELION to blood-related traits by integrating:
  - trans-regulatory association statistics,
  - whole-exome sequencing (WES)–based gene-level association results,
  - identification of candidate disease-proximal genes that mediate the effects of disease-associated loci on blood traits.
All real data analyses follow the same statistical definitions and analytical pipeline validated in the simulation studies.

- **Enrichment_analysis (enrichment_func.R)**

  Enrichment analyses applied to genes prioritized by **DANDELION** in real data.

  For each blood-related trait, we assess whether **DANDELION-identified disease-proximal genes** are enriched for **known Mendelian disease genes**, compared      with randomly sampled gene sets matched on gene length.  
  
  Enrichment is evaluated for multiple gene categories, including:
  - all disease-proximal genes identified by DANDELION,
  - disease-proximal genes that are significant in whole-exome sequencing (WES),
  - disease-proximal genes that are not WES-significant,
  - burden-only genes identified by WES but not prioritized by DANDELION.

This analysis provides biological validation of DANDELION by testing whether its prioritized genes are preferentially associated with known Mendelian disease mechanisms beyond standard burden-based approaches.

---

###  Methodological Consistency

Simulation studies and real data analyses are tightly coupled:

> Simulation studies establish the statistical validity and operating characteristics of DANDELION, while real data analyses demonstrate its practical utility.  
> All analyses are performed using the same core DANDELION statistic and modeling assumptions.



# Installation


You can install the development version of
`DANDELION` from Github via the `devtools` package. I suppose using
the `remotes` package would work as well.

Before installation of DANDELION, you are also requested the below packages:
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
setwd('/PATH_TO_FIG_FOLDER')

## load conservation score table
conv.table = read_excel('conservation_scores.xlsx')

## set saving directory of the figures
pic_dir = '/PATH_TO_FIG_FOLDER'


gen_fig(result.pair$gene.pair, p.wgs, eta.wgs=threshold, pic_dir)  



```
