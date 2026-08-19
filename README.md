
`DANDELION` is an R package for identifying **disease-proximal genes (DPGs)** that may mediate the effects of disease-associated loci on disease risk. It provides a causal mediation–inspired gene prioritization framework that integrates trans-regulatory association signals with gene-level disease association evidence.

This package accompanies the study:

> *Trans-regulatory gene mapping prioritizes disease drivers in asthma*.

The release version of DANDELION is archived on Zenodo at DOI: 10.5281/zenodo.19911607.

In the R functions and examples below:

* **gene1** denotes the disease-distal exposure, which can be either a regulatory gene or a SNP.
* **gene2** denotes the disease-proximal candidate mediator and trans-regulatory target.

## Overview of DANDELION R package

`DANDELION` combines two sources of evidence:

1. **Trans-regulatory association p-values**, representing the association between an upstream exposure and a downstream gene.
2. **Gene-level disease association p-values**, such as burden-test p-values from whole-exome sequencing (WES) or other gene-level association tests.

To assess mediation, `DANDELION` implements and adapts the **Divide-Aggregate Composite-null Test (DACT)** framework to the trans-gene regulation setting. Specifically, `DANDELION` combines trans-association p-values with gene-level burden-test p-values to prioritize candidate disease-proximal genes and distal-proximal gene pairs.

Details of the causal mediation assumptions, composite-null structure, p-value computation, and statistical theory are provided in the original DACT paper (Liu et al. 2022).


### 1. Gene-based DANDELION

In gene-based DANDELION, the exposure side contains disease-distal genes. The input trans-association matrix represents gene-to-gene trans-regulatory p-values.

A typical use case is:

```text
disease-distal gene -> trans target gene -> disease trait
```

When raw genotype and gene expression data are available, gene-based trans-regulatory effects can be estimated using approaches such as GBAT. In this setting, predicted genetically regulated expression of a distal gene can be tested against genome-wide target genes to obtain gene-to-gene trans-association p-values.

### 2. SNP-based DANDELION

In SNP-based DANDELION, the exposure side contains SNPs. The input trans-association matrix represents SNP-to-gene trans-regulatory p-values, such as trans-eQTL summary statistics.

A typical use case is:

```text
SNP -> trans target gene -> disease trait
```

This setting is useful when raw genotype and expression data are unavailable but trans-eQTL summary statistics are available.

## Features

* Supports both **gene-based** and **SNP-based** DANDELION analyses.
* Integrates trans-regulatory p-values with gene-level disease association p-values.
* Removes local cis genes using a user-defined genomic window before trans-mediation testing.
* Applies a DACT-style composite-null testing framework to prioritize candidate mediators.
* Returns structured R objects that users can inspect, print, and extract from.
* Provides post-processing functions for organizing significant gene pairs.
* Supports locus-level grouping of nearby exposure genes.
* Provides a simple network visualization function for DANDELION-prioritized relationships.

## Installation

You can install the development version from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("mxxptian/DANDELION")
```

After installation, load the package with:

```r
library(DANDELION)
```

### Optional dependency

The Bioconductor package `qvalue` is used for q-value estimation when available. If `qvalue` is not installed, `DANDELION` falls back to Benjamini-Hochberg adjusted p-values.

```r
# Optional
# install.packages("BiocManager")
BiocManager::install("qvalue")
```

## Main functions

| Function           | Description                                                                                                           |
| ------------------ | --------------------------------------------------------------------------------------------------------------------- |
| `med_gene()`       | Applies DANDELION to identify significant exposure-gene2 mediation pairs. The exposure can be a distal gene or a SNP. |
| `calc_pair.gene()` | Organizes significant gene-based DANDELION results and merges nearby exposure genes into locus-level regions.         |
| `calc_pair.snp()`  | Organizes significant SNP-based DANDELION results and maps SNP exposures to cis genes.                                |
| `gen_fig()`        | Generates a network visualization of DANDELION-prioritized gene-pair results.                                         |

## Required input data

### 1. `p.trans`: trans-association p-value matrix

`p.trans` is a numeric matrix of trans-regulatory association p-values.

* Rows correspond to candidate disease-proximal genes (**gene2**).
* Columns correspond to exposures (**gene1**), either genes or SNPs.
* Row names and column names must be provided.

For gene-based DANDELION:

```text
rows    = candidate proximal genes / trans target genes
columns = distal regulatory genes
```

For SNP-based DANDELION:

```text
rows    = candidate proximal genes / trans target genes
columns = SNPs
```

### 2. `p.wes`: gene-level disease association p-values

`p.wes` is a named numeric vector of gene-level disease association p-values.

Examples include:

* WES burden-test p-values,
* gene-level GWAS p-values,
* other gene-level association statistics.

The names of `p.wes` should match the row names of `p.trans`.

### 3. `ref.table`: gene annotation table

`ref.table` provides gene position and gene type information.

It must contain the following columns:

| Column       | Description                                      |
| ------------ | ------------------------------------------------ |
| `gene_name`  | Gene symbol or gene identifier                   |
| `type`       | Gene type, such as `protein_coding` or `lincRNA` |
| `Chromosome` | Chromosome in `chr1`, `chr2`, ... format         |
| `start`      | Gene start position                              |
| `end`        | Gene end position                                |

`DANDELION` internally keeps autosomal `protein_coding` and `lincRNA` genes and removes genes on `chrM`, `chrX`, and `chrY` for the main trans-mediation analysis.

### 4. `gene1.list`: candidate exposure list

`gene1.list` is a character vector of candidate exposures to analyze.

* For gene-based DANDELION, it should contain candidate distal genes.
* For SNP-based DANDELION, it should contain candidate SNPs.
* Values must be present in `colnames(p.trans)`.

### 5. `SNP.ref`: SNP annotation table for SNP-based DANDELION

When `gene1.type = "SNP"`, `SNP.ref` must be provided.

It must contain:

| Column   | Description                                 |
| -------- | ------------------------------------------- |
| `SNP`    | SNP rsID                                    |
| `SNPPos` | SNP genomic position                        |
| `SNPChr` | SNP chromosome, either `1` or `chr1` format |

### 6. `uniq_snp`: SNP-to-cis-gene mapping table

For `calc_pair.snp()`, `uniq_snp` is used to map SNP exposures to their corresponding cis genes.

It must contain:

| Column       | Description                   |
| ------------ | ----------------------------- |
| `SNP`        | SNP rsID                      |
| `GeneSymbol` | Corresponding cis gene symbol |

## Quick start: gene-based DANDELION

The following toy example demonstrates the gene-based workflow.

```r
library(DANDELION)

set.seed(1)

# Example trans-association p-value matrix.
# Rows are candidate disease-proximal genes; columns are exposure genes.
p.trans <- matrix(runif(60, 0.001, 0.9), nrow = 12, ncol = 5)
rownames(p.trans) <- paste0("G", 1:12)
colnames(p.trans) <- paste0("E", 1:5)

# Example gene-level disease association p-values.
p.wes <- runif(12, 0.001, 0.9)
names(p.wes) <- paste0("G", 1:12)

# Example gene annotation table.
ref.table <- data.frame(
  gene_name = c(paste0("G", 1:12), paste0("E", 1:5)),
  type = "protein_coding",
  Chromosome = "chr1",
  start = seq_len(17) * 1e7,
  end = seq_len(17) * 1e7 + 1000
)

# Run gene-based DANDELION.
res <- med_gene(
  p.trans = p.trans,
  p.wes = p.wes,
  ref.table = ref.table,
  gene1.list = colnames(p.trans),
  gene1.type = "Gene",
  target.fdr = 0.1,
  dist = 5e6
)

res
```

The returned object is a structured `dandelion_result` object. You can extract information directly:

```r
res$gene1
res$mat.sig[1:5, 1:3]
res$mat.p[1:5, 1:3]
```

Organize significant gene-gene pairs:

```r
ref.table.keep <- ref.table[
  ref.table$type %in% c("lincRNA", "protein_coding") &
    !(ref.table$Chromosome %in% c("chrM", "chrX", "chrY")),
]

pair.res <- calc_pair.gene(
  mat.sig = res$mat.sig,
  mat.p = res$mat.p,
  p.wes = p.wes,
  gene1 = res$gene1,
  ref.table.keep = ref.table.keep,
  eta.wgs = 1e-5
)

pair.res
```

## Quick start: SNP-based DANDELION

The following toy example demonstrates the SNP-based workflow.

```r
library(DANDELION)

set.seed(1)

# Example SNP-to-gene trans-association p-value matrix.
# Rows are candidate disease-proximal genes; columns are SNP exposures.
p.trans <- matrix(runif(60, 0.001, 0.9), nrow = 12, ncol = 5)
rownames(p.trans) <- paste0("G", 1:12)
colnames(p.trans) <- paste0("rs", 1:5)

p.wes <- runif(12, 0.001, 0.9)
names(p.wes) <- paste0("G", 1:12)

ref.table <- data.frame(
  gene_name = paste0("G", 1:12),
  type = "protein_coding",
  Chromosome = "chr1",
  start = seq_len(12) * 1e7,
  end = seq_len(12) * 1e7 + 1000
)

SNP.ref <- data.frame(
  SNP = paste0("rs", 1:5),
  SNPPos = seq_len(5) * 2e7,
  SNPChr = "1"
)

uniq_snp <- data.frame(
  SNP = paste0("rs", 1:5),
  GeneSymbol = paste0("G", 1:5)
)

# Run SNP-based DANDELION.
res <- med_gene(
  p.trans = p.trans,
  p.wes = p.wes,
  ref.table = ref.table,
  gene1.list = colnames(p.trans),
  gene1.type = "SNP",
  SNP.ref = SNP.ref,
  target.fdr = 0.1,
  dist = 5e6
)

res
```

Organize significant SNP-gene pairs:

```r
ref.table.keep <- ref.table[
  ref.table$type %in% c("lincRNA", "protein_coding") &
    !(ref.table$Chromosome %in% c("chrM", "chrX", "chrY")),
]

pair.res <- calc_pair.snp(
  mat.sig = res$mat.sig,
  mat.p = res$mat.p,
  p.wes = p.wes,
  gene1 = res$gene1,
  uniq_snp = uniq_snp,
  ref.table.keep = ref.table.keep,
  SNP.ref = SNP.ref,
  eta.wgs = 1e-5
)

pair.res
```

## Output objects

### Output from `med_gene()`

`med_gene()` returns an object of class `dandelion_result`.

| Output       | Description                                        |
| ------------ | -------------------------------------------------- |
| `gene1`      | Exposures with at least one valid DANDELION result |
| `mat.sig`    | Matrix encoding significant gene1-gene2 pairs      |
| `mat.p`      | Matrix of DANDELION p-values                       |
| `gene1.type` | Exposure type, either `Gene` or `SNP`              |
| `target.fdr` | FDR cutoff used for identifying significant pairs  |

### Output from `calc_pair.gene()` and `calc_pair.snp()`

These functions return objects of class `dandelion_pairs`.

| Output          | Description                                               |
| --------------- | --------------------------------------------------------- |
| `pairs_dact`    | Significant DANDELION pairs before locus-level merging    |
| `gene.pair`     | Pair table after organizing exposure genes into loci      |
| `sig_gene2`     | DPGs that are also significant in gene-level burden tests |
| `non_sig.gene2` | DPGs that are not significant in gene-level burden tests  |

## Network visualization

`gen_fig()` generates a network visualization from organized DANDELION gene-pair results.

```r
gen_fig(
  gene.pair = pair.res$gene.pair,
  p.wes = p.wes,
  eta.wgs = 1e-5,
  pic_dir = "path/to/output/folder"
)
```

The function saves `Network.pdf` to `pic_dir` and invisibly returns the output file path.

## Verbose output

By default, DANDELION functions do not print progress messages to the console. To show progress messages, set `verbose = TRUE`:

```r
res <- med_gene(
  p.trans = p.trans,
  p.wes = p.wes,
  ref.table = ref.table,
  gene1.list = colnames(p.trans),
  gene1.type = "Gene",
  verbose = TRUE
)
```

## Analysis directory

The `Analysis/` directory contains public templates and example scripts for running DANDELION simulation and real-data workflows.

During peer review, this repository provides lightweight, GitHub-friendly templates that demonstrate the expected input structure, DANDELION function calls, output format, and evaluation strategy. Full manuscript-scale data-generation scripts, internal HPC pipelines, real-data analysis scripts, and intermediate results are maintained privately during peer review to protect code provenance and avoid unauthorized redistribution.

### `Analysis/simulation/`

This folder contains simulation scripts used to evaluate the statistical performance of DANDELION under controlled mediation settings.

Available scripts include:

- `simulation_compare.R`  
  Implements the manuscript-scale simulation framework comparing DANDELION with ARCHIE. The simulation generates three independent samples to mimic realistic scenarios where genetic variants, gene expression, and complex traits are measured in non-overlapping cohorts. The framework evaluates gene-based DANDELION, SNP-based DANDELION, and ARCHIE in terms of statistical power and global false discovery rate (FDR) control.

- `sim_global_qvalue.R`  
  Evaluates the global FDR calibration of DANDELION under controlled simulation settings. This simulation generates disease-distal gene to target gene mediation structures and assesses gene-level and pair-level discovery performance across prespecified FDR thresholds using q-value based inference.

- `global_fdr_simulation_template.R`  
  Provides a simplified example demonstrating how to evaluate DANDELION at both the gene level and exposure-gene pair level, including performance metrics such as true positive proportion (TPP), FDR, and the number of selected genes/pairs.

- `simulation_compare_template.R`  
  Provides a lightweight example illustrating how to organize simulation workflows comparing different DANDELION analysis strategies using shared input/output structures.

The simulation framework evaluates:

- **Statistical power**: the proportion of true disease-proximal genes correctly identified.
- **False discovery rate control**: the ability of DANDELION to maintain the prespecified FDR levels.
- **Gene-level and pair-level performance**: evaluation of disease-proximal gene discovery and exposure-gene pair prioritization.
- **Method comparison**: comparison of DANDELION with alternative approaches under matched simulation settings.

The manuscript-scale simulation implements the DANDELION testing framework and downstream inference procedures used in the study, while the template scripts provide simplified examples for applying DANDELION in user-defined simulation settings.

---

### `Analysis/real_data/`

This folder contains scripts for applying DANDELION to real disease-trait datasets using SNP-based trans-regulatory summary statistics and gene-level disease association results.

Available scripts include:

- `Analysis/real_data/run_dandelion_snp_pipeline.R`  
  A general SNP-based DANDELION workflow illustrating the application of DANDELION using trans-eQTL summary statistics, SNP annotation, SNP-to-gene mapping information, and trait-specific gene-level disease association results.

- `Analysis/real_data/run_dandelion_asthma_pipeline.R`  
  The trait-specific SNP-based DANDELION workflow used for the asthma analysis reported in the manuscript.

- `Analysis/real_data/run_dandelion_real_data.R`  
  Implements the generalized real-data analysis workflow used for multi-trait analyses in the manuscript. The same pipeline was applied across multiple disease traits by changing trait-specific input files.

The real-data workflows include:

- loading trans-regulatory eQTL summary statistics;
- preparing SNP annotation and SNP-to-gene mapping information;
- preparing trait-specific gene-level disease association p-values;
- applying `med_gene()` for trans-regulatory mediation analysis;
- performing downstream SNP-gene pair prioritization;
- saving DANDELION results, significant SNP-gene pairs, and summary outputs.


## Citation

If you use `DANDELION`, please cite the following work.

The statistical testing method in `DANDELION` builds on the Divide-Aggregate Composite-null Test (DACT):

> Liu, Z., Shen, J., Barfield, R., Schwartz, J., Baccarelli, A. A., & Lin, X. (2022). Large-Scale Hypothesis Testing for Causal Mediation Effects with Applications in Genome-wide Epigenetic Studies. *Journal of the American Statistical Association*, 117(537), 67–81. [https://doi.org/10.1080/01621459.2021.1914634](https://doi.org/10.1080/01621459.2021.1914634)

The trans-gene regulation paper is here:

> Salamone, I.M.†, Tian, P.†, Qi, Z., Zhao, J., Zhang, L., Tan, Q., Li, J., Michael, A.N., Thornburg, A.G., Sakabe, N.J. and Minogue, M., 2026. Trans-regulatory gene mapping prioritizes disease drivers in asthma. Cell.

> Note: † indicates authors contributed equally to this work.

## Authorship and copyright

**R implementation:** Peixin Tian

**Copyright:** Copyright (C) 2025 Peixin Tian.

**License:** GPL-3.0. See the `LICENSE` file.

## License

This package is distributed under the GNU General Public License v3.0. See the `LICENSE` file for details.
