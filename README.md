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
