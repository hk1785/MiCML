# MiCML

**Title:** Microbiome Causal Machine Learning (MiCML)

**Version:** 1.0.1

**Description:** MiCML (microbiome causal machine learning) is a cloud computing platform that streamlines related data processing and analytic procedures for the analysis of treatment effects using microbiome profiles on user-friendly web environments. MiCML is in particular unique with the up-to-date features of (i) batch effect correction to mitigate systematic variation in collective large-scale microbiome data due to the differences in their underlying batches (e.g., lab or study environments), and (ii) causal machine learning to estimate treatment effects with consistency and then discern microbial taxa that enhance (or lower) the efficacy of the primary treatment. We also stress that MiCML can handle the data from either randomized controlled trials or observational studies. MiCML can be a useful analytic tool for microbiome-based personalized medicine. MiCML consists of three Data Processing modules, (i) Data Input, (ii) Batch Effect Correction & Quality Control, and (iii) Data Transformation and three Data Analysis modules, (i) Descriptive Analysis, (ii) Generalized Linear Models, and (iii) Causal Machine Learning.

**NeedsCompilation:** No

**Depends:** R(≥ 4.4.0)

**Imports:** 'ape', 'phyloseq', 'shiny', 'shinydashboard', 'shinyjs', 'shinyWidgets', 'session'

**License:** General Public License 3 (GPL3) 

**URLs:** Web Server (http://micml.micloud.kr), GitHub (http://github.com/hk1785/micml) 

**Maintainer:** Hyunwook Koh <hyunwook.koh@stonybrook.edu>

**Reference:** Koh H, Kim J, Jang H. MiCML: A causal machine learning cloud platform for the analysis of treatment effects using microbiome profiles (*In review*). 

<hr>

## GitHub Repository Contents

* **Data** - In this directory, example microbiome data are stored in a unified format, called phyloseq (see 'biom.Rdata' in 'Data/Phyloseq') as well as individual files (see feature table (otu.tab.txt), taxonomic table (tax.tab.txt), metadata/sample information (sam.dat.txt), and phylogenetic tree (tree.tre) in 'Data/Individual').
  
* **Source** - In this directory, all the R functions that are needed to run MiCML are stored.

* **www** - In this directory, some photos that are used to decorate the GUI of MiCML are stored.

* **app.R** - In this file, all the central codes to control for user-interfaces and server functions of MiCML are stored.

<hr>

## Required Data Components

MiCML requires four data components: feature table, taxonomic table, metadata/sample information, and phylogenetic tree. Details are as follows.

* **Feature table:** It should contain counts, where rows are features (OTUs or ASVs) and columns are subjects (row names are feature IDs and column names are subject IDs). 

* **Taxonomic table:** It should contain taxonomic names, where rows are features and columns are seven taxonomic ranks (row names are feature IDs and column names are 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species' or 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'). 

* **Metadata/Sample information:** It should contain variables for the subjects about host phenotypes, medical interventions, disease status or environmental/behavioral factors, where rows are subjects and columns are variables (row names are subject IDs, and column names are variable names). 

* **Phylogenetic tree:** It should be a rooted tree. Otherwise, MiCML automatically roots the tree through midpoint rooting (phangorn::midpoint). The tip labels of the phylogenetic tree are feature IDs. 

**Notice:** The features should be matched and identical across feature table, taxonomic table and phylogenetic tree. The subjects should be matched and identical between feature table and metadata/sample information. MiCML will analyze only the matched features and subjects.

<hr>

## Example Data

**(1) Phyloseq:** You can download example microbiome data (Limeta et al., 2020) in a unified format, called phyloseq, 'biom.Rdata' in the directory 'Data/Phyloseq'. For more details about 'phyloseq', see https://bioconductor.org/packages/release/bioc/html/phyloseq.html
```
library(phyloseq)

setwd('/yourdatadirectory/')

load(file = 'biom.Rdata')

otu.tab <- otu_table(biom)
tax.tab <- tax_table(biom)
sam.dat <- sample_data(biom)
tree <- phy_tree(biom)
```
You can check if the features are matched and identical across feature table, taxonomic table and phylogenetic tree, and the subjects are matched and identical between feature table and metadata/sample information using following code.
```
identical(rownames(otu.tab), rownames(tax.tab))
identical(rownames(otu.tab), tree$tip.label)
identical(colnames(otu.tab), rownames(sam.dat))
```
**(2) Individual Data:** You can download example microbiome data (Limeta et al., 2020) as individual files: feature table (otu.tab.txt), taxonomic table (tax.tab.txt), metadata/sample information (sam.dat.txt), and phylogenetic tree (tree.tre) in the directory 'Data/Individual'
```
library(ape)

setwd('/yourdatadirectory/')

otu.tab <- read.table(file = 'otu.tab.txt', check.names = FALSE)
tax.tab <- read.table(file = 'tax.tab.txt', check.names = FALSE)
sam.dat <- read.table(file = 'sam.dat.txt', check.names = FALSE)
tree <- read.tree(file = 'tree.tre')
```
You can check if the features are matched and identical across feature table, taxonomic table and phylogenetic tree, and the subjects are matched and identical between feature table and metadata/sample information using following code.
```
identical(rownames(otu.tab), rownames(tax.tab))
identical(rownames(otu.tab), tree$tip.label)
identical(colnames(otu.tab), rownames(sam.dat))
```
* Reference: Limeta A, Ji B, Levin M, Gatto F, Nielsen J. Meta-analysis of the gut microbiota in predicting response to cancer immunotherapy in metastatic melanoma. JCL Insight. 2020;5(23):e140940.

<hr>

## Prerequites

**Notice:** For the local implementation, you do not need to install all the pre-requite R packages individually. You only need to install the 'shiny' package, and then run a simple command in 'Launch App' below. Then, all the pre-requisite R packages will be installed and imported automatically. For Mac users, please make sure to dowload and install xQuartz before the local implementation. (https://www.xquartz.org/)

shiny
```
install.packages('shiny')
```

<hr>

## Launch App

```
library(shiny)

runGitHub('MiCML', 'hk1785', ref = 'main')
```

<hr>

## Troubleshooting Tips

If you have any problems for using MiCML, please report in issues (https://github.com/hk1785/micml/issues) or email Hyunwook Koh (hyunwook.koh@stonybrook.edu).
