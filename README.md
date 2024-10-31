# MiCML

Title: Microbiome-based causal machine learning (MiCML) for the analysis of treatment effects using microbiome profiles

Version: 1.0.0

Maintainer: Hyunwook Koh <hyunwook.koh@stonybrook.edu>

Description: The treatment effects are heterogenous across patients due to the differences in their microbiomes, which in turn implies that we can enhance the treatment effect by manipulating the patient’s microbiome profile. Then, the coadministration of microbiome-based dietary supplements (e.g., prebiotics, probiotics, dietary fiber) or therapeutics (e.g., antibiotics, pharmabiotics, phage therapy, microbiota transplantation) along with the primary treatment (e.g., immunotherapy) has been the subject of intensive investigation. However, for this, we first need to comprehend which microbes help (or prevent) the treatment to cure the patient’s disease, which is in principle the matter of interaction effects between treatment and microbiome on the patient’s recovery. MiCML (microbiome causal machine learning) is the first cloud computing platform that streamlines related data processing and analytic procedures for the analysis of treatment effects using microbiome profiles on user-friendly web environments. MiCML is in particular unique with the up-to-date features of (i) batch effect correction to mitigate systematic variation in collective large-scale microbiome data due to the differences in their underlying batches (e.g., lab or study environments), and (ii) causal machine learning to estimate treatment effects with consistency and then discern microbial taxa that enhance (or lower) the efficacy of the primary treatment. We also stress that MiCML can handle the data from either randomized controlled trials or observational studies. MiCML can be a useful analytic tool for microbiome-based personalized medicine. MiCML consists of three Data Processing modules, (i) Data Input, (ii) Batch Effect Correction & Quality Control, and (iii) Data Transformation and three Data Analysis modules, (i) Descriptive Analysis, (ii) Generalized Linear Models, and (iii) Causal Machine Learning.

NeedsCompilation: No

Depends: R(≥ 4.1.0)

Imports: Bioconductor ('BiocParallel', 'biomformat', 'phyloseq'); CRAN ('grf', 'betareg', 'BiasedUrn', 'BiocManager', 'bios2mds', 'CompQuadForm', 'dashboardthemes', 'devtools', 'DiagrammeR', 'dirmult', 'dplyr', 'DT', 'ecodist', 'edarf', 'entropart', 'erer', 'fBasics', 'forestplot', 'fossil', 'ggplot2', 'ggthemes', 'googleVis', 'gridExtra', 'gridGraphics', 'compositions', 'GUniFrac', 'htmltools', 'ICSNP', 'lme4', 'lmerTest', 'MiRKAT', 'mmpf', 'nlme', 'patchwork', 'phangorn', 'picante', 'plotly', 'PMCMRplus', 'quantreg', 'remotes', 'reticulate', 'rgl', 'rmarkdown', 'robCompositions', 'robustbase', 'seqinr', 'shiny', 'shinydashboard', 'shinyjs', 'shinyWidgets', 'stringr', 'tidyverse', 'vegan', 'xtable', 'zCompositions', 'zip', 'bda', 'mediation'); GitHub ('ConQuR')

License: GPL 1, GPL 2 

**URLs**: Web Server (http://micml.micloud.kr), GitHub (http://github.com/hk1785/micml) 

**Maintainer**: Hyunwook Koh <hyunwook.koh@stonybrook.edu>

**Reference**: Koh H, Kim J, Jang H. MiCML: A causal machine learning cloud platform for the analysis of treatment effects using microbiome profiles (*In review*). 

## GitHub Repository Contents

* **Data** - In this directory, example microbiome data are stored in a widely used unified format, called phyloseq, that can efficiently combine all essential microbiome data components as well as individual files.

* **Source** - In this directory, all the R functions that are needed to run MiMedSurv are stored.

* **www** - In this directory, some photos that are used to decorate the GUI of MiMedSurv are stored.

* **app.R** - In this file, all the central codes to control for user-interfaces and server functions of MiMedSurv are stored.

## Prerequites

* Notice: For the local implementation, you do not need to install all the pre-requite R packages individually. You only need to install the 'shiny' package, and then run a simple command in 'Launch App' below. Then, all the pre-requisite R packages will be installed and imported automatically. For Mac users, please make sure to dowload and install xQuartz before the local implementation. (https://www.xquartz.org/)

shiny
```
install.packages('shiny')
```

## Launch App

```
library(shiny)

runGitHub('MiCML', 'hk1785', ref = 'main')
```

## Troubleshooting Tips

If you have any problems for using MiCML, please report in issues (https://github.com/hk1785/micml/issues) or email Hyunwook Koh (hyunwook.koh@stonybrook.edu).
