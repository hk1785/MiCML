options(scipen=999)
options(warn=-1)

ls.pkg <- c('ape', 'BiocManager', 'bios2mds', 'caret', 'checkmate', 'compositions', 'data.table', 'doParallel', 'DT', 'ecodist', 'edarf', 'fossil', 'fontawesome', 
            'GUniFrac', 'googleVis', 'ggplot2', 'ggplotify', 'grid', 'grf', 'htmltools', 
            'Matrix', 'MiRKAT', 'mmpf', 'phangorn', 'picante', 'plotly', 'proxy', 
            'randomForest', 'remotes', 'reshape2', 'rpart', 'rpart.plot', 'rmarkdown', 
            'seqinr', 'shiny', 'shinydashboard', 'shinyjs', 'shinyWidgets', 'stringr', 
            'tidyverse', 'vegan', 'VGAM', 'xtable', 'zCompositions', 'zip')

new.pkg <- ls.pkg[!(ls.pkg %in% installed.packages()[,"Package"])]
if(length(new.pkg)) install.packages(new.pkg, repos = 'https://cloud.r-project.org/')

if(!require('phyloseq')) BiocManager::install('phyloseq')
if(!require('biomformat')) remotes::install_github('joey711/biomformat')
if(!require('dashboardthemes')) remotes::install_github('nik01010/dashboardthemes', force = TRUE)
if(!require('chatgpt')) remotes::install_github('jcrodriguez1989/chatgpt')
if(!require('mmpf')) install.packages('Source/mmpf_0.0.5.tar.gz', repos = NULL, type="source")
if(!require('edarf')) remotes::install_github('zmjones/edarf', subdir = "pkg", force = TRUE)
if(!require('ConQuR')) remotes::install_github('wdl2459/ConQuR', force = TRUE)
if(!require('sva')) BiocManager::install('sva')
if(!require('aVirtualTwins')) remotes::install_github("prise6/aVirtualTwins", build_vignettes = TRUE)
if(!require('MiVT')) remotes::install_github("hk1785/MiVT", force = TRUE)
if(!require('MiRKATMC')) remotes::install_github("Zhiwen-Owen-Jiang/MiRKATMC")

library(aVirtualTwins)
library(ape)
library(bios2mds)
library(BiocManager)
library(biomformat)
library(chatgpt)
library(ConQuR)
library(compositions)
library(caret)
library(doParallel)
library(dashboardthemes)
library(DT)
library(data.table)
library(dplyr)
library(ecodist)
library(edarf)
library(fontawesome)
library(fossil)
library(googleVis)
library(ggplot2)
library(grid)
library(grf)
library(ggplotify)
library(GUniFrac)
library(htmltools)
library(Matrix)
library(MiRKAT)
library(MiRKATMC)
library(MiVT)
library(mmpf)
library(plotly)
library(phangorn)
library(phyloseq)
library(proxy)
library(picante)
library(remotes)
library(rpart)
library(rpart.plot)
library(randomForest)
library(reshape2)
library(stringr)
library(shiny)
library(seqinr)
library(shinydashboard)
library(shinyWidgets)
library(shinyjs)
library(tidyverse)
library(vegan)
library(VGAM)
library(xtable)
library(zCompositions)
library(zip)

source("Source/MiDataProc.DataInput.R")
source("Source/MiDataProc.Descriptive.R")
source("Source/MiDataProc.GLM.R")
source("Source/MiDataProc.CML.R")

# COMMENTS ------

{
  TITLE = p("MiCML: Microbiome Causal Machine Learning for the Analysis of Treatment Effects Using Microbiome Profiles", style = "font-size:16pt")
  
  HOME_COMMENT_MV = p(strong("Motivation:"), "The treatment effects are heterogenous across patients due to the differences in their microbiomes, 
  which in turn implies that we can enhance the treatment effect by manipulating the patient’s microbiome profile. 
  Then, the coadministration of microbiome-based dietary supplements/therapeutics along with the primary treatment 
  has been the subject of intensive investigation. However, for this, we first need to comprehend which microbes help (or prevent) 
  the treatment to cure the patient’s disease. which is in principle the matter of interaction effects between treatment and microbiome on the patient’s recovery.", style = "font-size:12pt")
  
  HOME_COMMENT = p(strong("MiCML (Microbiome Causal Machine Learning)"), "is a cloud computing platform that streamlines related data processing and analytic procedures for the analysis of treatment effects 
  using microbiome profiles on user-friendly web environments. MiCML is in particular unique with the up-to-date features of (i) ", 
                   strong("batch effect correction"), "to mitigate systematic variation in collective large-scale microbiome data due to the differences in their underlying batches (e.g., lab or study environments), and (ii) ", 
                   strong("causal machine learning"), "to estimate treatment effects with consistency and then discern microbial taxa that enhance (or lower) the efficacy of the primary treatment. 
                   We also stress that MiCML can handle the data from either randomized controlled trials or observational studies. MiCML can be a useful analytic tool for microbiome-based personalized medicine. 
                   MiCML consists of three Data Processing modules, (i) Data Input, (ii) Batch Effect Correction & Quality Control, and (iii) Data Transformation 
                   and three Data Analysis modules, (i) Descriptive Analysis, (ii) Generalized Linear Models, and (iii) Causal Machine Learning.", style = "font-size:12pt")
  
  HOME_COMMENT2 = p(strong("URLs:"), "Web server (online implementation):", tags$a(href = "http://micml.micloud.kr", "http://micml.micloud.kr"), 
                    "; GitHub repository (local implementation):", 
                    tags$a(href = "https://github.com/hk1785/micml", "https://github.com/hk1785/micml"), style = "font-size:12pt")
  
  HOME_COMMENT3 = p(strong("Maintainers:"), "Hyunwook Koh (", tags$a(href = "hyunwook.koh@stonybrook.edu", "hyunwook.koh@stonybrook.edu"), ")", style = "font-size:12pt")
  
  HOME_COMMENT4 = p(strong("Reference:"), "Koh H, Kim J, Jang H. MiCML: A causal machine learning cloud platform for the analysis of treatment effects using microbiome profiles (In review)" 
                    , style = "font-size:12pt")
  
  INPUT_PHYLOSEQ_COMMENT1 = p(strong("Description:"), br(), br(), "This should be an '.rdata' or '.rds' file, and the data should be in 'phyloseq' format (see ", 
                              htmltools::a(tags$u("https://bioconductor.org/packages/release/bioc/html/phyloseq.html"), style = "color:red3"),
                              "). The phyloseq object should contain three necessary data, feature (OTU or ASV) table, taxonomic table and meta/sample information.",br(), br(),
                              strong("Details:"), br(), br(), 
                              strong("Feature table:"), "It should contain counts, where rows are features (OTUs or ASVs) and columns are subjects (row names are feature IDs and column names are subject IDs).", br(), br(),
                              strong("Taxonomic table:"),"It should contain taxonomic names, where rows are features and columns are seven taxonomic ranks (row names are feature IDs and column names are 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species' or 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species').", br(), br(),
                              strong("Metadata/Sample information:"),"It should contain variables for the subjects about host phenotypes, medical interventions, disease status or environmental/behavioral factors, where rows are subjects and columns are variables (row names are subject IDs, and column names are variable names).", br(), br(),
                              strong("Phylogenetic tree:"),"It should be a rooted tree. Otherwise, MiCML automatically roots the tree through midpoint rooting (phangorn::midpoint). The tip labels of the phylogenetic tree are feature IDs.", br(), br(),
                              "* The features should be matched and identical across feature table, taxonomic table and phylogenetic tree. 
                              The subjects should be matched and identical between feature table and metadata/sample information. MiCML will analyze only the matched features and subjects.", style = "font-size:11pt")
  
  INPUT_PHYLOSEQ_COMMENT2 = p("You can download example microbiome data 'biom.Rdata' in 'phyloseq' format. For more details about 'phyloseq', see ", 
                              htmltools::a(tags$u("https://bioconductor.org/packages/release/bioc/html/phyloseq.html"), style = "color:red3"), br(), br(), 
                              "> library(phyloseq)", br(), br(), 
                              "> setwd('/yourdatadirectory/')", br(), br(), 
                              "> load(file = 'biom.Rdata')", br(), br(), 
                              "> otu.tab <- otu_table(biom)", br(), 
                              "> tax.tab <- tax_table(biom)", br(), 
                              "> sam.dat <- sample_data(biom)", br(), 
                              "> tree <- phy_tree(biom)", br(),br(), 
                              "You can check if the features are matched and identical across feature table, taxonomic table and phylogenetic tree, and the subjects are matched and identical between feature table and metadata/sample information using following code.", br(), br(), 
                              "> identical(rownames(otu.tab), rownames(tax.tab))", br(), 
                              "> identical(rownames(otu.tab), tree$tip.label)", br(),
                              "> identical(colnames(otu.tab), rownames(sam.dat))", style = "font-size:11pt", br(), br(),
                              strong("Reference:"), "Limeta A, Ji B, Levin M, Gatto F, Nielsen J. Meta-analysis of the gut microbiota in predicting response to cancer immunotherapy in metastatic melanoma. JCL Insight. 2020;5(23):e140940.")
  
  INPUT_INDIVIDUAL_DATA_COMMENT = p(strong("Description:"), br(), br(), 
                                    strong("Feature table:"), "It should contain counts, where rows are features (OTUs or ASVs) and columns are subjects (row names are feature IDs and column names are subject IDs).", br(), br(),
                                    strong("Taxonomic table:"),"It should contain taxonomic names, where rows are features and columns are seven taxonomic ranks (row names are feature IDs and column names are 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species' or 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species').", br(), br(),
                                    strong("Metadata/Sample information:"),"It should contain variables for the subjects about host phenotypes, medical interventions, disease status or environmental/behavioral factors, where rows are subjects and columns are variables (row names are subject IDs, and column names are variable names).", br(), br(),
                                    strong("Phylogenetic tree:"),"It should be a rooted tree. Otherwise, MiCML automatically roots the tree through midpoint rooting (phangorn::midpoint). The tip labels of the phylogenetic tree are feature IDs.", br(), br(),
                                    "* The features should be matched and identical across feature table, taxonomic table and phylogenetic tree.
                                    The subjects should be matched and identical between feature table and metadata/sample information. MiCML will analyze only the matched features and subjects.", style = "font-size:11pt")
  
  INPUT_INDIVIDUAL_DATA_COMMENT2 = p("You can download example microbiome data 'Individual.zip'. This zip file contains four necessary data components: feature table (otu.tab.txt), taxonomic table (tax.tab.txt), metadata/sample information (sam.dat.txt), and phylogenetic tree (tree.tre).", br(), br(),
                                     "> library(ape)", br(), br(), 
                                     "> setwd('/yourdatadirectory/')", br(), br(), 
                                     "> otu.tab <- read.table(file = 'otu.tab.txt', check.names = FALSE)", br(), 
                                     "> tax.tab <- read.table(file = 'tax.tab.txt', check.names = FALSE)", br(), 
                                     "> sam.dat <- read.table(file = 'sam.dat.txt', check.names = FALSE)", br(), 
                                     "> tree <- read.tree(file = 'tree.tre)", br(),br(),
                                     "You can check if the features are matched and identical across feature table, taxonomic table and phylogenetic tree, and the subjects are matched and identical between feature table and metadata/sample information using following code.", br(), br(), 
                                     " > identical(rownames(otu.tab), rownames(tax.tab))", br(), 
                                     " > identical(rownames(otu.tab), tree$tip.label)", br(),
                                     " > identical(colnames(otu.tab), rownames(sam.dat))", style = "font-size:11pt", br(), br(),
                                     strong("Reference:"), "Limeta A, Ji B, Levin M, Gatto F, Nielsen J. Meta-analysis of the gut microbiota in predicting response to cancer immunotherapy in metastatic melanoma. JCL Insight. 2020;5(23):e140940.")
  
  QC_KINGDOM_COMMENT = p("A microbial kingdom to be analyzed. Default is 'Bacteria' for 16S data. Alternatively, you can type 'Fungi' for ITS data 
                         or any other kingdom of interest for shotgun metagenomic data.", 
                         style = "font-size:11pt")
  
  QC_LIBRARY_SIZE_COMMENT1 = p("Remove subjects that have low library sizes (total read counts). Default is 3,000.", 
                               style = "font-size:11pt")
  
  QC_LIBRARY_SIZE_COMMENT2 = p("Library size: The total read count per subject.", 
                               style = "font-size:11pt")
  
  QC_MEAN_PROP_COMMENT1 = p("Remove features (OTUs or ASVs) that have low mean relative abundances (Unit: %). Default is 0.02%.",
                            style = "font-size:11pt")
  
  QC_MEAN_PROP_COMMENT2 = p("Mean proportion: The average of relative abundances (i.e., proportions) per feature.", 
                            style = "font-size:11pt")
  
  QC_TAXA_NAME_COMMENT1 = p("Remove taxonomic names in the taxonomic table that are completely matched with the specified character strings. 
                            Multiple character strings should be separated by a comma. 
                            Default is \"\", \"metagenome\", \"gut metagenome\", \"mouse gut metagenome\".",
                            style = "font-size:11pt")
  QC_TAXA_NAME_COMMENT2 = p("Remove taxonomic names in the taxonomic table that are partially matched with the specified character strings 
                            (i.e., taxonomic names that contain the specified character strings). Multiple character strings should be separated by a comma. 
                            Default is \"uncultured\", \"incertae\", \"Incertae\", \"unidentified\", \"unclassified\", \"unknown\".",
                            style = "font-size:11pt")
  
  QC_BATCH_REFERENCE = p("1. Ling W, Lu J, Zhao N. et al. Batch effects removal for microbiome data via conditional quantile regression. Nat Commun. 2022;13(5418)", br(),
                         "2. Torgerson WS. Multidimensional scaling: I. Theorey and method. Psychometrika. 1952;17:401-419.", br(),
                         "3. Bray JR, Curtis JT. An ordination of the upland forest communities of southern Wisconsin. Ecol Monogr. 1957;27:325-349.", br(),
                         "4. Aitchison J. The statistical analysis of compositional data. J R Stat Soc Series B Stat Methodol. 1982;44(2):139-160.")
  
  QC_BATCH_CONQUR_REFERENCE = p("Ling W, Lu J, Zhao N. et al. Batch effects removal for microbiome data via conditional quantile regression. Nat Commun. 2022;13(5418).", style = "font-size:11pt")
  
  QC_BATCH_COMBAT_REFERENCE = p("Zhang Y. et al. ComBat-seq: batch effect adjustments for RNA-seq count data. NAR Genom Bioinform. 2020;2(3)lqaa078", style = "font-size:11pt")
  
  DATA_TRANSFORM_COMMENT = p("Transform the data into four different formats (1) CLR (centered log ratio) (Aitchison, 1982), (2) Count (Rarefied) (Sanders, 1968), (3) Proportion, (4) Arcsine-root 
                             for each taxonomic rank (phylum, class, order, familiy, genus, species).")
  
  DATA_TRANSFORM_REFERENCE = p("1. Aitchison J. The statistical analysis of compositional data. J R Stat Soc B. 1982;44(2):139-77.", br(),
                               "2. Sanders HL. Marine benthic diversity: A comparative study. Am Nat. 1968;102:243-282.")
  
  MANN_WHITNEY_REFERENCE = p("Mann HB, Whitney DR. On a test of whether one of two random variables is stochastically larger than the other. Ann Math Stat. 1947;18(1):50-60.")
  
  PEARSON_REFERENCE = p("Pearson K. On the criterion that a given system of deviations from the probable in the case of a correlated system of variables is such that I can be reasonably supposed to have arisen from random sampling. Philos Mag, Series 5. 1900;50(302):157-175.")
  
  WELCH_REFERENCE = p("Welch BL. The generalization of Student’s problem when several different population variances are involved. Biometrika. 1947;34(1-2):28-35.")
  
  FISHER_REFERENCE = p("Fisher RA. On the interpretation of Chi-square from contingency tables, and the calculation of P. J R Stat Soc. 1922;85(1):87-94.")
  
  GLM_REFERENCE_CLR = p("1. Aitchison J. The statistical analysis of compositional data. J R Stat Soc B. 1982;44(2):139-77.", br(), 
                        "2. Nelder J, Wedderburn R. Generalized linear models. J R Stat Soc A. 1972;135(3):370-384.", br(), 
                        "3. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B. 1995;57:289-300.")
  
  GLM_REFERENCE_NO_CLR = p("1. Nelder J, Wedderburn R. Generalized linear models. J R Stat Soc A. 1972;135(3):370-384.", br(), 
                           "2. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B. 1995;57:289-300.")
  
  LM_REFERENCE_CLR = p("1. Aitchison J. The statistical analysis of compositional data. J R Stat Soc B. 1982;44(2):139-77.", br(), 
                       "2. Rosenblatt M. Remarks on some nonparametric estimates of a density function. Ann Math Stat. 1956;27(3):832-837.", br(), 
                       "3. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B. 1995;57:289-300.")
  
  LM_REFERENCE_NO_CLR = p("1. Rosenblatt M. Remarks on some nonparametric estimates of a density function. Ann Math Stat. 1956;27(3):832-837.", br(), 
                          "2. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B. 1995;57:289-300.")
  
  DST_REFERENCE_CLR = p("1. Aitchison J. The statistical analysis of compositional data. J R Stat Soc B. 1982;44(2):139-77.", br(),
                        "2. Breiman L, Friedman JH, Olshen RA, Stone CJ. Classification and regression trees. CRC Press. 1984.", br(),
                        "3. Breiman L. Random forests. Mach Learn. 2001;45:5-32.", br(),
                        "4. Foster JC, Taylor JM, Ruberg SJ. Subgroup identification from randomized clinical trial data. Stat Med. 2011;30(24):2867-2880.", br(),
                        "5. Wager S, Athey S. Estimation and inference of heterogeneous treatment effects using random forests. J Am Stat Assoc. 2018;113(523):1228-42.")
  
  DST_REFERENCE_RARE = p("1. Sanders HL. Marine benthic diversity: A comparative study. Am Nat. 1968;102:243-282.", br(),
                         "2. Breiman L, Friedman JH, Olshen RA, Stone CJ. Classification and regression trees. CRC Press. 1984.", br(),
                         "3. Breiman L. Random forests. Mach Learn. 2001;45:5-32.", br(),
                         "4. Foster JC, Taylor JM, Ruberg SJ. Subgroup identification from randomized clinical trial data. Stat Med. 2011;30(24):2867-2880.", br(),
                         "5. Wager S, Athey S. Estimation and inference of heterogeneous treatment effects using random forests. J Am Stat Assoc. 2018;113(523):1228-42.")
  
  DST_REFERENCE = p("1. Breiman L, Friedman JH, Olshen RA, Stone CJ. Classification and regression trees. CRC Press. 1984.", br(),
                    "2. Breiman L. Random forests. Mach Learn. 2001;45:5-32.", br(),
                    "3. Foster JC, Taylor JM, Ruberg SJ. Subgroup identification from randomized clinical trial data. Stat Med. 2011;30(24):2867-2880.", br(),
                    "4. Wager S, Athey S. Estimation and inference of heterogeneous treatment effects using random forests. J Am Stat Assoc. 2018;113(523):1228-42.")
  
  PT_REFERENCE_CLR = p("1. Aitchison J. The statistical analysis of compositional data. J R Stat Soc B. 1982;44(2):139-77.", br(),
                       "2. Breiman L, Friedman JH, Olshen RA, Stone CJ. Classification and regression trees. CRC Press. 1984.", br(),
                       "3. Breiman L. Random forests. Mach Learn. 2001;45:5-32.", br(),
                       "4. Hirano K, Imbens GW, Rider G. Efficient estimation of average treatment effects using the estimated propensity score. Econometrica. 2003;71:1161-1189.")
  
  PT_REFERENCE_RARE = p("1. Sanders HL. Marine benthic diversity: A comparative study. Am Nat. 1968;102:243-282.", br(),
                        "2. Breiman L, Friedman JH, Olshen RA, Stone CJ. Classification and regression trees. CRC Press. 1984.", br(),
                        "3. Breiman L. Random forests. Mach Learn. 2001;45:5-32.", br(),
                        "4. Hirano K, Imbens GW, Rider G. Efficient estimation of average treatment effects using the estimated propensity score. Econometrica. 2003;71:1161-1189.")
  
  PT_REFERENCE = p("1. Breiman L, Friedman JH, Olshen RA, Stone CJ. Classification and regression trees. CRC Press. 1984.", br(),
                   "2. Breiman L. Random forests. Mach Learn. 2001;45:5-32.", br(),
                   "3. Foster JC, Taylor JM, Ruberg SJ. Subgroup identification from randomized clinical trial data. Stat Med. 2011;30(24):2867-2880.", br(),
                   "4. Hirano K, Imbens GW, Rider G. Efficient estimation of average treatment effects using the estimated propensity score. Econometrica. 2003;71:1161-1189.")
}

# UI ---------------------------------------------------------------------------

{
  ui = dashboardPage(
    title = "MiCML",
    dashboardHeader(title = span(TITLE, style = "float:left;font-size: 20px"), titleWidth = "100%"),
    dashboardSidebar(
      tags$script(JS("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';")),
      sidebarMenu(id = "side_menu",
                  
                  menuItem("Home", tabName = "home"),
                  
                  menuItem("Data Processing",
                           menuSubItem("Data Input", tabName = "step1", 
                                       icon = fontawesome::fa("upload", margin_left = "0.3em", margin_right = "0.1em")),
                           menuSubItem(span(span("Batch Effect Correction /", br()), 
                                            span("Quality Control", style = "margin-left: 23px")), tabName = "step2", 
                                       icon = fontawesome::fa("chart-bar", margin_left = "0.3em")),
                           menuSubItem("Data Transformation", tabName = "divdtCalculation", 
                                       icon = fontawesome::fa("calculator", margin_left = "0.3em", margin_right = "0.25em"))
                  ),
                  
                  menuItem("Data Analysis",
                           menuSubItem("Descriptive", tabName = "descriptive", 
                                       icon = fontawesome::fa("chart-simple", margin_right = "0.1em")),
                           menuSubItem("Generalized Linear Models", tabName = "glm", 
                                       icon = fontawesome::fa("chart-line", margin_right = "0em")),
                           menuSubItem("Causal Machine Learning", tabName = "causal_forest", 
                                       icon = fontawesome::fa("tree", margin_right = "0.1em"))
                  )
      )
    ),
    
    dashboardBody(
      
      # Custom Theme -----
      
      shinyDashboardThemeDIY(
        
        # general
        appFontFamily = "Arial"
        ,appFontColor = "rgb(0,0,0)"
        ,primaryFontColor = "rgb(0,0,0)"
        ,infoFontColor = "rgb(0,0,0)"
        ,successFontColor = "rgb(0,0,0)"
        ,warningFontColor = "rgb(0,0,0)"
        ,dangerFontColor = "rgb(0,0,0)"
        ,bodyBackColor = "rgb(255,255,255)"
        
        # header
        ,logoBackColor = "rgb(250, 60, 70)"
        
        ,headerButtonBackColor = "rgb(250, 60, 70)"
        ,headerButtonIconColor = "rgb(250, 60, 70)"
        ,headerButtonBackColorHover = "rgb(250, 60, 70)"
        ,headerButtonIconColorHover = "rgb(0,0,0)"
        
        ,headerBackColor = "rgb(250, 60, 70)"
        ,headerBoxShadowColor = "#aaaaaa"
        ,headerBoxShadowSize = "0px 0px 0px"
        
        # sidebar
        ,sidebarBackColor = "rgb(24,31,41)"
        ,sidebarPadding = 0
        
        ,sidebarMenuBackColor = "transparent"
        ,sidebarMenuPadding = 0
        ,sidebarMenuBorderRadius = 0
        
        ,sidebarShadowRadius = ""
        ,sidebarShadowColor = "0px 0px 0px"
        
        ,sidebarUserTextColor = "rgb(24,31,41)"
        
        ,sidebarSearchBackColor = "rgb(255, 255, 255)"
        ,sidebarSearchIconColor = "rgb(24,31,41)"
        ,sidebarSearchBorderColor = "rgb(24,31,41)"
        
        ,sidebarTabTextColor = "rgb(210,210,210)"
        ,sidebarTabTextSize = 14
        ,sidebarTabBorderStyle = "none"
        ,sidebarTabBorderColor = "none"
        ,sidebarTabBorderWidth = 0
        
        ,sidebarTabBackColorSelected = "rgb(45,52,63)"
        ,sidebarTabTextColorSelected = "rgb(252,255,255)"
        ,sidebarTabRadiusSelected = "0px"
        
        ,sidebarTabBackColorHover = "rgb(67,75,86)"
        ,sidebarTabTextColorHover = "rgb(252,255,255)"
        ,sidebarTabBorderStyleHover = "none"
        ,sidebarTabBorderColorHover = "none"
        ,sidebarTabBorderWidthHover = 0
        ,sidebarTabRadiusHover = "0px"
        
        # boxes
        ,boxBackColor = "rgb(245,245,245)"
        ,boxBorderRadius = 3
        ,boxShadowSize = "0px 0px 0px"
        ,boxShadowColor = "rgba(0,0,0,0)"
        ,boxTitleSize = 16
        ,boxDefaultColor = "rgb(210,214,220)"
        ,boxPrimaryColor = "rgb(35, 49, 64)"
        ,boxInfoColor = "rgb(250, 60, 70)"
        ,boxSuccessColor = "rgb(112,173,71)"
        ,boxWarningColor = "rgb(244,156,104)"
        ,boxDangerColor = "rgb(255,88,55)"
        
        ,tabBoxTabColor = "rgb(255,255,255)"
        ,tabBoxTabTextSize = 14
        ,tabBoxTabTextColor = "rgb(0,0,0)"
        ,tabBoxTabTextColorSelected = "rgb(35, 49, 64)"
        ,tabBoxBackColor = "rgb(255,255,255)"
        ,tabBoxHighlightColor = "rgb(250, 60, 70)"
        ,tabBoxBorderRadius = 0
        
        # inputs
        ,buttonBackColor = "rgb(245,245,245)"
        ,buttonTextColor = "rgb(0,0,0)"
        ,buttonBorderColor = "rgb(24,31,41)"
        ,buttonBorderRadius = 3
        
        ,buttonBackColorHover = "rgb(227,227,227)"
        ,buttonTextColorHover = "rgb(100,100,100)"
        ,buttonBorderColorHover = "rgb(200,200,200)"
        
        ,textboxBackColor = "rgb(255,255,255)"
        ,textboxBorderColor = "rgb(200,200,200)"
        ,textboxBorderRadius = 0
        ,textboxBackColorSelect = "rgb(245,245,245)"
        ,textboxBorderColorSelect = "rgb(200,200,200)"
        
        # tables
        ,tableBackColor = "rgb(255, 255, 255)"
        ,tableBorderColor = "rgb(245, 245, 245)"
        ,tableBorderTopSize = 1
        ,tableBorderRowSize = 1
        
      ),
      
      # Contents -----
      
      tags$head(tags$style(HTML(".content { padding-top: 2px;}"))),
      
      # Progress bar -----
      
      tags$head(tags$style(HTML('.progress-bar {background-color: rgb(22, 48, 91);}'))),
      
      # Pretty Radio Button -----
      
      tags$head(tags$style(HTML('.pretty input:checked~.state.p-primary label:after, .pretty.p-toggle .state.p-primary label:after {background-color: rgb(22, 48, 91)!important;}'))),
      
      # Slider -----
      
      setSliderColor(rep("rgb(22, 48, 91)", 100), seq(1, 100)),
      chooseSliderSkin("Flat"),
      
      # Tab Panel -----
      
      tags$style(HTML(".tabbable > .nav > li > a {color: #93063e}")),
      
      tags$script(src = "fileInput_text.js"),
      useShinyjs(),
      tabItems(
        
        ## Home -----
        
        tabItem(tabName = "home",
                div(id = "homepage", br(), HOME_COMMENT_MV, HOME_COMMENT, br(),
                    div(tags$img(src="home.png", height = 250), style = "text-align: center;"), br(),
                    HOME_COMMENT2, HOME_COMMENT3, HOME_COMMENT4)),
        
        ## 0. DATA INPUT -----
        
        tabItem(tabName = "step1", br(),
                fluidRow(column(width = 6,
                                box(width = NULL, status = "info", solidHeader = TRUE,
                                    title = strong("Data Input", style = "color:white"),
                                    selectInput("inputOption", h4(strong("Data type")), 
                                                choices = c("Choose one" = "", "Phyloseq", "Individual Data"), width = '30%'),
                                    div(id = "optionsInfo", 
                                        tags$p("You can choose phyloseq or individual data.", style = "font-size:11pt"), 
                                        tags$p("", style = "margin-bottom:-8px"), style = "margin-top: -15px"),
                                    uiOutput("moreOptions")
                                )
                ), column(width = 6, style='padding-left:0px', uiOutput("addDownloadinfo"))
                )
        ),
        
        ## 1-1. QC ----
        
        tabItem(tabName = "step2", br(),
                fluidRow(column(width = 3,  style = "padding-left:+15px",
                                
                                # Batch Effect Correction
                                
                                box(
                                  width = NULL, status = "info", solidHeader = TRUE, 
                                  title = strong("Batch Effect Correction", style = "color:white"),
                                  h4(strong("Batch effect correction?", style = "color:black")),
                                  p("Do you want to perform batch effect correction to balance the microbiome data across batches (e.g., labs, studies, locations, times) 
                                    while preserving the signals from other importanct primary and nuisance variables? (optional)", style = "font-size:10pt"),
                                  prettyRadioButtons("batch.yn", label = NULL,
                                                     animation = "jelly",c("No", "Yes"), selected = "No", 
                                                     icon = icon("check"), width = '80%'),
                                  uiOutput("batch.method.select"),
                                  uiOutput("batch"), 
                                  uiOutput("prim"), 
                                  uiOutput("covar")
                                ),
                                
                                # Quality Control
                                
                                box(
                                  width = NULL, status = "info", solidHeader = TRUE,
                                  title = strong("Quality Control", style = "color:white"),
                                  textInput("kingdom", h4(strong("Kingdom")), value = "Bacteria"),
                                  QC_KINGDOM_COMMENT,
                                  tags$style(type = 'text/css', '#slider1 .irs-grid-text {font-size: 1px}'),
                                  tags$style(type = 'text/css', '#slider2 .irs-grid-text {font-size: 1px}'),
                                  
                                  sliderInput("slider1", h4(strong("Library size")), 
                                              min=0, max=10000, value = 3000, step = 1000),
                                  QC_LIBRARY_SIZE_COMMENT1,
                                  QC_LIBRARY_SIZE_COMMENT2,
                                  
                                  sliderInput("slider2", h4(strong("Mean proportion")), 
                                              min = 0, max = 0.1, value = 0.02, step = 0.001,  post  = " %"),
                                  QC_MEAN_PROP_COMMENT1,
                                  QC_MEAN_PROP_COMMENT2,
                                  
                                  br(),
                                  p(" ", style = "margin-bottom: -20px;"),
                                  
                                  h4(strong("Errors in taxonomic names")),
                                  textInput("rem.str", label = "Complete match", value = ""),
                                  QC_TAXA_NAME_COMMENT1,
                                  
                                  textInput("part.rem.str", label = "Partial match", value = ""),
                                  QC_TAXA_NAME_COMMENT2,
                                  
                                  actionButton("run", (strong("Run!")), class = "btn-info"), 
                                  p(" ", style = "margin-bottom: +10px;"), 
                                  p(strong("Attention:"),"You have to click this Run button to perform 
                                    data transformation and further analyses.", style = "margin-bottom:-10px"), br()
                                ),
                                uiOutput("moreControls"),
                                uiOutput("qc_reference")
                ),
                
                column(width = 9, style = "padding-left:+10px",
                       box(
                         width = NULL, status = "info", solidHeader = TRUE,
                         fluidRow(width = 12,
                                  status = "info", solidHeader = TRUE, 
                                  valueBoxOutput("sample_Size", width = 3),
                                  valueBoxOutput("OTUs_Size", width = 3),
                                  valueBoxOutput("phyla", width = 3),
                                  valueBoxOutput("classes", width = 3)
                         ),
                         fluidRow(width = 12, 
                                  status = "info", solidHeader = TRUE,
                                  valueBoxOutput("orders", width = 3),
                                  valueBoxOutput("families", width = 3),
                                  valueBoxOutput("genera", width = 3),
                                  valueBoxOutput("species", width = 3)
                         ),
                         fluidRow(style = "position:relative",
                                  tabBox(width = 6, title = strong("Library Size", style = "color:black"), 
                                         tabPanel("Histogram",
                                                  plotlyOutput("hist"),
                                                  sliderInput("binwidth", "# of Bins:",
                                                              min = 0, max = 100, value = 50, width = "100%")),
                                         tabPanel("Box Plot", 
                                                  plotlyOutput("boxplot"))),
                                  tabBox(width = 6, title = strong("Mean Proportion", style = "color:black"), 
                                         tabPanel("Histogram",
                                                  plotlyOutput("hist2"),
                                                  sliderInput("binwidth2", "# of Bins:",
                                                              min = 0, max = 100, value = 50, width = "100%")),
                                         tabPanel("Box Plot",
                                                  plotlyOutput("boxplot2")))
                         )
                       ),
                       shinyjs::hidden(
                         shiny::div(id = "pcoa.area",
                                    box(width = NULL, status = "info", solidHeader = TRUE,
                                        title = strong("Batch Effect Correction", style = "color:white"),
                                        plotOutput("batch.pcoa", height = 600)))
                       )
                )
                )
        ),
        
        ## 1-2. Data Transformation -----
        
        tabItem(tabName = "divdtCalculation", br(),
                fluidRow(column(width = 6, style = "padding-left:+15px",
                                box(title = strong("Data Transformation", style = "color:white"), 
                                    width = NULL, status = "info", solidHeader = TRUE, DATA_TRANSFORM_COMMENT, 
                                    actionButton("dtRun", (strong("Run!")), class = "btn-info")),
                                uiOutput("dtDownload")
                ),
                
                column(width = 6, style='padding-left:0px',
                       box(title = strong("References", style = "color:white"), 
                           width = NULL, status = "info", solidHeader = TRUE, DATA_TRANSFORM_REFERENCE)
                )
                )
        ),
        
        ##2. Descriptive ------
        
        tabItem(tabName = "descriptive", br(),
                sidebarLayout(
                  sidebarPanel(width = 3,
                               uiOutput("de_response"),
                               p(" ", style = "margin-bottom: +10px;"),
                               uiOutput("de_response_rename"),
                               p(" ", style = "margin-bottom: +10px;"),
                               uiOutput("de_treat"),
                               uiOutput("de_treat_rename"),
                               uiOutput("de_chooseTest"),
                               p(" ", style = "margin-bottom: -10px;"),
                               uiOutput("de_legend"),
                               p(" ", style = "margin-bottom: -10px;"),
                               uiOutput("de_download"),
                               p(" ", style = "margin-bottom: +30px;"),
                               uiOutput("de_reference")),
                  mainPanel(width = 9,
                            fluidPage(width = NULL,
                                      uiOutput("de_display"))
                  ))),
        
        ## 3. Generalized Linear Models -----
        
        tabItem(tabName = "glm", br(),
                sidebarLayout(
                  sidebarPanel(width = 3,
                               uiOutput("int_format"),
                               uiOutput("int_response"),
                               p(" ", style = "margin-bottom: +10px;"),
                               uiOutput("int_response_rename"),
                               p(" ", style = "margin-bottom: +20px;"),
                               uiOutput("int_treat"),
                               uiOutput("int_treat_rename"),
                               uiOutput("int_cov"),
                               uiOutput("int_chooseTest"),
                               uiOutput("int_num_taxa"),
                               uiOutput("int_legend"),
                               p(" ", style = "margin-bottom: -10px;"),
                               uiOutput("int_download"),
                               p(" ", style = "margin-bottom: +30px;"),
                               uiOutput("int_reference")),
                  mainPanel(width = 9,
                            fluidPage(width = NULL,
                                      uiOutput("int_display"))
                  ))),
        
        ## 4. Causal Forest -----
        
        tabItem(tabName = "causal_forest", br(),
                sidebarLayout(
                  sidebarPanel(width = 3,
                               shinyjs::hidden(
                                 uiOutput("cf_format"),
                                 uiOutput("cf_response"),
                                 uiOutput("cf_response_rename"),
                                 uiOutput("cf_treat"),
                                 uiOutput("cf_treat_rename"),
                                 uiOutput('cf_method'),
                                 uiOutput("cf_cov"),
                                 uiOutput("cf_tax_rank"),
                                 uiOutput("cf_num_taxa"),
                                 uiOutput("cf_download"),
                                 uiOutput("cf_reference")
                               )),
                  mainPanel(width = 9,
                            fluidPage(width = NULL,
                                      uiOutput("cf_display"))
                  )))
      )
    )
  )
}

# SERVER -----------------------------------------------------------------------

server = function(input, output, session){
  options(shiny.maxRequestSize = 30*1024^2)
  
  env <- new.env()
  nm <- load(file = "Data/Phyloseq/biom.Rdata", env)[1]
  biom <- env[[nm]]
  
  otu.tab <- read.table(file = "Data/Individual/otu.tab.txt", header = TRUE, check.names = FALSE, sep = "\t")
  tax.tab <- read.table(file = "Data/Individual/tax.tab.txt", header = TRUE, check.names = FALSE, sep = "\t")
  sam.dat <- read.table(file = "Data/Individual/sam.dat.txt", header = TRUE, check.names = FALSE, sep = "\t")
  tree <- read.tree(file = "Data/Individual/tree.tre")
  
  output$downloadData.Immuno <- downloadHandler(
    filename = function() {
      paste("biom.Rdata", sep = "")
    },
    content = function(file1) {
      save(biom, file = file1)
    })
  output$downloadZip.Immuno <- downloadHandler(
    filename = function() {
      paste("individual",".zip", sep = "")
    },
    content <- function(fname) {
      temp <- setwd(tempdir())
      on.exit(setwd(temp))
      dataFiles = c("otu.tab.txt", "tax.tab.txt", "sam.dat.txt", "tree.tre")
      write.table(otu.tab, "otu.tab.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.table(tax.tab, "tax.tab.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.table(sam.dat, "sam.dat.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.tree(tree, "tree.tre")
      zip(zipfile=fname, files=dataFiles)
    })
  
  ## Reactive Variables -------
  infile = reactiveValues(biom = NULL, qc_biom = NULL, rare_biom = NULL, na_omit_biom = NULL)
  batch.infile = reactiveValues(batchid = NULL, bat_biom = NULL, qc_biom = NULL, rare_biom = NULL)
  chooseData = reactiveValues(sam.dat = NULL, mon.sin.rev.bin.con = NULL, prim_vars = NULL, alpha.div = NULL,
                              alpha.div.rare = NULL, alpha.div.qc = NULL, taxa.out = NULL, tax.tab = NULL)
  ds.Ks <- reactiveValues(res = NULL)
  taxa.results = reactiveValues(bin.var = NULL, cov.var = NULL, id.var = NULL, taxa = NULL, taxa.bin.sum.out = NULL,
                                con.var = NULL, taxa.con.sum.out = NULL, lib.size = NULL)
  beta_response <- reactiveValues(cat = NULL, cat.name = NULL)
  cf_response <- reactiveValues(cat = NULL)
  cf_treat <- reactiveValues(cat = NULL)
  response_is_binary <- reactiveValues(binary = NULL)
  rcol = reactiveValues(selected = "lightblue")
  
  ## 0. DATA INPUT -----------
  observeEvent(input$inputOption,{
    observe({
      if (input$inputOption == "Phyloseq") {
        
        shinyjs::hide(id = "optionsInfo")
        output$moreOptions <- renderUI({
          tagList(
            tags$style("
                       .btn-file {
                       border-top-left-radius: 5px !important; border-bottom-left-radius: 5px !important; border-left-style: solid !important; border-left-width: 1px !important;
                       border-top-right-radius: 0px !important; border-bottom-right-radius: 0px !important; border-right-width: 0px !important;
                       }"
            ),
            fileInput("phyloseqData", strong("Please upload your 'phyloseq' data (.rdata, .rds)", style = "color:black"), 
                      accept = c(".rdata", ".rds"), width = '80%'), div(style = "margin-top: -15px"),
            actionButton('Load_Phyloseq_Data', 'Upload', class = "btn-info"), 
            p(" ", style = "margin-bottom: +10px;"), 
            p(strong("Attention:"), "You have to click this Upload button to perform following data processing and downstream data analyses."),br(),
            shinyjs::hidden(
              shiny::div(id = "phyloseqUpload_error",
                         shiny::tags$p("Please upload a .rdata file!!",
                                       style = "color: red; font-weight: bold; padding-top: 5px;", class = "text-center"))
            ),
            INPUT_PHYLOSEQ_COMMENT1,
            p("", style = "margin-bottom:-8px")
          )
        })
        
        output$addDownloadinfo <- renderUI({
          tagList(
            box(title = strong("Example Data", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                downloadButton("downloadData.Immuno", "Immunotherapy", width = '30%', style = "color:black; background-color: red2"),
                br(),br(),
                INPUT_PHYLOSEQ_COMMENT2,
                p("", style = "margin-bottom:-8px")
            )
          )
        })
      } else if (input$inputOption == "Individual Data") {
        shinyjs::hide(id = "optionsInfo")
        output$moreOptions <- renderUI({
          tagList(
            tags$style("
                       .btn-file {
                       border-top-left-radius: 5px !important; border-bottom-left-radius: 5px !important; border-left-style: solid !important; border-left-width: 1px !important;
                       border-top-right-radius: 0px !important; border-bottom-right-radius: 0px !important; border-right-width: 0px !important;
                       }"
            ),
            fileInput("otuTable", strong("Please upload your feature (OTU or ASV) table (.txt, .csv, .biom)", style = "color:black"), 
                      accept = c(".txt", ".csv", ".biom"), width = '80%'), div(style = "margin-top: -15px"),
            fileInput("taxTable", strong("Please upload your taxonomic table (.txt, .tsv)", style = "color:black"), 
                      accept = c(".txt", ".tsv"), width = '80%'), div(style = "margin-top: -15px"),
            fileInput("samData", strong("Please upload your metadata/sample information (.txt, .csv)", style = "color:black"), 
                      accept = c(".txt", ".csv"), width = '80%'), div(style = "margin-top: -15px"),
            fileInput("tree", strong("Please upload your phylogenetic tree (.tre, .nwk)", style = "color:black"), 
                      accept = c(".tre", ".nwk"), width = '80%'), div(style = "margin-top: -15px"),
            actionButton('Load_Individual_Data', 'Upload', class = "btn-info"), br(),br(),
            shinyjs::hidden(
              shiny::div(id = "textfilesUpload_error",
                         shiny::tags$p("Please upload txt and tre files!!",
                                       style = "color: red; font-weight: bold; padding-top: 5px;", class = "text-center"))
            ),
            INPUT_INDIVIDUAL_DATA_COMMENT, 
            p("", style = "margin-bottom:-8px")
          )
        })
        
        output$addDownloadinfo <- renderUI({
          tagList(
            box(title = strong("Example Data", style = "color:white"), 
                width = NULL, status = "info", solidHeader = TRUE,
                downloadButton("downloadZip.Immuno", "Immunotherapy", 
                               width = '30%', style = "color:black; background-color: red2"),
                br(),br(),
                INPUT_INDIVIDUAL_DATA_COMMENT2,
                p("", style = "margin-bottom:-8px")
            )
          )
        })
      }
    })
  }, ignoreInit = TRUE, once = TRUE, ignoreNULL = TRUE)
  
  observe({
    toggleState("Load_Phyloseq_Data", !is.null(input$phyloseqData))
    toggleState("Load_Individual_Data", 
                !(is.null(input$otuTable) | is.null(input$taxTable) | is.null(input$samData)))
    toggleState("batch.yn", !is.null(infile$biom))
    toggleState("run", !is.null(infile$biom))
    toggleState("skip", !is.null(infile$biom))
    toggleState("slider1", !is.null(infile$biom))
    toggleState("slider2", !is.null(infile$biom))
    toggleState("kingdom", !is.null(infile$biom))
    toggleState("binwidth", !is.null(infile$biom))
    toggleState("binwidth2", !is.null(infile$biom))
    
    toggleState("dtRun", !is.null(infile$rare_biom))
    toggleState("datTransRun", !is.null(infile$rare_biom))
  })
  
  observeEvent(input$Load_Phyloseq_Data, {
    
    if (!is.null(input$phyloseqData)) {
      dataInfile  = reactive({
        phyloseq.data = input$phyloseqData
        
        ext <- tools::file_ext(phyloseq.data$datapath)
        
        req(phyloseq.data)
        if (ext == "Rdata") {
          
          phyloseq.dataPath = phyloseq.data$datapath
          e = new.env()
          name <- load(phyloseq.dataPath, envir = e)
          data <- e[[name]]
          
          otu.tab <- otu_table(data, taxa_are_rows = TRUE)
          tax.tab <- tax_table(data)
          sam.dat <- sample_data(data)
          tree <- rtree(nrow(otu.tab))
          tree$tip.label <- rownames(otu.tab)
          data <- merge_phyloseq(otu.tab, tax.tab, sam.dat, tree)
          
          if (sum(sapply(sample_data(data),is.factor))!=0) {
            sample_data(data)[,which(sapply(sample_data(data), is.factor))] = lapply(sample_data(data)[,which(sapply(sample_data(data), is.factor))], as.character)
          }
          
          colnames(tax_table(data)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
          
          if (sum(colnames(otu_table(data)) %in% rownames(sample_data(data))) < sum(rownames(otu_table(data)) %in% rownames(sample_data(data)))) {
            otu_table(data) = t(otu_table(data))
          }
          
          return(data)
        } else if (ext == "rds") {
          phyloseq.dataPath = phyloseq.data$datapath
          data <- readRDS(phyloseq.dataPath)
          
          otu.tab <- otu_table(data, taxa_are_rows = TRUE)
          tree <- rtree(nrow(otu.tab))
          tree$tip.label <- rownames(otu.tab)
          data <- merge_phyloseq(data, tree)
          
          if (sum(sapply(sample_data(data),is.factor))!=0) {
            sample_data(data)[,which(sapply(sample_data(data), is.factor))] = lapply(sample_data(data)[,which(sapply(sample_data(data), is.factor))], as.character)
          }
          colnames(tax_table(data)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
          
          if (sum(colnames(otu_table(data)) %in% rownames(sample_data(data))) < sum(rownames(otu_table(data)) %in% rownames(sample_data(data)))) {
            otu_table(data) = t(otu_table(data))
          }
          
          return(data)
        } else {
          shinyjs::toggle(id = "phyloseqUpload_error", anim = TRUE, time = 1, animType = "fade")
          shinyjs::delay(5000, shinyjs::toggle(id = "phyloseqUpload_error", anim = TRUE, time = 1, animType = "fade"))
          return(NULL)
        }
      })
    } else {
      return(NULL)
    }
    
    if (is.null(dataInfile)) {
      infile$biom <- NULL
      infile$qc_biom <- NULL
      infile$rare_biom = NULL
    } else {
      infile$biom <- dataInfile()
      infile$qc_biom <- dataInfile()
      infile$rare_biom = NULL
    }
    
    updateTabsetPanel(session, "side_menu",
                      selected = "step2")
    rcol$selected = "lightblue"
    
    if (!is.null(infile$biom)) QC$resume()
  })
  observeEvent(input$Load_Individual_Data, {
    shinyjs::disable("Load_Individual_Data")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        incProgress(3/10, message = "File Check")
        if (!is.null(input$otuTable) & !is.null(input$taxTable) & !is.null(input$samData)) {
          dataInfile  = reactive({
            otu.table = input$otuTable
            ext1 <- tools::file_ext(otu.table$datapath)
            
            tax.table = input$taxTable
            ext2 <- tools::file_ext(tax.table$datapath)
            
            sam.data = input$samData
            ext3 <- tools::file_ext(sam.data$datapath)
            
            req(otu.table, tax.table, sam.data)
            if ((ext1 == "txt"| ext1 == "csv" | ext1 == "biom") & (ext2 == "txt" | ext2 == "tsv") &
                (ext3 == "txt" | ext3 == "csv")) {
              otu.table.path = otu.table$datapath
              tax.table.path = tax.table$datapath
              sam.data.path = sam.data$datapath
              
              if (ext1 == "txt") {
                otu.tab <- read.table(otu.table.path, header=TRUE, check.names = FALSE, sep = "\t")
              } else if (ext1 == "csv") {
                otu.tab <- read.csv(otu.table.path, check.names = FALSE)
                rownames(otu.tab) = otu.tab[,1];otu.tab = otu.tab[,-1]
              } else if (ext1 == "biom") {
                biom <- read_biom(otu.table.path)
                otu.tab <- as.matrix(biom_data(biom))
              }
              
              if (ext2 == "txt") {
                tax.tab <- read.table(tax.table.path, header=TRUE, check.names = FALSE, sep = "\t")
              } else if (ext2 == "tsv") {
                tax.tab <- read.table(tax.table.path, header=TRUE, sep="\t")
                tax.tab = preprocess.tax.tab(tax.tab)
              }
              
              if (ext3 == "txt") {
                sam.dat <- read.table(sam.data.path, header=TRUE, check.names = FALSE, sep = "\t")
              } else if (ext3 == "csv") {
                sam.dat <- read.csv(sam.data.path, check.names = FALSE)
                rownames(sam.dat) = sam.dat[,1]
                sam.dat = sam.dat[,-1]
              }
              
              otu.tab <- otu_table(otu.tab, taxa_are_rows = TRUE)
              tax.tab <- tax_table(as.matrix(tax.tab))
              sam.dat <- sample_data(sam.dat)
              
              tree <- rtree(nrow(otu.tab))
              tree$tip.label <- rownames(otu.tab)
              
              if (sum(colnames(otu.tab) %in% rownames(sam.dat)) < sum(rownames(otu.tab) %in% rownames(sam.dat))) {
                otu.tab = t(otu.tab)
              }
              
              incProgress(3/10, message = "Validating")
              shiny::validate(
                if (biom.check.samples(otu.tab, sam.dat)) {
                  if (biom.check.otu(otu.tab, tax.tab)) {
                    showNotification(h4("Error: There is no common samples among OTU/feature table and Sample Data. And
                                        there is no common OTUs among OTU/feature table and taxonomic table."),
                                     type = "error")
                  } else {
                    showNotification(h4("Error: There is no common samples among OTU/feature table and Sample Data"),
                                     type = "error")
                  }
                } else if (biom.check.otu(otu.tab, tax.tab)) {
                  showNotification(h4("Error: There is no common OTUs among OTU/feature table and taxonomic table."),
                                   type = "error")
                } else {
                  NULL
                }
              )
              
              incProgress(1/10, message = "Merging")
              biomData <- merge_phyloseq(otu.tab, tax.tab, sam.dat, tree)
              return(biomData)
            } else {
              shinyjs::toggle(id = "textfilesUpload_error", anim = TRUE, time = 1, animType = "fade")
              shinyjs::delay(5000, shinyjs::toggle(id = "textfilesUpload_error", anim = TRUE, time = 1, animType = "fade"))
              return(NULL)
            }
          })
        } else {
          return(NULL)
        }
        
        if (is.null(dataInfile)) {
          infile$biom <- NULL
          infile$qc_biom <- NULL
          infile$rare_biom = NULL
        } else {
          infile$biom <- dataInfile()
          infile$qc_biom <- dataInfile()
          infile$rare_biom = NULL
        }
        
        updateTabsetPanel(session, "side_menu",
                          selected = "step2")
        rcol$selected = "lightblue"
        
        if (!is.null(infile$biom)) QC$resume()
      })
    shinyjs::enable("Load_Individual_Data")
  })
  
  ## 1-1. QC -----------
  # This reactive expression stores the input data from either the individual data or phyloseq data
  QC = observe(suspended = T,{
    taxa.results$lib.size <- lib.size.func(infile$biom)$lib.size
    
    # Plots graphs using example infile data
    output$hist <- renderPlotly({
      lib_size = lib.size.func(infile$qc_biom)$lib.size
      plot_ly(x = ~lib_size, nbinsx = input$binwidth,
              type = "histogram",
              marker = list(color = rcol$selected, line = list(color = "black", width = 2))) %>%
        layout(
          yaxis = list(title = "Frequency", zeroline = FALSE),
          xaxis = list(title = "Library Size", zeroline = FALSE))
    })
    
    output$hist2 <- renderPlotly({
      mean_prop = mean.prop.func(infile$qc_biom)$mean.prop
      plot_ly(x = ~mean_prop, nbinsx = input$binwidth2,
              type = "histogram",
              marker = list(color = rcol$selected, line = list(color = "black", width = 2))) %>%
        layout(
          yaxis = list(title = "Frequency", zeroline = FALSE),
          xaxis = list(title = "Mean Proportion", zeroline = FALSE))
    })
    
    output$boxplot<- renderPlotly({
      lib_size = lib.size.func(infile$qc_biom)$lib.size
      
      plot_ly(x = ~lib_size, type = "box", notched=TRUE, name = "Library Size",
              color = ~"lib_size", colors = rcol$selected, line = list(color = 'black'))%>%
        layout(
          yaxis = list(title = "", zeroline = FALSE),
          xaxis = list(title = "", zeroline = FALSE), showlegend = FALSE)
    })
    
    output$boxplot2<- renderPlotly({
      mean_prop = mean.prop.func(infile$qc_biom)$mean.prop
      
      plot_ly(x = ~mean_prop, type = "box", notched=TRUE, name = "Mean Proportion",
              color = ~"mean_prop", colors = rcol$selected, line = list(color = 'black'))%>%
        layout(
          yaxis = list(title = "", zeroline = FALSE),
          xaxis = list(title = "", zeroline = FALSE), showlegend = FALSE)
    })
    
    ## Number of Taxonomic Rank for biom before QC
    num_tax.rank = reactive({
      tax.tab = tax_table(infile$qc_biom)
      num.tax.rank(tax.tab)
    })
    
    ## Fills value boxes using example biom data
    output$sample_Size <- renderValueBox({
      valueBox(
        value = tags$p(paste0(lib.size.func(infile$qc_biom)$num.sams), style = "font-size: 75%;"),
        "Sample Size", icon = icon("user-circle"), color = "fuchsia")
    })
    
    output$OTUs_Size <- renderValueBox({
      valueBox(
        value = tags$p(paste0(lib.size.func(infile$qc_biom)$num.otus), style = "font-size: 75%;"),
        "Number of Features", icon = icon("dna"), color = "aqua")
    })
    
    output$phyla <- renderValueBox({
      num.phyla = num_tax.rank()[1]
      valueBox(
        value = tags$p(paste0(num.phyla), style = "font-size: 75%;"),
        "Number of Phyla", icon = icon("sitemap"), color = "orange")
    })
    
    output$classes <- renderValueBox({
      num.classes = num_tax.rank()[2]
      valueBox(
        value = tags$p(paste0(num.classes), style = "font-size: 75%;"),
        "Number of Classes", icon = icon("sitemap"), color = "purple")
    })
    
    output$orders <- renderValueBox({
      num.orders = num_tax.rank()[3]
      valueBox(
        value = tags$p(paste0(num.orders), style = "font-size: 75%;"),
        "Number of Orders", icon = icon("sitemap"), color = "blue")
    })
    
    output$families <- renderValueBox({
      num.families = num_tax.rank()[4]
      valueBox(
        value = tags$p(paste0(num.families), style = "font-size: 75%;"),
        "Number of Families", icon = icon("sitemap"), color = "red")
    })
    
    output$genera <- renderValueBox({
      num.genera = num_tax.rank()[5]
      valueBox(
        value = tags$p(paste0(num.genera), style = "font-size: 75%;"),
        "Number of Genera", icon = icon("sitemap"), color = "lime")
    })
    
    output$species <- renderValueBox({
      num.species = num_tax.rank()[6]
      valueBox(
        value = tags$p(paste0(num.species), style = "font-size: 75%;"),
        "Number of Species", icon = icon("sitemap"), color = "teal" )
    })
    
    ## This event handler checks whether there is an input file and updates the slider options accordingly
    maxi.slider1 = as.numeric(lib.size.func(infile$qc_biom)$lib.size.sum["3rd quartile"])
    updateSliderInput(session, "slider1", min = 0, max = round(maxi.slider1,-3))
  })
  
  observeEvent(input$batch.yn, {
    if(input$batch.yn == "Yes"){
      shinyjs::show("batch")
      shinyjs::show("prim")
      shinyjs::show("covar")
      shinyjs::show("batch.method")
      output$batch.method.select <- renderUI({
        tagList(
          prettyRadioButtons("batch.method", 
                             label = h4(strong("Method", style = "color:black")), 
                             animation = "jelly",
                             choices = c("ConQuR"), 
                             selected = "ConQuR", 
                             icon = icon("check"), 
                             width = '80%')
        )
      })
      
      output$batch <- renderUI({
        tagList(
          p(" ", style = "margin-top: 25px;"),
          h4(strong("Batch variable", style = "color:black")),
          p("A variable for batch IDs (e.g., labs, studies, states, locations, times)", style = "font-size:10pt"),
          p(" ", style = "margin-bottom: +15px;"),
          selectInput("batch.var", label = NULL,
                      choices = c(bmc.col.check(sample_data(infile$biom), type = "Multinomial"), bmc.col.check(sample_data(infile$biom), type = "Binary")),
                      selected = c(bmc.col.check(sample_data(infile$biom), type = "Multinomial"), bmc.col.check(sample_data(infile$biom), type = "Binary"))[1],
                      width = '70%')
        )
      })
    } else {
      shinyjs::hide("batch.method")
      shinyjs::hide("batch")
      shinyjs::hide("prim")
      shinyjs::hide("covar")
    }
  })
  
  observeEvent(input$batch.var, {
    observeEvent(input$batch.yn, {
      output$prim <- renderUI({
        tagList(
          p(" ", style = "margin-top: 25px;"),
          h4(strong("Primary variable", style = "color:black")),
          p("A binary or continuous variable on the patient's recovery (required). ConQuR maintains the variability in the microbiome data due to this primary variable.", style = "font-size:10pt"),
          p(" ", style = "margin-bottom: +15px;"),
          selectInput("prim.var", label = NULL,
                      choices = get.cov.col(sample_data(infile$biom))[!get.cov.col(sample_data(infile$biom)) == input$batch.var],
                      selected = get.cov.col(sample_data(infile$biom))[!get.cov.col(sample_data(infile$biom)) == input$batch.var][1],
                      width = '70%')
        )
      })
    })
  })
  
  observeEvent(input$batch.var, {
    output$covar <- renderUI({
      tagList(
        p(" ", style = "margin-top: 25px;"),
        h4(strong("Nuisance variable(s)", style = "color:black")),
        p("Other nuisance variable(s) on the patients' characteristics (optional). ConQuR maintains the variability in the microbiome data due to these nuisance variable(s).", style = "font-size:10pt"),
        p(" ", style = "margin-bottom: +15px;"),
        prettyCheckboxGroup("covar", label = NULL,
                            choices = get.cov.col(sample_data(infile$biom))[!get.cov.col(sample_data(infile$biom)) %in% c(input$batch.var, input$prim.var)],
                            width = '70%')
      )
    })
  })
  
  ## 2. Descriptive -------------------
  
  #0804 
  observeEvent(chooseData$sam.dat, {
    
    output$de_response <- renderUI({
      tagList(
        h4(strong("Response variable", style = "color:black")),
        p("A binary or continuous variable on the patient's health or disease status after treatment.", style = "font-size:10pt"), 
        selectInput("response_sel_de", label = NULL,
                    c("Choose one" = "", sort(names(chooseData$sam.dat))), selected = "Response", width = '70%'))
    })
    
    observeEvent(input$response_sel_de, {
      
      if (length(table(chooseData$sam.dat[[input$response_sel_de]])) == 2){
        
        shinyjs::show("de_response_rename")
        
        de.categos <<- sort(category.names(chooseData$sam.dat, input$response_sel_de))
        
        output$de_response_rename <- renderUI({
          tagList(
            h5(strong("Rename the response variable", style = "color:black")),
            p("Rename categories for response variable", style = "font-size:10pt"),
            textInput("de_cat_no_res", label = (paste0("Non-response (0): ", de.categos[1])), value = de.categos[1], width = '80%'),
            textInput("de_cat_res", label = (paste0("Response (1): ", de.categos[2])), value = de.categos[2], width = '80%')
          )
        })
        
        output$de_chooseTest <- renderUI({
          selectInput("de_method", label = h4(strong("Method", style = "color:black")), 
                      choices = c("Fisher's exact test (Default)", "Pearson's Chi-squared test"), width = '98%')
        })
        
        output$de_legend <- renderUI({
          tagList(
            prettyRadioButtons("de_legend_sel", label = h4(strong("Legend location", style = "color:black")),  choices = c("Top left", "Top right", "Top left on the left panel only", "Top right on the left panel only", "Top left on the right panel only", "Top right on the right panel only"), 
                               icon = icon("check"), animation = "jelly", selected = "Top left", width = '80%'),
            p(" ", style = "margin-bottom: -10px;"),
            p("You can change the location of the figure legend to minimize overlaps with the other parts in the graph."),
            p(" ", style = "margin-bottom: +15px;"),
            actionButton("run_des", (strong("Run!")), class = "btn-info")
          )
        })
      } else {
        shinyjs::hide("de_response_rename")
        
        output$de_chooseTest <- renderUI({
          selectInput("de_method", label = h4(strong("Method", style = "color:black")), 
                      choices = c("Mann-Whitney test (Default)", "Welch's t-test"), width = '98%')})
        
        output$de_legend <- renderUI({
          tagList(
            prettyRadioButtons("de_legend_sel", label = h4(strong("Legend location", style = "color:black")),  choices = c("Top left", "Top right"), shape = c("round"), selected = "Top left", width = '80%'),
            p(" ", style = "margin-bottom: -10px;"),
            p("You can change the location of the figure legend to minimize overlaps with the other parts in the graph."),
            p(" ", style = "margin-bottom: +15px;"),
            actionButton("run_des", (strong("Run!")), class = "btn-info")
          )
        }) 
      }
      
      treat_vars_options <<- sel.var.binary(chooseData$sam.dat, input$response_sel_de)
      
      output$de_treat <- renderUI({
        tagList(
          p(" ", style = "margin-bottom: +15px;"),
          h4(strong("Treatment variable", style = "color:black")),
          p("A binary variable on the treatment status (e.g., placebo vs. drug, old drug vs. new drug).", style = "font-size:10pt"), 
          selectInput("treat_sel_de", label = NULL,
                      choices = sort(treat_vars_options), selected = "Treatment", width = '70%'))
        
      })
      
      observeEvent(input$treat_sel_de, {
        
        if (length(table(chooseData$sam.dat[[input$treat_sel_de]])) == 2){
          
          shinyjs::show("de_treat_rename")
          
          treat.categos <<- sort(category.names(chooseData$sam.dat, input$treat_sel_de))
          
          output$de_treat_rename <- renderUI({
            
            tagList(
              h5(strong("Rename the treatment variable", style = "color:black")),
              p("Rename categories for the treatment variable", style = "font-size:10pt"),
              textInput("treat_con", label = (paste0("Control (0): ", treat.categos[1])), value = treat.categos[1], width = '80%'),
              textInput("treat_tr", label = (paste0("Treatment (1): ", treat.categos[2])), value = treat.categos[2], width = '80%') 
            )
          })
          
        } else {
          shinyjs::hide("de_treat_rename")
        }
      }) 
    }) 
  })
  
  ## 3. Generalized Linear Models -------------------
  observeEvent(chooseData$taxa.out, {
    observeEvent(chooseData$sam.dat, {
      
      output$int_format <- renderUI({
        prettyRadioButtons("int_dataType", label = h4(strong("Data format", style = "color:black")), icon = icon("check"), animation = "jelly",
                           c("CLR (Default)", "Count (Rarefied)", "Proportion", "Arcsine-root"), selected = "CLR (Default)", width = '70%')
      })
      
      output$int_response <- renderUI({
        tagList(
          h4(strong("Response variable", style = "color:black")),
          p("A binary or continuous variable on the post-treatment status of the patients' health or disease.", style = "font-size:10pt"), 
          selectInput("response_sel_int", label = NULL,
                      sort(names(chooseData$sam.dat)), selected = "Response", width = '70%'))
      })
      
      observeEvent(input$response_sel_int, {
        
        if (length(table(chooseData$sam.dat[[input$response_sel_int]])) == 2){
          
          shinyjs::show("int_response_rename")
          
          int.categos <<- sort(category.names(chooseData$sam.dat, input$response_sel_int))
          
          output$int_response_rename <- renderUI({
            tagList(
              h5(strong("Rename the response variable", style = "color:black")),
              p("Rename categories for response variable", style = "font-size:10pt"),
              textInput("int_cat_no_res", label = (paste0("Non-response (0): ", int.categos[1])), value = int.categos[1], width = '80%'),
              textInput("int_cat_res", label = (paste0("Response (1): ", int.categos[2])), value = int.categos[2], width = '80%')
            )
          })
          
          output$int_chooseTest <- renderUI({
            selectInput("int_method", label = h4(strong("Method", style = "color:black")), 
                        choices = c("Logistic regression"), width = '98%')
          })
          
          output$int_legend <- renderUI({
            tagList(
              prettyRadioButtons("int_legend_sel", label = h4(strong("Legend location", style = "color:black")),  choices = c("Top left", "Top right", "Bottom left", "Bottom right"), icon = icon("check"), animation = "jelly", selected = "Top left", width = '80%'),
              p(" ", style = "margin-bottom: -10px;"),
              p("You can change the location of the figure legend to minimize overlaps with the other parts in the graph."),
              p(" ", style = "margin-bottom: +20px;"),
              h4(strong("Taxonomic ranks", style = "color:black")),
              p(" ", style = "margin-bottom: -10px;"),
              prettyRadioButtons("tax_rank", label = NULL,  animation = "jelly",
                                 c("Phylum - Genus (Default)", "Phylum - Species"), selected = "Phylum - Genus (Default)",
                                 icon = icon("check"), width = '80%'),
              p(" ", style = "margin-bottom: -10px;"),
              p("The taxonomic ranks to be surveyed. 'Phylum to Genus (Default)' for 16S rRNA amplicon sequencing; 'Phylum to Species' for shotgun metagenomic sequencing.", style = "font-size:10pt"),
              
              actionButton("run_int", (strong("Run!")), class = "btn-info")
              
            )
          })
        } else {
          
          shinyjs::hide("int_response_rename")
          
          output$int_chooseTest <- renderUI({
            selectInput("int_method", label = h4(strong("Method", style = "color:black")), 
                        choices = c("Linear regression"), width = '98%')})
          
          output$int_legend <- renderUI({
            tagList(
              prettyRadioButtons("int_legend_sel", label = h4(strong("Legend location", style = "color:black")),  choices = c("Top left", "Top right", "Bottom left", "Bottom right"), shape = c("round"), selected = "Top left", width = '80%'),
              p(" ", style = "margin-bottom: -10px;"),
              p("You can change the location of the figure legend to minimize overlaps with the other parts in the graph."),
              p(" ", style = "margin-bottom: +20px;"),
              h4(strong("Taxonomic ranks", style = "color:black")),
              p(" ", style = "margin-bottom: -10px;"),
              prettyRadioButtons("tax_rank", label = NULL,  animation = "jelly",
                                 c("Phylum - Genus (Default)", "Phylum - Species"), selected = "Phylum - Genus (Default)",
                                 icon = icon("check"), width = '80%'),
              p(" ", style = "margin-bottom: -10px;"),
              p("The taxonomic ranks to be surveyed. 'Phylum to Genus (Default)' for 16S rRNA amplicon sequencing; 'Phylum to Species' for shotgun metagenomic sequencing.", style = "font-size:10pt"),
              actionButton("run_int", (strong("Run!")), class = "btn-info")
            )
          })
        }
        
        treat_vars_options <<- sel.var.binary(chooseData$sam.dat, input$response_sel_int)
        
        output$int_treat <- renderUI({
          tagList(
            h4(strong("Treatment variable", style = "color:black")),
            p("A binary variable on the treatment status (e.g., placebo vs. drug, old drug vs. new drug).", style = "font-size:10pt"), 
            selectInput("treat_sel_int", label = NULL,
                        choices = sort(treat_vars_options), selected = "Treatment", width = '70%'))
        })
        
        observeEvent(input$treat_sel_int, {
          
          if (length(table(chooseData$sam.dat[[input$treat_sel_int]])) == 2){
            shinyjs::show("int_treat_rename")
            treat.categos <<- sort(category.names(chooseData$sam.dat, input$treat_sel_int))
            
            output$int_treat_rename <- renderUI({
              tagList(
                h5(strong("Rename the treatment variable", style = "color:black")),
                p("Rename categories for the treatment variable", style = "font-size:10pt"),
                textInput("treat_con_int", label = (paste0("Control (0): ", treat.categos[1])), value = treat.categos[1], width = '80%'),
                textInput("treat_tr_int", label = (paste0("Treatment (1): ", treat.categos[2])), value = treat.categos[2], width = '80%') 
              )
            })
          } else {
            shinyjs::hide("int_treat_rename")
          }
          
          int_cov_options <- select.covariates.func(chooseData$sam.dat, input$treat_sel_int, input$response_sel_int)
          
          output$int_cov <- renderUI({
            tagList(
              p(" ", style = "margin-top: 25px;"),
              h4(strong("Covariate(s)", style = "color:black")),
              p("Potential confounders (e.g., age, gender) to be adjusted for.", style = "font-size:10pt"),
              p(" ", style = "margin-bottom: +15px;"),
              prettyRadioButtons("int_covariate",label = NULL, 
                                 icon = icon("check"),
                                 animation = "jelly", c("None", "Covariate(s)"), selected = "None", width = '70%'),
              p(" ", style = "margin-bottom: -10px;"),
              
              shinyjs::hidden(
                shiny::div(id = "int_covariateOptions", style = "margin-left: 2%",
                           prettyCheckboxGroup("int_cov_sel"," Please select covariate(s)", 
                                               choices = int_cov_options , width = '70%'))))
          })
          
          observeEvent(input$int_covariate,{
            if (input$int_covariate == "Covariate(s)") {
              shinyjs::show("int_covariateOptions")
            } 
            else if (input$int_covariate == "None") {
              shinyjs::hide("int_covariateOptions")
            }
          })
          output$int_num_taxa <- renderUI({
          })
        }) 
      }) 
    })
  })
  
  ## 4. Causal Forest -------------------
  output$cf_format <- renderUI({
    tagList(
      prettyRadioButtons("cf_dataType", label = h4(strong("Data format", style = "color:black")), icon = icon("check"), animation = "jelly",
                         c("CLR (Default)", "Count (Rarefied)", "Proportion", "Arcsine-root"), selected = "CLR (Default)",width = '70%'))
  })
  
  output$cf_response <- renderUI({
    tagList(
      p(" ", style = "margin-top: 25px;"),
      h4(strong("Response variable", style = "color:black")),
      p("A binary or continuous variable on the post-treatment status of the patients’ health or disease.", style = "font-size:10pt"),
      p(" ", style = "margin-bottom: +15px;"),
      selectInput("cf_response", label = NULL, 
                  choices = c(bmc.col.check(data.frame(chooseData$sam.dat), type = "Binary"), bmc.col.check(data.frame(chooseData$sam.dat), type = "Continuous")), 
                  selected = c(bmc.col.check(data.frame(chooseData$sam.dat), type = "Binary"), bmc.col.check(data.frame(chooseData$sam.dat), type = "Continuous"))[1], width = '70%')
    )
  })
  
  observeEvent(input$cf_response, {
    if(length(table(chooseData$sam.dat[[input$cf_response]])) == 2){
      cf_response$cat = category.names(chooseData$sam.dat, input$cf_response)
      cf_response_length = length(cf_response$cat)
      response_is_binary$binary = TRUE
      
      output$cf_response_rename <- renderUI({
        tagList(
          h5(strong("Rename categories for response variable", style = "color:black")),
          p("You can rename the categories of response variable. MiCML keeps up to 8 characters in the output graphs.", style = "font-size:10pt"),
          lapply(1:cf_response_length, function(i){
            textInput(paste0("cf_response_label", i), label = (paste0("Group / Category ", i, " : ",cf_response$cat[i])), value = cf_response$cat[i], width = '80%')
          })
        )
      })
    } else {
      response_is_binary$binary = FALSE
      
      output$cf_response_rename <- renderUI({
        tagList(
          h5(strong("Rename the response variable", style = "color:black")),
          p("You can rename the response variable. MiCML keeps up to 8 characters in the output graph.", style = "font-size:10pt"),
          textInput(paste0("cf_response_label"), label = NULL, value = input$cf_response, width = '80%')
          
        )
      })
    }
  })
  
  output$cf_treat <- renderUI({
    tagList(
      p(" ", style = "margin-top: 25px;"),
      h4(strong("Treatment variable", style = "color:black")),
      p("A binary variable on the primary treatment status (e.g,. placebo vs. treatment, old treatment vs. new treatment).", style = "font-size:10pt"),
      p(" ", style = "margin-bottom: +15px;"),
      selectInput("cf_treat", label = NULL, 
                  choices = setdiff(bmc.col.check(chooseData$sam.dat, type = "Binary"), input$cf_response), 
                  selected = setdiff(bmc.col.check(chooseData$sam.dat, type = "Binary"), input$cf_response)[1], width = '70%')
    )
  })
  
  observeEvent(input$cf_treat, {
    if(length(table(chooseData$sam.dat[[input$cf_treat]])) == 2){
      cf_treat$cat = category.names(chooseData$sam.dat, input$cf_treat)
      cf_treat_length = length(cf_treat$cat)
      
      output$cf_treat_rename <- renderUI({
        tagList(
          h5(strong("Rename categories for treatment variable", style = "color:black")),
          p("You can rename the categories of treatment variable. MiCML keeps up to 8 characters in the output graphs.", style = "font-size:10pt"),
          lapply(1:cf_treat_length, function(i){
            textInput(paste0("cf_treat_label", i), label = (paste0("Group / Category ", i, " : ",cf_treat$cat[i])), value = cf_treat$cat[i], width = '80%')
          })
        )
      })
    } else {
      output$cf_treat_rename <- renderUI({
        tagList(
          h5(strong("Rename the treatment variable", style = "color:black")),
          p("You can rename the treatment variable. MiCML keeps up to 8 characters in the output graph.", style = "font-size:10pt"),
          textInput(paste0("cf_treat_label"), label = NULL, value = input$cf_treat, width = '80%')
        )
      })
    }
  })
  
  output$cf_method <- renderUI({
    tagList(
      p(" ", style = "margin-top: 25px;"),
      h4(strong("Method", style = "color:black")),
      p("You can choose 'Double-sample tree' for a randomized controlled trial or 'Propensity tree' for an observational study.", style = "font-size:10pt"),
      p(" ", style = "margin-bottom: +15px;"),
      prettyRadioButtons(inputId = "cf_method",
                         label = NULL,
                         icon = icon("check"), animation = "jelly",
                         choices = c("Double-sample tree", "Propensity tree"), selected = "Propensity tree", width = '70%'),)
  })
  
  observeEvent(input$cf_method, {
    output$cf_cov <- renderUI({
      tagList(
        p(" ", style = "margin-top: 25px;"),
        h4(strong("Covariate(s)", style = "color:black")),
        p("Potential confounders (e.g., age, gender) to be adjusted for.", style = "font-size:10pt"),
        p(" ", style = "margin-bottom: +15px;"),
        prettyRadioButtons("cf_cov",label = NULL, icon = icon("check"),
                           animation = "jelly", c("None", "Covariate(s)"), selected = "None", width = '70%'),
        shinyjs::hidden(
          shiny::div(id = "cf_cov_list", style = "margin-left: 2%",
                     prettyCheckboxGroup("cf_cov_options"," Please select covariate(s)", choices = setdiff(colnames(chooseData$sam.dat), c(input$cf_response, input$cf_treat)), width = '70%')))
      )
    })
    
    observeEvent(input$cf_method, {
      if(input$cf_method == "Double-sample tree") shinyjs::hide("cf_cov")
      else shinyjs::show("cf_cov")
    })
    
    observeEvent(input$cf_cov, {
      if(input$cf_cov == "None") shinyjs::hide("cf_cov_list")
      else shinyjs::show("cf_cov_list")
    })
    
    output$cf_tax_rank <- renderUI({
      tagList(
        p(" ", style = "margin-top: 25px;"),
        h4(strong("Taxonomic ranks", style = "color:black")),
        p("The taxonomic ranks to be surveyed. 'Phylum to Genus (Default)' for 16S rRNA amplicon sequencing; 'Phylum to Species' for shotgun metagenomic sequencing.", style = "font-size:10pt"),
        p(" ", style = "margin-bottom: +15px;"),
        prettyRadioButtons("cf_tax_rank", label = NULL, icon = icon("check"), animation = "jelly",
                           c("Phylum - Genus (Default)", "Phylum - Species"), selected = "Phylum - Genus (Default)", width = '70%')
      )
    })
    
    output$cf_num_taxa <- renderUI({
      tagList(
        p(" ", style = "margin-top: 25px;"),
        h4(strong("# taxa to be displayed", style = "color:black")),
        p("The maximum number of taxa to be displayed in variable importance and partial dependence plots (Default: 20).", style = "font-size:10pt"),
        p(" ", style = "margin-bottom: +15px;"),
        sliderInput("cf_num_taxa", label = NULL, min = 5, max = 20, value = 20, step = 5),
        
        actionButton("cf_runButton", strong("Run!"), class = "btn-info"),
      )
    })
  })
  
  ## Run Buttons ------
  
  ## 1-1. QC -----
  observeEvent(input$run, {
    if(input$batch.yn == "Yes"){
      withProgress(
        message = 'Calculation in progress', 
        detail = 'This may take a while...', value = 0, {
          
          incProgress(1/10, message = "Data Trimming in progress")
          
          if (nchar(input$part.rem.str) == 0) {
            rem.tax.complete <- rem.tax.d
            rem.tax.partial <- rem.tax.str.d
          } else {
            rem.tax.complete <- unique(c(unlist(strsplit(input$rem.str, split = ",")), rem.tax.d))
            rem.tax.partial <- unique(c(unlist(strsplit(input$part.rem.str, split = ",")), rem.tax.str.d))
          }
          
          tax.tab <- tax_table(infile$biom)
          
          if (input$kingdom != "all") {
            ind <- is.element(tax.tab[,1], input$kingdom)
            shiny::validate(
              if (sum(ind) == 0) {
                showNotification(h4(paste("Error: Please select valid Kingdom. Available kingdoms are:", 
                                          paste(c(na.omit(unique(tax.tab[,1])) ,"and all"), collapse = ", "))), type = "error")
              } else {
                NULL
              }
            )
          }
          
          shinyjs::disable("batch.yn")
          shinyjs::disable("batch.method.select")
          shinyjs::disable("batch.var")
          shinyjs::disable("prim.var")
          shinyjs::disable("covar")
          shinyjs::disable("run")
          shinyjs::disable("slider1")
          shinyjs::disable("slider2")
          shinyjs::disable("kingdom")
          shinyjs::disable("skip")
          shinyjs::disable("binwidth")
          shinyjs::disable("binwidth2")
          shinyjs::disable("rem.str")
          shinyjs::disable("part.rem.str")
          
          rcol$selected = "rgba(255, 0, 0, 0.6)"
          
          tree.exists <- !is.null(access(infile$biom, "phy_tree"))
          
          # 1. QC
          
          infile$qc_biom = biom.clean(infile$biom, 
                                      input$kingdom, 
                                      lib.size.cut.off = input$slider1, 
                                      mean.prop.cut.off = input$slider2/100,
                                      rem.tax = rem.tax.complete, rem.tax.str = rem.tax.partial,
                                      tree.exists = tree.exists)
          
          incProgress(2/10, message = "Batch effect correction in progress")
          
          # 2. Remove NA
          
          otu.tab <- as.data.frame(t(otu_table(infile$qc_biom)))
          sam.dat <- sample_data(infile$qc_biom)
          
          remain.samples <- sam.dat[,c(input$batch.var, input$prim.var, input$covar)] %>% na.omit %>% rownames
          
          infile$na_omit_biom <- prune_samples(remain.samples, infile$qc_biom)
          if(input$batch.method == "ConQuR"){
            otu.tab <- as.data.frame(t(otu_table(infile$na_omit_biom)))
          } else {
            otu.tab <- as.matrix(otu_table(infile$na_omit_biom))
          }
          tax.tab <- tax_table(infile$na_omit_biom)
          tree <- phy_tree(infile$na_omit_biom)
          sam.dat <- sample_data(infile$na_omit_biom)
          
          # 3. Check numeric / binary factor / multicategory factor
          
          batchid <<- as.factor(sam.dat[[input$batch.var]])
          
          prim.var <- input$prim.var
          if(!is.null(input$covar)){
            cov <- input$covar
            df <- data.frame(sam.dat[,prim.var], sam.dat[,cov])
          } else {
            df <- data.frame(sam.dat[,prim.var])
          }
          
          for(name in colnames(df)){
            type <- col.str.check(sam.dat, name)
            if(type == "factor"){
              df[[name]] <- as.factor(df[[name]])
            }
            else if(type == "numeric"){
              df[[name]] <- as.numeric(df[[name]])
            } else {
              df[[name]] <- df[[name]]
            }
          }
          
          f1 <- as.formula(paste("~" ,paste(names(df), collapse = "+"), collapse = ""))
          covar <- model.matrix(f1, data = df)[,-1]
          
          if(input$batch.method == "ConQuR"){
            ref.bat <- names(which.max(table(batchid)))
            set.seed(578)
            try(adjusted.otu.tab <- ConQuR(tax_tab = otu.tab, batchid = batchid,
                                           covariates = covar, batch_ref = ref.bat,
                                           logistic_lasso = T, quantile_type = "lasso", interplt = T, num_core = 1), silent = TRUE)
            bat.otu.tab <- otu_table(t(as.data.frame(as.matrix(adjusted.otu.tab))), taxa_are_rows = TRUE)
          }
          else{
            set.seed(578)
            adjusted.otu.tab <- sva::ComBat_seq(otu.tab, batch=batchid, group=NULL, covar_mod = covar)
            bat.otu.tab <- otu_table(as.data.frame(as.matrix(adjusted.otu.tab)), taxa_are_rows = TRUE)
          }
          
          bat.sam.dat <- sample_data(infile$na_omit_biom)
          bat.tax.tab <- tax_table(infile$na_omit_biom)
          bat.tree <- phy_tree(infile$na_omit_biom)
          
          batch.infile$bat_biom <- merge_phyloseq(bat.otu.tab, bat.tax.tab, bat.sam.dat, bat.tree)
          
          # 4. Rarefying
          incProgress(2/10, message = "Rarefying in progress")
          
          # Rarefying original biom data
          lib_size.sum = lib.size.func(infile$qc_biom)$lib.size.sum
          infile$rare_biom = rarefy.func(infile$qc_biom, 
                                         cut.off = lib_size.sum["Minimum"],
                                         multi.rarefy = 1,
                                         tree.exists = tree.exists)
          
          # Rarefying corrected biom data
          batch.infile$qc_biom <- otu.tab.clean(batch.infile$bat_biom,
                                                lib.size.cut.off = input$slider1,
                                                mean.prop.cut.off = input$slider2/100,
                                                tree.exists = tree.exists)
          
          lib_size.sum = lib.size.func(batch.infile$qc_biom)$lib.size.sum
          batch.infile$rare_biom <- rarefy.func(batch.infile$qc_biom,
                                                cut.off = lib_size.sum["Minimum"],
                                                multi.rarefy = 1,
                                                tree.exists = tree.exists)
          
          # 5. PCoA Plot code
          qc.biom <- infile$qc_biom
          rare.biom <- infile$rare_biom
          qc.biom.bat <- batch.infile$qc_biom
          rare.biom.bat <- batch.infile$rare_biom
          
          shinyjs::show("pcoa.area")
          
          incProgress(2/10, message = "Plotting PCoA in progress")
          output$batch.pcoa <- renderPlot({
            try(batch.correct.PCoA(qc.biom = qc.biom, rare.biom = rare.biom,
                                   qc.biom.bat = qc.biom.bat, rare.biom.bat = rare.biom.bat,
                                   batch.var = input$batch.var), silent = TRUE)
          })
          
          incProgress(2/10, message = "Saving File in progress")
          
          chooseData$sam.dat = sample_data(batch.infile$qc_biom)
          chooseData$mon.sin.rev.bin.con = is.mon.sin.rev.bin.con(chooseData$sam.dat)
          chooseData$prim_vars = pri.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con)
          chooseData$tax.tab = tax_table(batch.infile$rare_biom)
          
          output$moreControls <- renderUI({
            tagList(
              box(title = strong("Download Data", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                  span(textOutput("text"), style="font-size:12pt"),
                  
                  h5("Data after Quality Control and Batch Effect Correction"),
                  downloadButton("downloadData2", "Download", width = '50%', 
                                 style = "color:black; background-color: red3"),br(),
                  
                  h5("Data after Quality Control, Batch Effect Correction and Rarefaction"),
                  downloadButton("downloadData3", "Download", width = '50%', 
                                 style = "color:black; background-color: red3"),br(),
                  
                  p("For your reference, you can download the data files above for 
                  (1) the phyloseq object (biom.after.qc.bat) after QC and batch effect correction, and 
                  (2) the phyloseq object (rare.biom.after.qc.bat) after QC, batch effect correction, and rarefaction.", 
                    style = "font-size:11pt")))
          })
          
          output$text <- renderText({"You are all set! You can proceed to data analysis!"})
          
          observeEvent(input$batch.yn, {
            if(input$batch.yn == "Yes"){
              output$qc_reference <- renderUI({
                tagList(
                  box(title = strong("Reference", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE, QC_BATCH_REFERENCE)
                )
              })
            }
          })
          
          biom.after.qc.bat = batch.infile$qc_biom
          output$downloadData2 <- downloadHandler(
            filename = function() {
              paste("biom.after.qc.bat.Rdata")
            },
            content = function(file1) {
              save(biom.after.qc.bat, file = file1)
            })
          
          rare.biom.after.qc.bat = batch.infile$rare_biom
          output$downloadData3 <- downloadHandler(
            filename = function() {
              paste("rare.biom.after.qc.bat.Rdata")
            },
            content = function(file1) {
              save(rare.biom.after.qc.bat, file = file1)
            })
          
          incProgress(1/10, message = "Done")
          shinyjs::enable("batch.yn")
          shinyjs::enable("batch.method.select")
          shinyjs::enable("batch.var")
          shinyjs::enable("prim.var")
          shinyjs::enable("covar")
          shinyjs::enable("run")
          shinyjs::enable("slider1")
          shinyjs::enable("slider2")
          shinyjs::enable("kingdom")
          shinyjs::enable("skip")
          shinyjs::enable("binwidth")
          shinyjs::enable("binwidth2")
          shinyjs::enable("rem.str")
          shinyjs::enable("part.rem.str")
          
          infile$qc_biom <- batch.infile$qc_biom
          infile$rare_biom <- batch.infile$rare_biom
        })
    } else {
      withProgress(
        message = 'Calculation in progress', 
        detail = 'This may take a while...', value = 0, {
          
          incProgress(1/10, message = "Data Trimming in progress")
          
          if (nchar(input$part.rem.str) == 0) {
            rem.tax.complete <- rem.tax.d
            rem.tax.partial <- rem.tax.str.d
          } else {
            rem.tax.complete <- unique(c(unlist(strsplit(input$rem.str, split = ",")), rem.tax.d))
            rem.tax.partial <- unique(c(unlist(strsplit(input$part.rem.str, split = ",")), rem.tax.str.d))
          }
          
          tax.tab <- tax_table(infile$biom)
          
          if (input$kingdom != "all") {
            ind <- is.element(tax.tab[,1], input$kingdom)
            shiny::validate(
              if (sum(ind) == 0) {
                showNotification(h4(paste("Error: Please select valid Kingdom. Available kingdoms are:", 
                                          paste(c(na.omit(unique(tax.tab[,1])) ,"and all"), collapse = ", "))), type = "error")
              } else {
                NULL
              }
            )
          }
          
          shinyjs::disable("batch.yn")
          shinyjs::disable("batch.method.select")
          shinyjs::disable("batch.var")
          shinyjs::disable("prim.var")
          shinyjs::disable("covar")
          shinyjs::disable("run")
          shinyjs::disable("slider1")
          shinyjs::disable("slider2")
          shinyjs::disable("kingdom")
          shinyjs::disable("skip")
          shinyjs::disable("binwidth")
          shinyjs::disable("binwidth2")
          shinyjs::disable("rem.str")
          shinyjs::disable("part.rem.str")
          
          rcol$selected = "rgba(255, 0, 0, 0.6)"
          
          tree.exists <- !is.null(access(infile$biom, "phy_tree"))
          
          infile$qc_biom = biom.clean(infile$biom, 
                                      input$kingdom, 
                                      lib.size.cut.off = input$slider1, 
                                      mean.prop.cut.off = input$slider2/100,
                                      rem.tax = rem.tax.complete, rem.tax.str = rem.tax.partial,
                                      tree.exists = tree.exists)
          
          incProgress(3/10, message = "Rarefying in progress")
          lib_size.sum = lib.size.func(infile$qc_biom)$lib.size.sum
          infile$rare_biom = rarefy.func(infile$qc_biom, 
                                         cut.off = lib_size.sum["Minimum"],
                                         multi.rarefy = 1,
                                         tree.exists = tree.exists)
          
          shinyjs::hide("pcoa.area")
          
          incProgress(2/10, message = "Saving File in progress")
          
          chooseData$sam.dat = sample_data(infile$qc_biom)
          chooseData$mon.sin.rev.bin.con = is.mon.sin.rev.bin.con(chooseData$sam.dat)
          chooseData$prim_vars = pri.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con)
          chooseData$tax.tab = tax_table(infile$rare_biom)
          
          output$moreControls <- renderUI({
            tagList(
              box(title = strong("Download Data", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                  span(textOutput("text"), style="font-size:12pt"),
                  h5("Data after Quality Control"),
                  downloadButton("downloadData2", "Download", width = '50%', style = "color:black; background-color: red3"),br(),
                  h5("Data after Quality Control and Rarefaction"),
                  downloadButton("downloadData3", "Download", width = '50%', style = "color:black; background-color: red3"),br(),
                  p("For your reference, you can download the data files above for the phyloseq object (biom.after.qc) after QC and
                      (rare.biom.after.qc) after QC and rarefaction.",
                    style = "font-size:11pt")
              )
            )
          })
          
          output$text <- renderText({"You are all set! You can proceed to data analysis!"})
          
          biom.after.qc = infile$qc_biom
          output$downloadData2 <- downloadHandler(
            filename = function() {
              paste("biom.after.qc.Rdata")
            },
            content = function(file1) {
              save(biom.after.qc, file = file1)
            })
          
          rare.biom.after.qc = infile$rare_biom
          output$downloadData3 <- downloadHandler(
            filename = function() {
              paste("rare.biom.after.qc.Rdata")
            },
            content = function(file1) {
              save(rare.biom.after.qc, file = file1)
            })
          
          incProgress(1/10, message = "Done")
          shinyjs::enable("batch.yn")
          shinyjs::disable("batch.method.select")
          shinyjs::enable("batch.var")
          shinyjs::enable("prim.var")
          shinyjs::enable("covar")
          shinyjs::enable("run")
          shinyjs::enable("slider1")
          shinyjs::enable("slider2")
          shinyjs::enable("kingdom")
          shinyjs::enable("skip")
          shinyjs::enable("binwidth")
          shinyjs::enable("binwidth2")
          shinyjs::enable("rem.str")
          shinyjs::enable("part.rem.str")
        })
    }
  })
  
  ## 1-2. Data Transformation -----
  observeEvent(input$dtRun, {
    withProgress(
      message = 'Calculation in progress', 
      detail = 'This may take a while...', value = 0, {
        shinyjs::disable("dtRun")
        
        incProgress(2/10, message = "Data transformation in progress")
        rare.otu.tab <- otu_table(infile$rare_biom)
        rare.tax.tab <- tax_table(infile$rare_biom)
        no.rare.otu.tab <- otu_table(infile$qc_biom)
        no.rare.tax.tab <- tax_table(infile$qc_biom)
        
        chooseData$taxa.out = tax.trans(no.rare.otu.tab, no.rare.tax.tab, rare.otu.tab, rare.tax.tab)
        chooseData$taxa.names.out = taxa.names.rank(chooseData$taxa.out[[1]])
        chooseData$tax.tab = rare.tax.tab
        
        incProgress(2/10, message = "Saving")
        
        output$dtDownload <- renderUI({
          tagList(
            box(title = strong("Download Data", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                span(textOutput("text"), style="font-size:15pt"),
                # strong(p("Taxonomic Abundance Data",style = "font-size:11pt")),
                h5("Count", HTML('&emsp;'), HTML('&emsp;'), HTML('&emsp;'), HTML('&emsp;'), HTML('&emsp;'), "Count (Rarefied)"),
                downloadButton("taxadataCount", "Download", width = '50%', style = " background-color: red3"), HTML('&emsp;'),
                downloadButton("taxadataRareCount", "Download", width = '50%', style = " background-color: red3"), br(), 
                h5("Proportion", HTML('&emsp;'), HTML('&emsp;'), HTML('&emsp;'), HTML('&nbsp;'), "CLR"),
                downloadButton("taxadataProp", "Download", width = '50%', style = " background-color: red3"), HTML('&emsp;'),
                downloadButton("taxadataCLR", "Download", width = '50%', style = " background-color: red3"), br(),
                h5("Arcsine-root"),
                downloadButton("taxadataArc", "Download", width = '50%', style = " background-color: red3"), br(), p("", style = "margin-bottom:5px")
            )
          )
        })
        
        count_biom = chooseData$taxa.out$count
        rare_biom = chooseData$taxa.out$rare.count
        prop_biom = chooseData$taxa.out$prop
        clr_biom = chooseData$taxa.out$clr
        arc_biom = chooseData$taxa.out$arcsin
        
        output$taxadataCount <- downloadHandler(
          filename = function() {
            paste("Count.Data.zip")},
          content = function(count.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(count_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=count.file, files=dataFiles)
          }
        )
        
        output$taxadataRareCount <- downloadHandler(
          filename = function() {
            paste("Rarefied.Count.Data.zip")},
          content = function(rare.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(rare_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=rare.file, files=dataFiles)
          }
        )
        output$taxadataProp <- downloadHandler(
          filename = function() {
            paste("Proportion.Data.zip")},
          content = function(prop.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(prop_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=prop.file, files=dataFiles)
          }
        )
        
        output$taxadataCLR <- downloadHandler(
          filename = function() {
            paste("CLR.Transformed.Data.zip")},
          content = function(clr.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(clr_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=clr.file, files=dataFiles)
          }
        )
        
        output$taxadataArc <- downloadHandler(
          filename = function() {
            paste("Arcsin.Transformed.Data.zip")},
          content = function(arc.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(arc_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=arc.file, files=dataFiles)
          }
        )
        incProgress(1/10, message = "Done")
        shinyjs::enable("dtRun")
      })
    
    shinyjs::show("cf_format")
    shinyjs::show("cf_response")
    shinyjs::show("cf_response_rename")
    shinyjs::show("cf_treat")
    shinyjs::show("cf_treat_rename")
    shinyjs::show('cf_method')
    shinyjs::show("cf_cov")
    shinyjs::show("cf_tax_rank")
    shinyjs::show("cf_num_taxa")
    shinyjs::show("cf_download")
    shinyjs::show("cf_reference")
  })
  
  ## 2. Descriptive -------------------
  
  observeEvent(input$run_des, {
    
    shinyjs::disable("run_des")
    shinyjs::disable("de_response")
    shinyjs::disable("de_response_rename")
    shinyjs::disable("de_treat")
    shinyjs::disable("de_treat_rename")
    shinyjs::disable("de_chooseTest")
    shinyjs::disable("de_legend")
    
    withProgress(
      message = "Calculation in progress",
      detail = "This may take a while...", value = 0, {
        
        treat.categos <- sort(category.names(chooseData$sam.dat, input$treat_sel_de))
        new.treat.categos <- c(input$treat_con, input$treat_tr)
        
        Treatment <- try(rename.bin.var (chooseData$sam.dat[[input$treat_sel_de]], treat.categos, new.treat.categos), silent = TRUE) 
        
        treat.na.ind <- which(is.na(Treatment))
        
        legend.loc <- input$de_legend_sel
        
        if (length(table(chooseData$sam.dat[[input$response_sel_de]])) == 2){
          
          res.categos <- sort(category.names(chooseData$sam.dat, input$response_sel_de))
          new.res.categos <- c(input$de_cat_no_res, input$de_cat_res)
          
          Response <- try(rename.bin.var (chooseData$sam.dat[[input$response_sel_de]], res.categos, new.res.categos), silent = TRUE) 
          
          response.na.ind <- which(is.na(Response))
          na.ind <- unique(c(treat.na.ind, response.na.ind))
          
          if (length(na.ind) != 0){
            Treatment <- Treatment[-na.ind]
            Response <- Response[-na.ind]
          }
          
          if(input$de_method == "Fisher's exact test (Default)"){
            
            fisher.result <- tryCatch(bin.fisher.test(Response, Treatment), error = function(e) {
              message ("No outcome is available!")
              showModal(modalDialog(div("No outcome is available!")))
              return(NULL)
            })
            
            output$de_display = renderUI({
              tagList(
                box(title = strong("Fisher's exact test", style = "color:white"),  
                    align = "center", width = NULL, solidHeader = TRUE,  status = "info", 
                    plotOutput("fisher_plot", height = 800, width = 650)
                )
              )
            })
            
            output$de_reference <- renderUI({
              tagList(
                box(title = strong("Reference", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                    p(FISHER_REFERENCE, style = "font-size:11pt"))
              )
            })
            
            output$fisher_plot = try(renderPlot({bin.des.plot(Response, Treatment, legend.loc, fisher.result, "Fisher's exact test (Default)")}), silent = TRUE)
            
          } else if (input$de_method == "Pearson's Chi-squared test"){
            
            chisq.result <- tryCatch(bin.chisq.test(Response, Treatment), error = function(e) {
              message ("No outcome is available!")
              showModal(modalDialog(div("No outcome is available!")))
              return(NULL)
            })
            
            output$de_reference <- renderUI({
              tagList(
                box(title = strong("Reference", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                    p(PEARSON_REFERENCE, style = "font-size:11pt"))
              )
            })
            
            output$de_display = renderUI({
              tagList(
                box(title = strong("Pearson's Chi-squared test", style = "color:white"),  
                    align = "center", width = NULL, solidHeader = TRUE,  status = "info", 
                    plotOutput("chisq_plot", height = 800, width = 650)
                )
              )
            })
            
            output$chisq_plot = try(renderPlot({bin.des.plot(Response, Treatment, legend.loc, chisq.result, "Pearson's Chi-squared test")}), silent = TRUE)
            
          }
          
        } else {
          Response <- as.numeric (chooseData$sam.dat[[input$response_sel_de]]) 
          response.na.ind <- which(is.na(Response))
          na.ind <- unique(c(treat.na.ind, response.na.ind))
          
          if (length(na.ind) != 0){
            Treatment <- Treatment[-na.ind]
            Response <- Response[-na.ind]
          }
          
          if (input$de_method == "Mann-Whitney test (Default)"){
            
            wilcox.result <- tryCatch(con.wilcox.test(Response, Treatment), error = function(e) {
              message ("No outcome is available!")
              showModal(modalDialog(div("No outcome is available!")))
              return(NULL)
            })
            output$de_reference <- renderUI({
              tagList(
                box(title = strong("Reference", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                    p(MANN_WHITNEY_REFERENCE, style = "font-size:11pt"))
              )
            })
            
            output$de_display = renderUI({
              tagList(
                box(title = strong("Mann-Whitney test", style = "color:white"),  
                    align = "center", width = NULL, solidHeader = TRUE,  status = "info", 
                    plotOutput("wilcox_plot", height = 800, width = 650)
                )
              )
            })
            output$wilcox_plot <- try(renderPlot({con.des.plot(Response, Treatment, new.treat.categos, legend.loc, wilcox.result, input$de_method)}), silent = TRUE) 
            
          } else if (input$de_method == "Welch's t-test"){
            
            welch.result <- tryCatch(con.welch.test(Response, Treatment), error = function(e) {
              message ("No outcome is available!")
              showModal(modalDialog(div("No outcome is available!")))
              return(NULL)
            })
            
            output$de_reference <- renderUI({
              tagList(
                box(title = strong("Reference", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                    p(WELCH_REFERENCE, style = "font-size:11pt"))
              )
            })
            
            output$de_display = renderUI({
              tagList(
                box(title = strong("Welch's t-test", style = "color:white"),  
                    align = "center", width = NULL, solidHeader = TRUE,  status = "info", 
                    plotOutput("welch_plot", height = 800, width = 650)
                )
              )
            })
            output$welch_plot <- try(renderPlot({con.des.plot(Response, Treatment, new.treat.categos, legend.loc, welch.result, input$de_method)}), silent = TRUE) 
          }
        }
      }) 
    
    shinyjs::enable("run_des")
    shinyjs::enable("de_response")
    shinyjs::enable("de_response_rename")
    shinyjs::enable("de_treat")
    shinyjs::enable("de_treat_rename")
    shinyjs::enable("de_chooseTest")
    shinyjs::enable("de_legend")
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  ## 3. Interaction Effect -------------------
  
  
  observeEvent(input$run_int, {
    
    shinyjs::disable("run_int")
    shinyjs::disable("int_format")
    shinyjs::disable("int_response")
    shinyjs::disable("int_response_rename")
    shinyjs::disable("int_treat")
    shinyjs::disable("int_treat_rename")
    shinyjs::disable("int_cov")
    shinyjs::disable("int_chooseTestv")
    shinyjs::disable("int_num_taxa")
    shinyjs::disable("int_legend")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        
        if (input$int_dataType == "Count (Rarefied)") {
          taxa.dataType = "rare.count"
        } else if (input$int_dataType == "CLR (Default)") {
          taxa.dataType = "clr"
        } else if (input$int_dataType == "Arcsine-root") {
          taxa.dataType = "arcsin"  
        } else if (input$int_dataType == "Proportion"){
          taxa.dataType = "imp.prop"
        }
        
        if (input$tax_rank == "Phylum - Genus (Default)"){
          include = FALSE 
        } else {
          include = TRUE
        }
        
        taxa.out <- chooseData$taxa.out[[taxa.dataType]]
        sam.dat <- chooseData$sam.dat
        
        treat.categos <- sort(category.names(chooseData$sam.dat, input$treat_sel_int))
        new.treat.categos <- c(input$treat_con_int, input$treat_tr_int)
        
        Treatment <- try(rename.bin.var (chooseData$sam.dat[[input$treat_sel_int]], treat.categos, new.treat.categos), silent = TRUE) 
        
        legend.loc.int <- input$int_legend_sel
        num.taxa <- input$num_taxa_sel
        
        if (length(table(chooseData$sam.dat[[input$response_sel_int]])) == 2){
          
          res.categos <- sort(category.names(chooseData$sam.dat, input$response_sel_int))
          new.res.categos <- c(input$int_cat_no_res, input$int_cat_res)
          
          Response <- try(rename.bin.var (chooseData$sam.dat[[input$response_sel_int]], res.categos, new.res.categos), silent = TRUE) 
          
          treat.na.ind <- which(is.na(Treatment))
          response.na.ind <- which(is.na(Response))
          
          na.ind <- unique(c(treat.na.ind, response.na.ind))
          
          if (length(na.ind) != 0){
            Treatment <- Treatment[-na.ind]
            Response <- Response[-na.ind]
            
            for (i in 1:(5+include)){
              taxa.out[[i]] <- taxa.out[[i]][-na.ind, ]
            }
          }
          
          if(input$int_method == "Logistic regression"){ 
            
            if (input$int_covariate == "None"){
              
              Treatment.Levels <- substr(levels(Treatment), 1, 8)
              Treatment <- as.numeric(Treatment) - 1
              
              reg_result <- tryCatch(logistic.regression.no.cov(taxa.out, include, Treatment, Treatment.Levels, Response), error = function(e) {
                message ("No outcome is available!")
                showModal(modalDialog(div("No outcome is available!")))
                return(NULL)
              })
              
              height_plot <- try(height_output(reg_result, num.taxa, include), silent = TRUE) 
              
              if (include) {
                
                incProgress(3/10, message = "Displaying Results in progress")
                
                output$int_display = renderUI({
                  tagList(
                    tabBox(title = strong("Logistic Regression", style = "color:black", side = "right"), width = NULL,
                           tabPanel("Phylum", align = "center",
                                    plotOutput("rank1", height = height_plot[[1]]),  
                           )
                           ,
                           tabPanel("Class", align = "center",
                                    plotOutput("rank2", height = height_plot[[2]]),
                           )
                           ,tabPanel("Order", align = "center",
                                     plotOutput("rank3", height = height_plot[[3]]),
                           )
                           ,tabPanel("Family", align = "center",
                                     plotOutput("rank4", height = height_plot[[4]]),
                           )
                           ,tabPanel("Genus", align = "center",
                                     plotOutput("rank5", height = height_plot[[5]]),
                           )
                           ,tabPanel("Species", align = "center",
                                     plotOutput("rank6", height = height_plot[[6]]),
                           )
                    )
                  )
                })
                
                output$rank1 = try(renderPlot({ 
                  output.to.show (taxa.out, reg_result, Treatment, Treatment.Levels, Response, 1, num.taxa, legend.loc.int)
                }), silent = TRUE)
                
                incProgress(4/10, message = "Displaying Results in progress: Phylum")
                
                output$rank2 = try(renderPlot({ 
                  output.to.show (taxa.out, reg_result, Treatment, Treatment.Levels, Response, 2, num.taxa, legend.loc.int)
                }), silent = TRUE)
                
                incProgress(5/10, message = "Displaying Results in progress: Class")
                
                output$rank3 = try(renderPlot({ 
                  output.to.show (taxa.out, reg_result, Treatment, Treatment.Levels, Response, 3, num.taxa, legend.loc.int)
                }), silent = TRUE)
                
                incProgress(6/10, message = "Displaying Results in progress: Order")
                
                output$rank4 = try(renderPlot({ 
                  output.to.show (taxa.out, reg_result, Treatment, Treatment.Levels, Response, 4, num.taxa, legend.loc.int)
                }), silent = TRUE)
                
                incProgress(7/10, message = "Displaying Results in progress: Family")
                
                output$rank5 = try(renderPlot({ 
                  output.to.show (taxa.out, reg_result, Treatment, Treatment.Levels, Response, 5, num.taxa, legend.loc.int)
                }), silent = TRUE)
                
                incProgress(8/10, message = "Displaying Results in progress: Genus")
                
                output$rank6 = try(renderPlot({ 
                  output.to.show (taxa.out, reg_result, Treatment, Treatment.Levels, Response, 6, num.taxa, legend.loc.int)
                }), silent = TRUE)
                
                incProgress(9/10, message = "Displaying Results in progress: Species")
                
              } else {
                incProgress(3/10, message = "Displaying Results in progress")
                
                output$int_display = renderUI({
                  tagList(
                    tabBox(title = strong("Logistic Regression", style = "color:black", side = "right"), width = NULL,
                           tabPanel("Phylum", align = "center",
                                    plotOutput("rank1", height = height_plot[[1]]),
                           )
                           ,
                           tabPanel("Class", align = "center",
                                    plotOutput("rank2", height = height_plot[[2]]),
                           )
                           ,tabPanel("Order", align = "center",
                                     plotOutput("rank3", height = height_plot[[3]]),
                           )
                           ,tabPanel("Family", align = "center",
                                     plotOutput("rank4", height = height_plot[[4]]),
                           )
                           ,tabPanel("Genus", align = "center",
                                     plotOutput("rank5", height = height_plot[[5]]),
                           )
                    )
                  )
                })
                
                output$rank1 = try(renderPlot({ 
                  output.to.show (taxa.out, reg_result, Treatment, Treatment.Levels, Response, 1, num.taxa, legend.loc.int)
                }), silent = TRUE)
                
                incProgress(4/10, message = "Displaying Results in progress: Phylum")
                
                output$rank2 = try(renderPlot({ 
                  output.to.show (taxa.out, reg_result, Treatment, Treatment.Levels, Response, 2, num.taxa, legend.loc.int)
                }), silent = TRUE)
                
                incProgress(5/10, message = "Displaying Results in progress: Class")
                
                output$rank3 = try(renderPlot({ 
                  output.to.show (taxa.out, reg_result, Treatment, Treatment.Levels, Response, 3, num.taxa, legend.loc.int)
                }), silent = TRUE)
                
                incProgress(6/10, message = "Displaying Results in progress: Order")
                
                output$rank4 = try(renderPlot({ 
                  output.to.show (taxa.out, reg_result, Treatment, Treatment.Levels, Response, 4, num.taxa, legend.loc.int)
                }), silent = TRUE)
                
                incProgress(7/10, message = "Displaying Results in progress: Family")
                
                output$rank5 = try(renderPlot({ 
                  output.to.show (taxa.out, reg_result, Treatment, Treatment.Levels, Response, 5, num.taxa, legend.loc.int)
                }), silent = TRUE)
                
                incProgress(8/10, message = "Displaying Results in progress: Genus")
                
              }
            } else {
              
              Treatment.Levels <- substr(levels(Treatment), 1, 8)
              Treatment <- as.numeric(Treatment) - 1
              
              Covariates <- try(cov.factorize.func(chooseData$sam.dat, input$int_cov_sel), silent = TRUE) 
              
              treat.na.ind <- which(is.na(Treatment))
              response.na.ind <- which(is.na(Response))
              cov.na.ind <- which(rowSums(is.na(Covariates)) > 0)
              
              na.ind <- unique(c(treat.na.ind, response.na.ind, cov.na.ind))
              
              if (length(na.ind) != 0){
                Treatment <- Treatment[-na.ind]
                Response <- Response[-na.ind]
                Covariates <- Covariates[-na.ind, ]
                
                for (i in 1:(5+include)){
                  taxa.out[[i]] <- taxa.out[[i]][-na.ind, ]
                }
              }
              
              reg_result <<- tryCatch(logistic.regression.with.cov(taxa.out, include, Treatment, Treatment.Levels, Covariates, Response), error = function(e) {
                message ("No outcome is available!")
                showModal(modalDialog(div("No outcome is available!")))
                return(NULL)
              })
              
              height_plot <- try(height_output(reg_result, num.taxa, include), silent = TRUE) 
              
              if (include) {
                incProgress(3/10, message = "Displaying Results in progress")
                
                output$int_display = renderUI({
                  tagList(
                    tabBox(title = strong("Logistic Regression", style = "color:black", side = "right"), width = NULL,
                           tabPanel("Phylum", align = "center",
                                    plotOutput("rank1", height = height_plot[[1]]),  
                           )
                           ,
                           tabPanel("Class", align = "center",
                                    plotOutput("rank2", height = height_plot[[2]]),
                           )
                           ,tabPanel("Order", align = "center",
                                     plotOutput("rank3", height = height_plot[[3]]),
                           )
                           ,tabPanel("Family", align = "center",
                                     plotOutput("rank4", height = height_plot[[4]]),
                           )
                           ,tabPanel("Genus", align = "center",
                                     plotOutput("rank5", height = height_plot[[5]]),
                           )
                           ,tabPanel("Species", align = "center",
                                     plotOutput("rank6", height = height_plot[[6]]),
                           )
                    )
                  )
                })
                
                output_result_phylum <<- result_for_output(taxa.out, reg_result, Treatment, Treatment.Levels, Response, 1, num.taxa, covariates = Covariates, response.type = "binary")
                output_result_class <<- result_for_output(taxa.out, reg_result, Treatment, Treatment.Levels, Response, 2, num.taxa, covariates = Covariates, response.type = "binary")
                output_result_order <<- result_for_output(taxa.out, reg_result, Treatment, Treatment.Levels, Response, 3, num.taxa, covariates = Covariates, response.type = "binary")
                output_result_family <<- result_for_output(taxa.out, reg_result, Treatment, Treatment.Levels, Response, 4, num.taxa, covariates = Covariates, response.type = "binary")
                output_result_genus <<- result_for_output(taxa.out, reg_result, Treatment, Treatment.Levels, Response, 5, num.taxa, covariates = Covariates, response.type = "binary")
                output_result_species <<- result_for_output(taxa.out, reg_result, Treatment, Treatment.Levels, Response, 6, num.taxa, covariates = Covariates, response.type = "binary")
                
                output$rank1 = try(renderPlot({ 
                  output.to.show.2 (taxa.out, reg_result, output_result_phylum, Treatment, Treatment.Levels, Response, 1, num.taxa, legend.loc.int, Covariates, "binary")
                }), silent = TRUE)
                
                incProgress(4/10, message = "Displaying Results in progress: Phylum")
                
                
                output$rank2 = try(renderPlot({ 
                  output.to.show.2 (taxa.out, reg_result, output_result_class, Treatment, Treatment.Levels, Response, 2, num.taxa, legend.loc.int, Covariates, "binary")
                }), silent = TRUE)
                
                incProgress(5/10, message = "Displaying Results in progress: Class")
                
                output$rank3 = try(renderPlot({ 
                  output.to.show.2 (taxa.out, reg_result, output_result_order, Treatment, Treatment.Levels, Response, 3, num.taxa, legend.loc.int, Covariates, "binary")
                }), silent = TRUE)
                
                incProgress(6/10, message = "Displaying Results in progress: Order")
                
                output$rank4 = try(renderPlot({ 
                  output.to.show.2 (taxa.out, reg_result, output_result_family, Treatment, Treatment.Levels, Response, 4, num.taxa, legend.loc.int, Covariates, "binary")
                }), silent = TRUE)
                
                incProgress(7/10, message = "Displaying Results in progress: Family")
                
                output$rank5 = try(renderPlot({ 
                  output.to.show.2 (taxa.out, reg_result, output_result_genus, Treatment, Treatment.Levels, Response, 5, num.taxa, legend.loc.int, Covariates, "binary")
                }), silent = TRUE)
                
                
                output$rank6 = try(renderPlot({ 
                  output.to.show.2 (taxa.out, reg_result, output_result_species, Treatment, Treatment.Levels, Response, 6, num.taxa, legend.loc.int, Covariates, "binary")
                }), silent = TRUE)
                
                
              } else {
                
                incProgress(3/10, message = "Displaying Results in progress")
                
                output$int_display = renderUI({
                  tagList(
                    tabBox(title = strong("Logistic Regression", style = "color:black", side = "right"), width = NULL,
                           tabPanel("Phylum", align = "center",
                                    plotOutput("rank1", height = height_plot[[1]]),
                           )
                           ,
                           tabPanel("Class", align = "center",
                                    plotOutput("rank2", height = height_plot[[2]]),
                           )
                           ,tabPanel("Order", align = "center",
                                     plotOutput("rank3", height = height_plot[[3]]),
                           )
                           ,tabPanel("Family", align = "center",
                                     plotOutput("rank4", height = height_plot[[4]]),
                           )
                           ,tabPanel("Genus", align = "center",
                                     plotOutput("rank5", height = height_plot[[5]]),
                           )
                    )
                  )
                })
                
                output_result_phylum <<- result_for_output(taxa.out, reg_result, Treatment, Treatment.Levels, Response, 1, num.taxa, covariates = Covariates, response.type = "binary")
                output_result_class <<- result_for_output(taxa.out, reg_result, Treatment, Treatment.Levels, Response, 2, num.taxa, covariates = Covariates, response.type = "binary")
                output_result_order <<- result_for_output(taxa.out, reg_result, Treatment, Treatment.Levels, Response, 3, num.taxa, covariates = Covariates, response.type = "binary")
                output_result_family <<- result_for_output(taxa.out, reg_result, Treatment, Treatment.Levels, Response, 4, num.taxa, covariates = Covariates, response.type = "binary")
                output_result_genus <<- result_for_output(taxa.out, reg_result, Treatment, Treatment.Levels, Response, 5, num.taxa, covariates = Covariates, response.type = "binary")
                
                output$rank1 = try(renderPlot({ 
                  output.to.show.2 (taxa.out, reg_result, output_result_phylum, Treatment, Treatment.Levels, Response, 1, num.taxa, legend.loc.int, Covariates, "binary")
                }), silent = TRUE)
                
                output$rank2 = try(renderPlot({ 
                  output.to.show.2 (taxa.out, reg_result, output_result_class, Treatment, Treatment.Levels, Response, 2, num.taxa, legend.loc.int, Covariates, "binary")
                }), silent = TRUE)
                
                output$rank3 = try(renderPlot({ 
                  output.to.show.2 (taxa.out, reg_result, output_result_order, Treatment, Treatment.Levels, Response, 3, num.taxa, legend.loc.int, Covariates, "binary")
                }), silent = TRUE)
                
                output$rank4 = try(renderPlot({ 
                  output.to.show.2 (taxa.out, reg_result, output_result_family, Treatment, Treatment.Levels, Response, 4, num.taxa, legend.loc.int, Covariates, "binary")
                }), silent = TRUE)
                
                
                output$rank5 = try(renderPlot({ 
                  output.to.show.2 (taxa.out, reg_result, output_result_genus, Treatment, Treatment.Levels, Response, 5, num.taxa, legend.loc.int, Covariates, "binary")
                }), silent = TRUE)
                
              }
            }
            
            if(taxa.dataType =="clr"){
              output$int_reference <- renderUI({
                tagList(
                  box(title = strong("Reference", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                      p(GLM_REFERENCE_CLR, style = "font-size:11pt"))
                )
              })
            } else {
              output$int_reference <- renderUI({
                tagList(
                  box(title = strong("Reference", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                      p(GLM_REFERENCE_NO_CLR, style = "font-size:11pt"))
                )
              })
            }
          }
        } else if (length(table(chooseData$sam.dat[[input$response_sel_int]])) > 2) {
          
          Response <- chooseData$sam.dat[[input$response_sel_int]]
          
          if (input$int_covariate == "None"){
            
            Treatment.Levels <- substr(levels(Treatment), 1, 8)
            Treatment <- as.numeric(Treatment) - 1
            
            treat.na.ind <- which(is.na(Treatment))
            response.na.ind <- which(is.na(Response))
            
            na.ind <- unique(c(treat.na.ind, response.na.ind))
            
            if (length(na.ind) != 0){
              Treatment <- Treatment[-na.ind]
              Response <- Response[-na.ind]
              
              for (i in 1:(5+include)){
                taxa.out[[i]] <- taxa.out[[i]][-na.ind, ]
              }
            }
            
            reg_result <- tryCatch(gaussian.regression.no.cov(taxa.out, include, Treatment, Treatment.Levels, Response), error = function(e) {
              message ("No outcome is available!")
              showModal(modalDialog(div("No outcome is available!")))
              return(NULL)
            })
            
            height_plot <- try(height_output(reg_result, num.taxa, include), silent = TRUE) 
            
            if (include) {
              incProgress(2/10, message = "Displaying Results in progress")
              
              output$int_display = renderUI({
                tagList(
                  tabBox(title = NULL, width = NULL,
                         tabPanel("Phylum", align = "center",
                                  plotOutput("rank1", height = height_plot[[1]]),  
                         )
                         ,
                         tabPanel("Class", align = "center",
                                  plotOutput("rank2", height = height_plot[[2]]),
                         )
                         ,tabPanel("Order", align = "center",
                                   plotOutput("rank3", height = height_plot[[3]]),
                         )
                         ,tabPanel("Family", align = "center",
                                   plotOutput("rank4", height = height_plot[[4]]),
                         )
                         ,tabPanel("Genus", align = "center",
                                   plotOutput("rank5", height = height_plot[[5]]),
                         )
                         ,tabPanel("Species", align = "center",
                                   plotOutput("rank6", height = height_plot[[6]]),
                         )
                  )
                )
              })
              
              output$rank1 = try(renderPlot({ 
                output.to.show (taxa.out, reg_result, Treatment, Treatment.Levels, Response, 1, num.taxa, legend.loc.int, covariates = NULL, response.type = "continuous")
              }), silent = TRUE)
              
              incProgress(3/10, message = "Displaying Results in progress: Phylum")
              
              output$rank2 = try(renderPlot({ 
                output.to.show (taxa.out, reg_result, Treatment, Treatment.Levels, Response, 2, num.taxa, legend.loc.int, covariates = NULL, response.type = "continuous")
              }), silent = TRUE)
              
              incProgress(4/10, message = "Displaying Results in progress: Class")
              
              output$rank3 = try(renderPlot({ 
                output.to.show (taxa.out, reg_result, Treatment, Treatment.Levels, Response, 3, num.taxa, legend.loc.int, covariates = NULL, response.type = "continuous")
              }), silent = TRUE)
              
              incProgress(5/10, message = "Displaying Results in progress: Order")
              
              output$rank4 = try(renderPlot({ 
                output.to.show (taxa.out, reg_result, Treatment, Treatment.Levels, Response, 4, num.taxa, legend.loc.int, covariates = NULL, response.type = "continuous")
              }), silent = TRUE)
              
              incProgress(6/10, message = "Displaying Results in progress: Family")
              
              output$rank5 = try(renderPlot({ 
                output.to.show (taxa.out, reg_result, Treatment, Treatment.Levels, Response, 5, num.taxa, legend.loc.int, covariates = NULL, response.type = "continuous")
              }), silent = TRUE)
              
              incProgress(7/10, message = "Displaying Results in progress: Genus")
              
              output$rank6 = try(renderPlot({ 
                output.to.show (taxa.out, reg_result, Treatment, Treatment.Levels, Response, 6, num.taxa, legend.loc.int, covariates = NULL, response.type = "continuous")
              }), silent = TRUE)
              
              incProgress(8/10, message = "Displaying Results in progress: Species")
              
            } else {
              incProgress(3/10, message = "Displaying Results in progress")
              
              output$int_display = renderUI({
                tagList(
                  tabBox(title = NULL, width = NULL,
                         tabPanel("Phylum", align = "center",
                                  plotOutput("rank1", height = height_plot[[1]]),
                         )
                         ,
                         tabPanel("Class", align = "center",
                                  plotOutput("rank2", height = height_plot[[2]]),
                         )
                         ,tabPanel("Order", align = "center",
                                   plotOutput("rank3", height = height_plot[[3]]),
                         )
                         ,tabPanel("Family", align = "center",
                                   plotOutput("rank4", height = height_plot[[4]]),
                         )
                         ,tabPanel("Genus", align = "center",
                                   plotOutput("rank5", height = height_plot[[5]]),
                         )
                  )
                )
              })
              
              output$rank1 = try(renderPlot({ 
                output.to.show (taxa.out, reg_result, Treatment, Treatment.Levels, Response, 1, num.taxa, legend.loc.int, covariates = NULL, response.type = "continuous")
              }), silent = TRUE)
              
              incProgress(3/10, message = "Displaying Results in progress: Phylum")
              
              output$rank2 = try(renderPlot({ 
                output.to.show (taxa.out, reg_result, Treatment, Treatment.Levels, Response, 2, num.taxa, legend.loc.int, covariates = NULL, response.type = "continuous")
              }), silent = TRUE)
              
              incProgress(4/10, message = "Displaying Results in progress: Class")
              
              output$rank3 = try(renderPlot({ 
                output.to.show (taxa.out, reg_result, Treatment, Treatment.Levels, Response, 3, num.taxa, legend.loc.int, covariates = NULL, response.type = "continuous")
              }), silent = TRUE)
              
              incProgress(5/10, message = "Displaying Results in progress: Order")
              
              output$rank4 = try(renderPlot({ 
                output.to.show (taxa.out, reg_result, Treatment, Treatment.Levels, Response, 4, num.taxa, legend.loc.int, covariates = NULL, response.type = "continuous")
              }), silent = TRUE)
              
              incProgress(6/10, message = "Displaying Results in progress: Family")
              
              output$rank5 = try(renderPlot({ 
                output.to.show (taxa.out, reg_result, Treatment, Treatment.Levels, Response, 5, num.taxa, legend.loc.int, covariates = NULL, response.type = "continuous")
              }), silent = TRUE)
              
              incProgress(7/10, message = "Displaying Results in progress: Genus")
              
            }
          } else {
            
            Treatment.Levels <- substr(levels(Treatment), 1, 8)
            Treatment <- as.numeric(Treatment) - 1
            Covariates <- try(cov.factorize.func(chooseData$sam.dat, input$int_cov_sel), silent = TRUE) 
            
            treat.na.ind <- which(is.na(Treatment))
            response.na.ind <- which(is.na(Response))
            cov.na.ind <- which(rowSums(is.na(Covariates)) > 0)
            
            na.ind <- unique(c(treat.na.ind, response.na.ind, cov.na.ind))
            
            if (length(na.ind) != 0){
              Treatment <- Treatment[-na.ind]
              Response <- Response[-na.ind]
              Covariates <- Covariates[-na.ind, ]
              
              for (i in 1:(5+include)){
                taxa.out[[i]] <- taxa.out[[i]][-na.ind, ]
              }
            }
            
            reg_result <- tryCatch(gaussian.regression.with.cov(taxa.out, include, Treatment, Treatment.Levels, Covariates, Response), error = function(e) {
              message ("No outcome is available!")
              showModal(modalDialog(div("No outcome is available!")))
              return(NULL)
            })
            
            height_plot <- height_output(reg_result, num.taxa, include)
            
            if (include) {
              incProgress(3/10, message = "Displaying Results in progress")
              
              output$int_display = renderUI({
                tagList(
                  tabBox(title = NULL, width = NULL,
                         tabPanel("Phylum", align = "center",
                                  plotOutput("rank1", height = height_plot[[1]]),  
                         )
                         ,
                         tabPanel("Class", align = "center",
                                  plotOutput("rank2", height = height_plot[[2]]),
                         )
                         ,tabPanel("Order", align = "center",
                                   plotOutput("rank3", height = height_plot[[3]]),
                         )
                         ,tabPanel("Family", align = "center",
                                   plotOutput("rank4", height = height_plot[[4]]),
                         )
                         ,tabPanel("Genus", align = "center",
                                   plotOutput("rank5", height = height_plot[[5]]),
                         )
                         ,tabPanel("Species", align = "center",
                                   plotOutput("rank6", height = height_plot[[6]]),
                         )
                  )
                )
              })
              
              output_result_phylum <<- result_for_output(taxa.out, reg_result, Treatment, Treatment.Levels, Response, 1, num.taxa, covariates = Covariates, response.type = "continuous")
              output_result_class <<- result_for_output(taxa.out, reg_result, Treatment, Treatment.Levels, Response, 2, num.taxa, covariates = Covariates, response.type = "continuous")
              output_result_order <<- result_for_output(taxa.out, reg_result, Treatment, Treatment.Levels, Response, 3, num.taxa, covariates = Covariates, response.type = "continuous")
              output_result_family <<- result_for_output(taxa.out, reg_result, Treatment, Treatment.Levels, Response, 4, num.taxa, covariates = Covariates, response.type = "continuous")
              output_result_genus <<- result_for_output(taxa.out, reg_result, Treatment, Treatment.Levels, Response, 5, num.taxa, covariates = Covariates, response.type = "continuous")
              output_result_species <<- result_for_output(taxa.out, reg_result, Treatment, Treatment.Levels, Response, 6, num.taxa, covariates = Covariates, response.type = "continuous")
              
              output$rank1 = try(renderPlot({ 
                output.to.show.2 (taxa.out, reg_result, output_result_phylum, Treatment, Treatment.Levels, Response, 1, num.taxa, legend.loc.int, Covariates, "continuous")
              }), silent = TRUE)
              
              output$rank2 = try(renderPlot({ 
                output.to.show.2 (taxa.out, reg_result, output_result_class, Treatment, Treatment.Levels, Response, 2, num.taxa, legend.loc.int, Covariates, "continuous")
              }), silent = TRUE)
              
              output$rank3 = try(renderPlot({ 
                output.to.show.2 (taxa.out, reg_result, output_result_order, Treatment, Treatment.Levels, Response, 3, num.taxa, legend.loc.int, Covariates, "continuous")
              }), silent = TRUE)
              
              output$rank4 = try(renderPlot({ 
                output.to.show.2 (taxa.out, reg_result, output_result_family, Treatment, Treatment.Levels, Response, 4, num.taxa, legend.loc.int, Covariates, "continuous")
              }), silent = TRUE)
              
              
              output$rank5 = try(renderPlot({ 
                output.to.show.2 (taxa.out, reg_result, output_result_genus, Treatment, Treatment.Levels, Response, 5, num.taxa, legend.loc.int, Covariates, "continuous")
              }), silent = TRUE)
              
              output$rank6 = try(renderPlot({ 
                output.to.show.2 (taxa.out, reg_result, output_result_genus, Treatment, Treatment.Levels, Response, 6, num.taxa, legend.loc.int, Covariates, "continuous")
              }), silent = TRUE)
              
              
            } else {
              incProgress(3/10, message = "Displaying Results in progress")
              
              output$int_display = renderUI({
                tagList(
                  tabBox(title = strong("Logistic Regression", style = "color:black", side = "right"), width = NULL,
                         tabPanel("Phylum", align = "center",
                                  plotOutput("rank1", height = height_plot[[1]]),
                         )
                         ,
                         tabPanel("Class", align = "center",
                                  plotOutput("rank2", height = height_plot[[2]]),
                         )
                         ,tabPanel("Order", align = "center",
                                   plotOutput("rank3", height = height_plot[[3]]),
                         )
                         ,tabPanel("Family", align = "center",
                                   plotOutput("rank4", height = height_plot[[4]]),
                         )
                         ,tabPanel("Genus", align = "center",
                                   plotOutput("rank5", height = height_plot[[5]]),
                         )
                  )
                )
              })
              output_result_phylum <<- result_for_output(taxa.out, reg_result, Treatment, Treatment.Levels, Response, 1, num.taxa, covariates = Covariates, response.type = "continuous")
              output_result_class <<- result_for_output(taxa.out, reg_result, Treatment, Treatment.Levels, Response, 2, num.taxa, covariates = Covariates, response.type = "continuous")
              output_result_order <<- result_for_output(taxa.out, reg_result, Treatment, Treatment.Levels, Response, 3, num.taxa, covariates = Covariates, response.type = "continuous")
              output_result_family <<- result_for_output(taxa.out, reg_result, Treatment, Treatment.Levels, Response, 4, num.taxa, covariates = Covariates, response.type = "continuous")
              output_result_genus <<- result_for_output(taxa.out, reg_result, Treatment, Treatment.Levels, Response, 5, num.taxa, covariates = Covariates, response.type = "continuous")
              
              output$rank1 = try(renderPlot({ 
                output.to.show.2 (taxa.out, reg_result, output_result_phylum, Treatment, Treatment.Levels, Response, 1, num.taxa, legend.loc.int, Covariates, "continuous")
              }), silent = TRUE)
              
              output$rank2 = try(renderPlot({ 
                output.to.show.2 (taxa.out, reg_result, output_result_class, Treatment, Treatment.Levels, Response, 2, num.taxa, legend.loc.int, Covariates, "continuous")
              }), silent = TRUE)
              
              output$rank3 = try(renderPlot({ 
                output.to.show.2 (taxa.out, reg_result, output_result_order, Treatment, Treatment.Levels, Response, 3, num.taxa, legend.loc.int, Covariates, "continuous")
              }), silent = TRUE)
              
              output$rank4 = try(renderPlot({ 
                output.to.show.2 (taxa.out, reg_result, output_result_family, Treatment, Treatment.Levels, Response, 4, num.taxa, legend.loc.int, Covariates, "continuous")
              }), silent = TRUE)
              
              
              output$rank5 = try(renderPlot({ 
                output.to.show.2 (taxa.out, reg_result, output_result_genus, Treatment, Treatment.Levels, Response, 5, num.taxa, legend.loc.int, Covariates, "continuous")
              }), silent = TRUE)
              
              
            }
          }
          
          if(taxa.dataType =="clr"){
            output$int_reference <- renderUI({
              tagList(
                box(title = strong("Reference", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                    p(LM_REFERENCE_CLR, style = "font-size:11pt"))
              )
            })
          } else {
            output$int_reference <- renderUI({
              tagList(
                box(title = strong("Reference", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                    p(LM_REFERENCE_NO_CLR, style = "font-size:11pt"))
              )
            })
          }
        }
      }) 
    
    shinyjs::enable("run_int")
    shinyjs::enable("int_format")
    shinyjs::enable("int_response")
    shinyjs::enable("int_response_rename")
    shinyjs::enable("int_treat")
    shinyjs::enable("int_treat_rename")
    shinyjs::enable("int_cov")
    shinyjs::enable("int_chooseTestv")
    shinyjs::enable("int_num_taxa")
    shinyjs::enable("int_legend")
    
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  ## 4. Causal Forest -------------------
  observeEvent(input$cf_runButton, {
    shinyjs::disable("cf_format")
    shinyjs::disable("cf_response")
    shinyjs::disable("cf_response_rename")
    shinyjs::disable("cf_treat")
    shinyjs::disable("cf_treat_rename")
    shinyjs::disable("cf_method")
    shinyjs::disable("cf_cov")
    shinyjs::disable("cf_cov_list")
    shinyjs::disable("cf_tax_rank")
    shinyjs::disable("cf_num_taxa")
    shinyjs::disable("cf_download")
    shinyjs::disable("cf_reference")
    
    withProgress(
      message = "Calculation in progress",
      detail = "This may take a while...", value = 0, {
        if (input$cf_dataType == "Count (Rarefied)") {
          type = "rare.count"
          Feature <- as.data.frame(t(otu_table(infile$rare_biom)))
          if(input$cf_method == "Double-sample tree"){
            ref = DST_REFERENCE_RARE
          } else{
            ref = PT_REFERENCE_RARE
          }
        } else if (input$cf_dataType == "Proportion") {
          type = "prop"
          Feature <- as.data.frame(t(otu_table(infile$qc_biom))/colSums(otu_table(infile$qc_biom)))
          if(input$cf_method == "Double-sample tree"){
            ref = DST_REFERENCE
          } else {
            ref = PT_REFERENCE
          }
        } else if (input$cf_dataType == "CLR (Default)") {
          type = "clr"
          Feature <- as.data.frame(clr(t(otu_table(infile$qc_biom))+0.1))
          if(input$cf_method == "Double-sample tree"){
            ref = DST_REFERENCE_CLR
          } else {
            ref = PT_REFERENCE_CLR
          }
        } else if(input$cf_dataType == "Arcsine-root"){
          type = "arcsin"
          Feature <- as.data.frame(asin(sqrt(t(otu_table(infile$qc_biom))/colSums(otu_table(infile$qc_biom)))))
          if(input$cf_method == "Double-sample tree"){
            ref = DST_REFERENCE
          } else {
            ref = PT_REFERENCE
          }
        }
        colnames(Feature) <- paste("Feature", 1:ncol(Feature), sep = "")
        
        if(input$cf_tax_rank == "Phylum - Genus (Default)"){
          level.names = get.level.names(include = FALSE)
        } else {
          level.names = get.level.names(include = TRUE)
        }
        
        taxa.out <- chooseData$taxa.out
        sam.dat <- chooseData$sam.dat
        
        data <- taxa.out[[type]]
        colnames.list.mult <- colnames.to.ind(data)
        
        ### a. Double Sample Tree ----------
        if(input$cf_method == "Double-sample tree") {
          
          # When the response variable is binary
          if(response_is_binary$binary == TRUE) {
            Response <- as.factor(sam.dat[[input$cf_response]]) 
            Response <- as.numeric(Response) - 1
          }
          
          # When the response variable is continuous
          else{
            Response <- sam.dat[[input$cf_response]]
          }
          
          Treatment <- as.factor(sam.dat[[input$cf_treat]]) 
          Treatment.Levels <- substr(levels(Treatment), 1, 8)
          Treatment <- as.numeric(Treatment) - 1
          
          dt.fit.list <- list()
          rf.fit.list <- list()
          imp.plot.list <- list()
          pd.plot.list <- list()
          var.names.list <- list()
          dt.var.list <- list()
          rf.var.list <- list()
          width.list <- list()
          
          # Step 1. Treatment Effect Prediction
          incProgress(1/20, message = "Double Sample Tree: Treatment Effect Prediction in progress")
          set.seed(578)
          step1.result <- try(double.sample.treatment.pred(Feature, Response, Treatment, n.tree = 12000), silent = TRUE)
          
          for(name in level.names){
            # Step 2. Subgroup Identification
            incProgress(1/40, message = paste0(str_to_title(name), ": Subgroup Identification in progress"))
            step2.result <- try(subgroup.id(Treat.Effect = step1.result$Treat.Effect, taxa.out = taxa.out, type = type, level.name = name), silent = TRUE)
            
            var.names.list[[name]] <- data.frame(sub = colnames(step2.result$Taxa), ori = step2.result$taxa.names)
            
            # Step 3. BoRT
            # incProgress(1/40, message = paste0(str_to_title(name), ": BoRT in progress"))
            set.seed(578)
            step3.result <- try(bort.func(step2.result, name, n.tree = 100), silent = TRUE)
            
            # Step 4. Treatment Effect Prediction (randomForest)
            incProgress(1/40, message = paste0(str_to_title(name), ": Treatment Effect Prediction in progress"))
            set.seed(578)
            step4.result <- try(bort.treatment.pred(step2.result, name, n.tree = 100), silent = TRUE)
            
            dt.fit.list[[name]] <- step2.result$best.dt.fit
            # bort.out.list[[name]] <- step3.result$BoRT.out
            rf.fit.list[[name]] <- step4.result
            
            # DT used variable
            dt.var.list[[name]] <- try(cf.dt.used.var(step2.result), silent = TRUE)
            
            # RF used variable
            rf.var.list[[name]] <- try(rf.used.var(step2.result, step4.result), silent = TRUE)
            
            # Importance Plot
            incProgress(1/30, message = paste0(str_to_title(name), ": Preparing Feature Importance Plot"))
            imp.plot.list[[name]] <- try(cf.imp.plot(fit = step4.result,
                                                     type = 1,
                                                     subgroup.id.result = step2.result,
                                                     n = as.numeric(input$cf_num_taxa),
                                                     data.type = type, 
                                                     level.name = name), silent = TRUE)
            
            # Partial Dependence Plot
            incProgress(1/30, message = paste0(str_to_title(name), ": Preparing Partial Dependence Plot"))
            pd.plot.list[[name]] <- try(cf.pdp.reg(fit = step4.result,
                                                   X = step2.result$Taxa,
                                                   n = as.numeric(input$cf_num_taxa),
                                                   data.type = type), silent = TRUE)
            
            # PDP width setting
            if(as.numeric(input$cf_num_taxa) > ncol(step2.result$Taxa)){
              n <- ncol(step2.result$Taxa)
              if(n %% 5 == 0){
                n <- n / 5
              } else {
                n <- (n / 5) + 1
              }
            } else {
              n <- as.numeric(input$cf_num_taxa)
              n <- n / 5
            }
            width.list[[name]] <- n * 160
            
          }
        }
        
        ## b. Propensity Tree ----------
        else{
          
          # When the response variable is binary
          if(response_is_binary$binary == TRUE){
            Response <- as.factor(sam.dat[[input$cf_response]]) 
            Response <- as.numeric(Response) - 1
          }
          
          # When the response variable is continuous
          else{
            Response <- sam.dat[[input$cf_response]]
          }
          
          Treatment <- as.factor(sam.dat[[input$cf_treat]]) 
          Treatment.Levels <- substr(levels(Treatment), 1, 8)
          Treatment <- as.numeric(Treatment) - 1
          
          dt.fit.list <- list()
          rf.fit.list <- list()
          imp.plot.list <- list()
          pd.plot.list <- list()
          var.names.list <- list()
          dt.var.list <- list()
          rf.var.list <- list()
          imp.var.list <- list()
          width.list <- list()
          
          # Step 1-1. Treatment Effect Prediction with Covariate(s)
          if(!is.null(input$cf_cov_options)){
            cov <- input$cf_cov_options
            df <- data.frame(sam.dat[,cov])
            
            for(name in colnames(df)){
              dtype <- col.str.check(sam.dat, name)
              if(dtype == "factor"){
                df[[name]] <- as.factor(df[[name]])
              } else if(dtype == "numeric"){
                df[[name]] <- as.numeric(df[[name]])
              } else {
                df[[name]] <- df[[name]]
              }
            }
            
            f1 <- as.formula(paste("~" ,paste(names(df), collapse = "+"), collapse = ""))
            Covariate <- model.matrix(f1, data = df)[,-1]
            
            incProgress(1/20, message = "Propensity Tree with Covariate(s): Treatment Effect Prediction in progress")
            set.seed(578)
            step1.result <- try(propensity.treatment.pred(Feature, Response, Treatment, Covariate, n.tree = 12000), silent = TRUE)
          }
          
          # Step 1-2. Treatment Effect Prediction without Covariate
          
          else {
            incProgress(1/20, message = "Propensity Tree without Covariate: Treatment Effect Prediction in progress")
            set.seed(578)
            step1.result <- try(propensity.treatment.pred(Feature, Response, Treatment, n.tree = 12000), silent = TRUE)
          }
          
          for(name in level.names){
            
            # Step 2. Subgroup Identification
            incProgress(1/40, message = paste0(str_to_title(name), ": Subgroup Identification in progress"))
            set.seed(578)
            step2.result <- try(subgroup.id(Treat.Effect = step1.result$Treat.Effect, taxa.out = taxa.out, type = type, level.name = name), silent = TRUE)
            
            var.names.list[[name]] <- data.frame(sub = colnames(step2.result$Taxa), ori = step2.result$taxa.names)
            
            # Step 3. BoRT
            
            set.seed(578)
            step3.result <- try(bort.func(step2.result, name, n.tree = 100), silent = TRUE)
            
            # Step 4. Treatment Effect Prediction (randomForest)
            
            incProgress(1/40, message = paste0(str_to_title(name), ": Treatment Effect Prediction in progress"))
            set.seed(578)
            step4.result <- try(bort.treatment.pred(step2.result, name, n.tree = 100), silent = TRUE)
            
            dt.fit.list[[name]] <- step2.result$best.dt.fit
            
            rf.fit.list[[name]] <- step4.result
            
            # DT used variable
            
            dt.var.list[[name]] <- try(cf.dt.used.var(step2.result), silent = TRUE)
            
            # RF used variable
            
            rf.var.list[[name]] <- try(rf.used.var(step2.result, step4.result), silent = TRUE)
            
            # RF important variable
            
            imp.var.list[[name]] <- try(cf.imp.df(fit = step4.result, type = 0), silent = TRUE)
            
            # Importance Plot
            
            incProgress(1/30, message = paste0(str_to_title(name), ": Preparing Feature Importance Plot"))
            imp.plot.list[[name]] <- try(cf.imp.plot(fit = step4.result,
                                                     type = 1,
                                                     subgroup.id.result = step2.result,
                                                     n = as.numeric(input$cf_num_taxa),
                                                     data.type = type, 
                                                     level.name = name), silent = TRUE)
            
            # Partial Dependence Plot
            
            incProgress(1/30, message = paste0(str_to_title(name), ": Preparing Partial Dependence Plot"))
            pd.plot.list[[name]] <- try(cf.pdp.reg(fit = step4.result,
                                                   X = step2.result$Taxa,
                                                   n = as.numeric(input$cf_num_taxa),
                                                   data.type = type), silent = TRUE)
            
            # PDP width setting
            
            if(as.numeric(input$cf_num_taxa) > ncol(step2.result$Taxa)){
              n <- ncol(step2.result$Taxa)
              if(n %% 5 == 0){
                n <- n / 5
              }
              else{
                n <- (n / 5) + 1
              }
            }
            else{
              n <- as.numeric(input$cf_num_taxa)
              n <- n / 5
            }
            width.list[[name]] <- n * 160
          }
        }
        
        output$cf_display <- renderUI({
          tagList(
            tabBox(title = strong("Causal Machine Learning", style = "color:black", side = "right"), width = NULL,
                   tabPanel("Subgroup Identification", do.call(tabsetPanel, lapply(1:length(level.names), function(i){
                     tabPanel(str_to_title(level.names[i]), br(), 
                              box(title = strong("Subgroup Identification", style = "color:white"),  
                                  align = "center", width = NULL, solidHeader = TRUE,  status = "info", br(),
                                  plotOutput(paste0("DT", i), height = 700, width = 700), br(),
                                  dataTableOutput(paste0("DT_var_name", i), height = "auto", width = 600), br())
                     )
                   }))),
                   tabPanel("Treatment Effect Prediction", do.call(tabsetPanel, lapply(1:length(level.names), function(i){
                     tabPanel(str_to_title(level.names[i]), br(),
                              box(title = strong("Variable Importance", style = "color:white"),  
                                  align = "center", width = NULL, solidHeader = TRUE,  status = "info", br(),
                                  plotOutput(paste0("imp_plot", i), height = 700, width = 700)),
                              box(title = strong("Partial Dependence", style = "color:white"),  
                                  align = "center", width = NULL, solidHeader = TRUE,  status = "info", br(),
                                  plotOutput(paste0("pd_plot", i), height = 750, width = width.list[[i]]), br(),
                                  dataTableOutput(paste0("pd_var_name", i), height = "auto", width = 600))
                     )
                   })))
            )
          )
        })
        
        lapply(1:length(level.names), function(j){
          output[[paste0("DT",j)]] <- renderPlot({
            tryCatch(subgroup.id.vis(dt.fit.list[[j]]), error = function(e){
              message("Visualization not available! Check the input.")
              showModal(modalDialog(div("Visualization not available! Check the input.")))
              return(NULL)
            })
          })
          
          output[[paste0("DT_var_name",j)]] <- renderDataTable({
            tryCatch(dt.var.list[[j]], error = function(e){
              message("Visualization not available! Check the input.")
              showModal(modalDialog(div("Visualization not available! Check the input.")))
              return(NULL)
            })
          })
          
          output[[paste0("imp_plot", j)]] <- renderPlot({
            tryCatch(imp.plot.list[[j]], error = function(e){
              message("Visualization not available! Check the input.")
              showModal(modalDialog(div("Visualization not available! Check the input.")))
              return(NULL)
            })
          })
          
          output[[paste0("pd_plot", j)]] <- renderPlot({
            tryCatch(pd.plot.list[[j]], error = function(e){
              message("Visualization not available! Check the input.")
              showModal(modalDialog(div("Visualization not available! Check the input.")))
              return(NULL)
            })
          })
          
          output[[paste0("pd_var_name", j)]] <- renderDataTable({
            tryCatch(rf.var.list[[j]] %>% head(as.numeric(input$cf_num_taxa)), error = function(e){
              message("Visualization not available! Check the input.")
              showModal(modalDialog(div("Visualization not available! Check the input.")))
              return(NULL)
            })
          })
        })
        
        output$cf_download <- renderUI({
          tagList(
            p(" ", style = "margin-top: 20px;"),
            box(title = strong("Download Output", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                p("You can download the data analysis outputs.",
                  style = "font-size:11pt"),
                h5("Feature Importance"),
                downloadButton("rf_downloadTable", "Download", width = '50%', style = "background-color: red3")
            )
          )
        })
        
        output$rf_downloadTable <- downloadHandler(
          filename = function() {
            paste("Feature_Importance.zip")
          },
          content = function(DA.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            if (length(level.names) == 5) {
              dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt")
              for(i in 1:length(dataFiles)){
                write.table(try(cf.imp.df(fit = rf.fit.list[[i]], type = 0), silent = TRUE), file = dataFiles[i], sep = "\t")
              }
            }
            else if(length(level.names) == 6) {
              dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
              for(i in 1:length(dataFiles)){
                write.table(try(cf.imp.df(fit = rf.fit.list[[i]], type = 0), silent = TRUE), file = dataFiles[i], sep = "\t")
              }
            }
            zip(zipfile=DA.file, files=dataFiles)
          }
        )
        
        output$cf_reference <- renderUI({
          tagList(
            p(" ", style = "margin-top: 20px;"),
            box(title = strong("Reference", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                p(ref, style = "font-size:11pt"))
          )
        })
        
      }
    )
    
    shinyjs::enable("cf_format")
    shinyjs::enable("cf_response")
    shinyjs::enable("cf_response_rename")
    shinyjs::enable("cf_treat")
    shinyjs::enable("cf_treat_rename")
    shinyjs::enable("cf_method")
    shinyjs::enable("cf_cov")
    shinyjs::enable("cf_cov_list")
    shinyjs::enable("cf_tax_rank")
    shinyjs::enable("cf_num_taxa")
    shinyjs::enable("cf_download")
    shinyjs::enable("cf_reference")
  })
}

# RUN ------

shinyApp(ui = ui, server = server)
