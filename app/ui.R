library(shiny)
library(BiocManager)
library(seqinr)
library(shinydashboard)
library(tidyverse)
library(plotly)
library(shinyWidgets)
library(shinyjs)
library(googleVis)
library(xtable)
library(DT)
library(htmltools)
library(phangorn)
library(bios2mds)
library(zip)
library(ape)
library(zCompositions)
library(compositions)
library(stringr)
library(caret)
library(ggplot2)
library(data.table)
library(fontawesome)
library(grid)
library(ggplotify)
library(remotes)
library(doParallel)

source("setSliderColor.R")
source("MiDataProc.Data.Upload.R")
source("MiDataProc.Data.Input.R")
source("MiDataProc.Descriptive.R")
source("MiDataProc.GLM.R")
source("MiDataProc.Beta.Diversity.R")
source("MiDataProc.ML.Models.R")
source("MiDataProc.Causal.R")

options(scipen=999)

# COMMENTS ------

{
  TITLE = p("MiCML: Microbiome Causal Machine Learning for the Analysis of Treatment Effects Using Microbial Profiles", style = "font-size:16pt")
  HOME_COMMENT_MV = p(strong("Importance:"), "The treatment effects are heterogenous by patients due to the differences 
                      in their microbiomes, which in turn implies that we can enhance the treatment effect by manipulating 
                      the patient’s microbiome profile. Then, the coadministration of microbiome-based dietary 
                      supplements (e.g., prebiotics, probiotics, dietary fiber) or therapeutics (e.g., antibiotics, pharmabiotics, 
                      phage therapy, microbiota transplantation) along with the primary treatment (e.g., immunotherapy) 
                      has been the subject of intensive investigation. For this, we first need to comprehend which 
                      microbes help (or prevent) the treatment to cure the patient’s disease, which is in principle the matter 
                      of interaction effects between treatment and microbiome on the patient’s recovery.", strong("MiCML (microbiome causal machine learning)"), "is the first 
                      cloud computing platform that streamlines related data processing and analytic procedures for the 
                      interaction effects on user-friendly web environments. MiCML is in particular unique with the up-to-date features 
                      of (i)", strong("batch effect correction"),  "to mitigate systematic variation in collective large-scale microbiome data 
                      due to the differences in their underlying batches (e.g., lab or study environments), 
                      and (ii)", strong("causal machine learning"), "to estimate treatment effects with consistency and then discern microbial 
                      taxa that enhance (or lower) the efficacy of the primary treatment. We also stress that MiCML 
                      can handle the data from either randomized controlled trials or observational studies. 
                      MiCML can be a useful analytic tool for microbiome-based personalized medicine to enhance patient 
                      well-beings while lowering medical expenses.", style = "font-size:12pt")
  
  HOME_COMMENT = p(strong("Description:"), "MiCML is a web cloud computing platform for the comprehensive analysis of treatment effects using microbiome profiles. MiCML consists of three", strong("Data Processing"), "modules, (i) Data Input, (ii) Batch Effect Correction & Quality Control, and (iii) Data Transformation and three", strong("Data Analysis"), "modules, 
                   (i) Descriptive Analysis, (ii) Generalized Linear Models, and (iii) Causal Machine Learning.", style = "font-size:12pt")
  
  HOME_COMMENT2 = p(strong("URLs:"), "Web server (online implementation):", tags$a(href = "http://micml.micloud.kr", "http://micml.micloud.kr"), 
                    "; GitHub repository (local implementation):", 
                    tags$a(href = "https://github.com/hk1785/micmlgit", "https://github.com/hk1785/micmlgit"), style = "font-size:12pt")
  
  HOME_COMMENT3 = p(strong("Maintainers:"), "Hyunwook Koh (", tags$a(href = "hyunwook.koh@stonybrook.edu", "hyunwook.koh@stonybrook.edu"), ")", style = "font-size:12pt")
  
  HOME_COMMENT4 = p(strong("Reference:"), "Koh H, Kim J, Jang H. 
                    Microbiome causal machine learning MiCML for the analysis of treatment effects using microbial profiles on user-friendly web environments (Submitted)", style = "font-size:12pt")
  
  INPUT_PHYLOSEQ_COMMENT1 = p(strong("Description:"), br(), br(), "This should be an '.rdata' or '.rds' file, and the data should be in 'phyloseq' format (see ", 
                              htmltools::a(tags$u("https://bioconductor.org/packages/release/bioc/html/phyloseq.html"), style = "color:red3"),
                              "). The phyloseq object should contain three necessary data, feature (OTU or ASV) table, taxonomic table and meta/sample information.",br(), br(),
                              strong("Details:"), br(), br(), 
                              strong("Feature table:"), "It should contain counts, where rows are features (OTUs or ASVs) and columns are units (row names are feature IDs and column names are unit IDs).", br(), br(),
                              strong("Taxonomic table:"),"It should contain taxonomic names, where rows are features and columns are seven taxonomic ranks (row names are feature IDs and column names are 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species' or 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species').", br(), br(),
                              strong("Metadata/Sample information:"),"It should contain variables for the units about host phenotypes, medical interventions, disease status or environmental/behavioral factors, where rows are units and columns are variables (row names are unit IDs, and column names are variable names).", br(), br(),
                              strong("Phylogenetic tree:"),"It should be a rooted tree. Otherwise, MiCML automatically roots the tree through midpoint rooting (phangorn::midpoint). The tip labels of the phylogenetic tree are feature IDs.", br(), br(),
                              "* The features should be matched and identical across feature table, taxonomic table and phylogenetic tree. 
                              The subjects should be matched and identical between feature table and metadata/sample information. MiCML will analyze only the matched features and subjects.", style = "font-size:11pt")
  
  INPUT_PHYLOSEQ_COMMENT2 = p("You can download example microbiome data 'Immuno.Metagenome.Char.Rdata' in 'phyloseq' format. For more details about 'phyloseq', see ", 
                              htmltools::a(tags$u("https://bioconductor.org/packages/release/bioc/html/phyloseq.html"), style = "color:red3"), br(), br(), 
                              "> setwd('/yourdatadirectory/')", br(), br(), 
                              "> load(file = 'Immuno.Metagenome.Char.Rdata')", br(), br(), 
                              "> library(phyloseq)", br(), br(), 
                              " > otu.tab <- otu_table(biom)", br(), 
                              " > tax.tab <- tax_table(biom)", br(), 
                              " > sam.dat <- sample_data(biom)", br(), 
                              " > tree <- phy_tree(biom)", br(),br(), 
                              "You can check if the features are matched and identical across feature table and taxonomic table, and the units are matched and identical between feature table and metadata/sample information using following code.", br(), br(), 
                              " > identical(rownames(otu.tab), rownames(tax.tab))", br(), 
                              " > identical(rownames(otu.tab), tree$tip.label)", br(),
                              " > identical(colnames(otu.tab), rownames(sam.dat))", style = "font-size:11pt", br(), br(),
                              strong("Reference:"), "Limeta A, Ji B, Levin M, Gatto F, Nielsen J. Meta-analysis of the gut microbiota in predicting response to cancer immunotherapy in metastatic melanoma. JCL Insight. 2020;5(23):e140940.")
  
  INPUT_INDIVIDUAL_DATA_COMMENT = p(strong("Description:"), br(), br(), 
                                    strong("Feature table:"), "It should contain counts, where rows are features (OTUs or ASVs) and columns are units (row names are feature IDs and column names are unit IDs).", br(), br(),
                                    strong("Taxonomic table:"),"It should contain taxonomic names, where rows are features and columns are seven taxonomic ranks (row names are feature IDs and column names are 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species' or 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species').", br(), br(),
                                    strong("Metadata/Sample information:"),"It should contain variables for the units about host phenotypes, medical interventions, disease status or environmental/behavioral factors, where rows are units and columns are variables (row names are unit IDs, and column names are variable names).", br(), br(),
                                    strong("Phylogenetic tree:"),"It should be a rooted tree. Otherwise, MiCML automatically roots the tree through midpoint rooting (phangorn::midpoint). The tip labels of the phylogenetic tree are feature IDs.", br(), br(),
                                    "* The features should be matched and identical across feature table, taxonomic table and phylogenetic tree.
                                    The subjects should be matched and identical between feature table and metadata/sample information. MiCML will analyze only the matched features and subjects.", style = "font-size:11pt")
  
  INPUT_INDIVIDUAL_DATA_COMMENT2 = p("You can download example microbiome data 'Immuno.Metagenome.Char.zip'. This zip file contains three necessary data components, feature table (otu.tab.txt), taxonomic table (tax.tab.txt), and metadata/sample information (sam.dat.txt).", br(), br(),
                                     "> setwd('/yourdatadirectory/')", br(), br(), 
                                     "> otu.tab <- read.table(file = 'otu.tab.txt', check.names = FALSE)", br(), 
                                     "> tax.tab <- read.table(file = 'tax.tab.txt', check.names = FALSE)", br(), 
                                     "> sam.dat <- read.table(file = 'sam.dat.txt', check.names = FALSE)", br(), 
                                     "> tree <- read.tree(file = 'tree.tre)", br(),br(),
                                     "You can check if the features are matched and identical across feature table and taxonomic table, 
                                     and the units are matched and identical between feature table and metadata/sample information using following code.", br(), br(), 
                                     " > identical(rownames(otu.tab), rownames(tax.tab))", br(), 
                                     " > identical(rownames(otu.tab), tree$tip.label)", br(),
                                     " > identical(colnames(otu.tab), rownames(sam.dat))", style = "font-size:11pt", br(), br(),
                                     strong("Reference:"), "Limeta A, Ji B, Levin M, Gatto F, Nielsen J. Meta-analysis of the gut microbiota in predicting response to cancer immunotherapy in metastatic melanoma. JCL Insight. 2020;5(23):e140940.")
  
  QC_KINGDOM_COMMENT = p("A microbial kingdom to be analyzed. Default is 'Bacteria' for 16S data. Alternatively, you can type 'Fungi' for ITS data 
                         or any other kingdom of interest for shotgun metagenomic data.", 
                         style = "font-size:11pt")
  
  QC_LIBRARY_SIZE_COMMENT1 = p("Remove units that have low library sizes (total read counts). Default is 3,000.", 
                               style = "font-size:11pt")
  
  QC_LIBRARY_SIZE_COMMENT2 = p("Library size: The total read count per unit.", 
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
  
  QC_BATCH_REFERENCE = p(h5(strong("ConQuR"), style = "margin-bottom: -15px"), br(),
                         "Ling W, Lu J, Zhao N. et al. Batch effects removal for microbiome data via conditional quantile regression. Nat Commun. 2022;13(5418)", br(),
                         h5(strong("PCoA"), style = "margin-bottom: -15px"), br(),
                         "Torgerson WS. Multidimensional scaling: I. Theorey and method. Psychometrika. 1952;17:401-419.", br(),
                         h5(strong("Bray-Curtis"), style = "margin-bottom: -15px"), br(),
                         "Bray JR, Curtis JT. An ordination of the upland forest communities of southern Wisconsin. Ecol Monogr. 1957;27:325-349.", br(),
                         h5(strong("CLR"), style = "margin-bottom: -15px"), br(),
                         "Aitchison J. The statistical analysis of compositional data. J R Stat Soc Series B Stat Methodol. 1982;44(2):139-160.")
  
  QC_BATCH_CONQUR_REFERENCE = p("1. Ling W, Lu J, Zhao N. et al. Batch effects removal for microbiome data via conditional quantile regression. Nat Commun. 2022;13(5418).", style = "font-size:11pt")
  
  QC_BATCH_COMBAT_REFERENCE = p("1. Zhang Y. et al. ComBat-seq: batch effect adjustments for RNA-seq count data. NAR Genom Bioinform. 2020;2(3)lqaa078", style = "font-size:11pt")
  
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
        ,logoBackColor = "rgb(230, 70, 50)"
        
        ,headerButtonBackColor = "rgb(230, 70, 50)"
        ,headerButtonIconColor = "rgb(230, 70, 50)"
        ,headerButtonBackColorHover = "rgb(230, 70, 50)"
        ,headerButtonIconColorHover = "rgb(0,0,0)"
        
        ,headerBackColor = "rgb(230, 70, 50)"
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
        ,boxInfoColor = "rgb(230, 70, 50)"
        ,boxSuccessColor = "rgb(112,173,71)"
        ,boxWarningColor = "rgb(244,156,104)"
        ,boxDangerColor = "rgb(255,88,55)"
        
        ,tabBoxTabColor = "rgb(255,255,255)"
        ,tabBoxTabTextSize = 14
        ,tabBoxTabTextColor = "rgb(0,0,0)"
        ,tabBoxTabTextColorSelected = "rgb(35, 49, 64)"
        ,tabBoxBackColor = "rgb(255,255,255)"
        ,tabBoxHighlightColor = "rgb(230, 70, 50)"
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
                div(id = "homepage", br(), HOME_COMMENT_MV, HOME_COMMENT, 
                    p(" ", style = "margin-bottom: 10px;"),
                    div(tags$img(src="Home2.png", height = 450, width = 480), style = "text-align: center;"), br(),
                    HOME_COMMENT2, HOME_COMMENT3, HOME_COMMENT4)),
        
        ## 0. DATA INPUT -----
        tabItem(tabName = "step1", br(),
                fluidRow(column(width = 6,
                                box(width = NULL, status = "info", solidHeader = TRUE,
                                    title = strong("Data Input", style = "color:white"),
                                    selectInput("inputOption", h5(strong("Data type")), 
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
                                  h5(strong("Batch effect correction?", style = "color:black")),
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
                                  textInput("kingdom", h5(strong("Kingdom")), value = "Bacteria"),
                                  QC_KINGDOM_COMMENT,
                                  tags$style(type = 'text/css', '#slider1 .irs-grid-text {font-size: 1px}'),
                                  tags$style(type = 'text/css', '#slider2 .irs-grid-text {font-size: 1px}'),
                                  
                                  sliderInput("slider1", h5(strong("Library size")), 
                                              min=0, max=10000, value = 3000, step = 1000),
                                  QC_LIBRARY_SIZE_COMMENT1,
                                  QC_LIBRARY_SIZE_COMMENT2,
                                  
                                  sliderInput("slider2", h5(strong("Mean proportion")), 
                                              min = 0, max = 0.1, value = 0.02, step = 0.001,  post  = " %"),
                                  QC_MEAN_PROP_COMMENT1,
                                  QC_MEAN_PROP_COMMENT2,
                                  
                                  br(),
                                  p(" ", style = "margin-bottom: -20px;"),
                                  
                                  h5(strong("Errors in taxonomic names")),
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
