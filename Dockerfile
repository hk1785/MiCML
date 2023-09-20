FROM openanalytics/r-base

MAINTAINER Hyunwook Koh "hyunwook.koh@stonybrook.edu"

RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libgsl-dev \
    libssh2-1-dev \
    libssl1.1 \
    libxml2-dev \
    build-essential \
    r-base-dev \
    pkg-config \
    cmake \
    libtiff5-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libgdal-dev \
    && rm -rf /var/lib/apt/lists/*
    
RUN apt-get update && apt-get install -y \
    libmpfr-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('ape', 'bios2mds', 'BiocManager', 'caret', 'compositions', 'data.table', 'doParallel', 'DT', 'ecodist', 'fontawesome', 'fossil', 'ggplot2', 'grid', 'ggplotify', 'googleVis', 'grf', 'GUniFrac', 'htmltools','MiRKAT', 'picante', 'plotly', 'phangorn', 'proxy', 'rmarkdown', 'randomForest', 'remotes', 'rpart', 'rpart.plot', 'reshape2', 'seqinr', 'shinydashboard', 'shiny', 'shinyWidgets', 'shinyjs', 'stringr', 'tidyverse', 'vegan', 'VGAM', 'xtable', 'zip', 'zCompositions'), repos='https://cloud.r-project.org/')"

RUN R -e "remotes::install_github('prise6/aVirtualTwins', build_vignettes = TRUE)"
RUN R -e "remotes::install_github('wdl2459/ConQuR')
RUN R -e "BiocManager::install('phyloseq')"
RUN R -e "remotes::install_github('joey711/biomformat')"
RUN R -e "remotes::install_github('jcrodriguez1989/chatgpt')"
RUN R -e "remotes::install_github('zmjones/edarf', subdir = 'pkg')"
RUN R -e "remotes::install_github('nik01010/dashboardthemes', force = TRUE)"
RUN R -e "remotes::install_github('hk1785/MiVT', force = TRUE)
RUN R -e "remotes::install_github('Zhiwen-Owen-Jiang/MiRKATMC', force = TRUE)
RUN R -e "BiocManager::install('sva')"

RUN mkdir /root/app
COPY app /root/app
COPY Rprofile.site /usr/lib/R/etc/

COPY app/Data/Immuno.Metagenome.Char.Rdata /root/app
COPY app/Data/otu_tab.txt /root/app
COPY app/Data/sam_dat.txt /root/app
COPY app/Data/tax_tab.txt /root/app
COPY app/Data/tree.tre /root/app

COPY app/www/Home2.png /root/app

COPY app/MiDataProc.Data.Input.R /root/app
COPY app/MiDataProc.Data.Upload.R /root/app
COPY app/MiDataProc.Description.R /root/app
COPY app/MiDataProc.GLM.R /root/app
COPY app/MiDataProc.ML.Models.R /root/app
COPY app/MiDataProc.Causal.R /root/app
COPY app/MiDataProc.Beta.Diversity.R /root/app

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/root/app')"]
