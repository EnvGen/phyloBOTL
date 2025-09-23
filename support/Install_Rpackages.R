#!/usr/bin/env Rscript
#Usage: Rscript --vanilla Install_Rpackages.R <Ncpus>
args = commandArgs(trailingOnly=TRUE)

if (!is.na(args[1]) && args[1] > 1 ) CPUS=args[1] else CPUS=1

options(Ncpus = CPUS)
#Packages to be installed
CRAN_packages <- c("argparser","logistf","ape","tidyverse", "igraph","tidytree","gridExtra","future", "FactoMineR","phylolm", "factoextra","ggnewscale","randomcoloR",
                   "cluster","umap","tidyr","phangorn","pheatmap", "vegan","ggrepel", "data.table", "RColorBrewer", "devtools", "readr", "BiocManager","randomForest","pROC")

bioconduc_packages=c("ggtree")

github_packages= c("gggenomes") 

all_packages=c(CRAN_packages,bioconduc_packages,github_packages)


#Installation
new.packages <- CRAN_packages[!(CRAN_packages %in% installed.packages()[,"Package"])]
if(length(new.packages) >0 ) {
  for (s in new.packages) {
    if (!s %in% installed.packages()[,"Package"]) {
      
      cat("INFO: Installing :", s, "\n")
      suppressPackageStartupMessages(install.packages(s, dependencies=TRUE, 
                                                      Ncpus = CPUS,
                                                      quiet = T,
                                                      repos='http://cran.rstudio.com/'))}
  }
}
new.packages <- bioconduc_packages[!(bioconduc_packages %in% installed.packages()[,"Package"])]
if(length(new.packages) >0 ) {
  for (s in new.packages) {
    if (!s %in% installed.packages()[,"Package"]) {
      cat("INFO: Installing :", s, "\n")
      suppressPackageStartupMessages(BiocManager::install(s,
                                                          Ncpus = CPUS,
                                                          quiet = T
      ))}
  }
}

new.packages <- github_packages[!(github_packages %in% installed.packages()[,"Package"])]
if(length(new.packages) >0 ) {
  cat("INFO: Installing : gggenomes \n")
  remotes::install_github("thackl/thacklr",Ncpus = CPUS,
                          quiet = T)
  remotes::install_github("thackl/gggenomes",Ncpus = CPUS,
                          quiet = T)
}

# Final check 
new.packages <- all_packages[!(all_packages %in% installed.packages()[,"Package"])]
if(length(new.packages) >0 ) {
  cat("The following package(s) was(were) not installed:", new.packages, "\n")
} else {cat("All packages are installed \n")}
