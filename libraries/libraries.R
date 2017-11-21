#libraries from CRAN
library(data.table)
library(igraph)
library(stringr)
library(parallel)
#libraries from bioconductor
library(fgsea)
library(rhdf5)
#this should be downloaded from the cmap L1000 repo and added  to the libs directory
source(file = "l1ktools/R/cmap/io.R") # from https://github.com/cmap/l1ktools
#These files are provided in this repo
source(file = "libraries/cMap_functions.R")
source(file = "libraries/RankMatrix_annotation_functions.R")
source(file = "libraries/functions_LINCS.R")
source(file = "libraries/functions_spearman.R")
source(file = "libraries/functions_rankedexp_fgsea_enrichment_graph.R")
source(file = "libraries/optimal_cut_function.R")
source(file = "libraries/handy_network_functions.R")
load(  file = "libraries/GOs_and_pathways.RData") #list of pathways 