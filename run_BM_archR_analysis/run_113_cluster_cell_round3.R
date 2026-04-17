##################################################
## Script: run_n_cluster_cell_round3.R
## Purpose: Perform clustering and UMAP dimensionality reduction on ArchR projects for H3K4me1, and H3K27me3 datasets
## INPUT:
## - ArchR projects (directories):
## * ./tmp/H3K27me3_coCnT_final2/
## * ./tmp/H3K4me1_coCnT_final2/
## OUTPUT:
## - Updated ArchR projects with clustering and UMAP embeddings:
## * ./tmp/H3K27me3_coCnT_final2/
## * ./tmp/H3K4me1_coCnT_final2/
## 

library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(ggsci)

library(ArchR)
addArchRGenome("hg38")
addArchRThreads(threads = 11)
color = ArchRPalettes$stallion

source("../../lib/R/lib_umap.R")

## Load project
proj_k27 = loadArchRProject("./tmp/H3K27me3_coCnT_final2/")
proj_h3kme1 = loadArchRProject("./tmp/H3K4me1_coCnT_final2/")

proj_k27 = run_cluster(proj_k27)
proj_h3kme1 = run_cluster(proj_h3kme1)

saveArchRProject(proj_k27, outputDirectory = "./tmp/H3K27me3_coCnT_final2/", load = T)
saveArchRProject(proj_h3kme1, outputDirectory = "./tmp/H3K4me1_coCnT_final2/", load = T)


