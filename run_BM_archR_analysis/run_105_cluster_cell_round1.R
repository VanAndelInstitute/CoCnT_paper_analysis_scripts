#############################################################
### Script: run_n_cluster_cell_round1.R
### Purpose: Perform clustering and UMAP dimensionality reduction on ArchR projects for H3K27me3, H3K4me1, H3K4me2, and H3K4me3 datasets after doublet removal to visualize overall structure and identify potential clusters
### INPUT:
###  - ArchR projects (directories):
###  * ./tmp/H3K27me3_coCnT_doubletRemoved15perc/
###  * ./tmp/H3K4me1_coCnT_doubletRemoved15perc
###  * ./tmp/H3K4me2_coCnT_doubletRemoved15perc/
###  * ./tmp/H3K4me3_coCnT_doubletRemoved15perc
### OUTPUT:
###  - Updated ArchR projects with clustering and UMAP embeddings:
###  * ./tmp/H3K27me3_coCnT_doubletRemoved15perc/
###  * ./tmp/H3K4me1_coCnT_doubletRemoved15perc
###  * ./tmp/H3K4me2_coCnT_doubletRemoved15perc/
###  * ./tmp/H3K4me3_coCnT_doubletRemoved15perc
###

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
proj_k27 = loadArchRProject("./tmp/H3K27me3_coCnT_doubletRemoved15perc/")
proj_h3kme1 = loadArchRProject("./tmp/H3K4me1_coCnT_doubletRemoved15perc/")
proj_h3kme2 = loadArchRProject("./tmp/H3K4me2_coCnT_doubletRemoved15perc/")
proj_h3kme3 = loadArchRProject("./tmp/H3K4me3_coCnT_doubletRemoved15perc/")

## LSI, UMAP
proj_k27 = run_cluster(proj_k27)
proj_h3kme1 = run_cluster(proj_h3kme1)
proj_h3kme2 = run_cluster(proj_h3kme2)
proj_h3kme3 = run_cluster(proj_h3kme3)

## Save project
saveArchRProject(proj_k27, outputDirectory = "./tmp/H3K27me3_coCnT_doubletRemoved15perc/", load = T)
saveArchRProject(proj_h3kme1, outputDirectory = "./tmp/H3K4me1_coCnT_doubletRemoved15perc/", load = T)
saveArchRProject(proj_h3kme2, outputDirectory = "./tmp/H3K4me2_coCnT_doubletRemoved15perc/", load = T)
saveArchRProject(proj_h3kme3, outputDirectory = "./tmp/H3K4me3_coCnT_doubletRemoved15perc/", load = T)

