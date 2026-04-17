###############################################################
### Script: run_n_cluster_cell_round0.R
### Purpose: Perform initial clustering and UMAP dimensionality reduction on ArchR projects for H3K27me3, H3K4me1, H3K4me2, and H3K4me3 datasets without doublet removal to visualize overall structure and identify potential clusters
### INPUT:
###   - ArchR projects (directories):
###   * ./tmp/H3K27me3_coCnT_all/
###   * ./tmp/H3K4me1_coCnT/
###   * ./tmp/H3K4me2_coCnT/
###   * ./tmp/H3K4me3_coCnT/
###
###  OUTPUT:
###    - Updated ArchR projects with clustering and UMAP embeddings:
###    * ./tmp/H3K27me3_coCnT_all/
###    * ./tmp/H3K4me1_coCnT/
###    * ./tmp/H3K4me2_coCnT/
###    * ./tmp/H3K4me3_coCnT/
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

source("../lib/R/lib_umap.R")

## Load project
proj_k27 = loadArchRProject("./tmp/H3K27me3_coCnT_all/")
proj_h3kme1 = loadArchRProject("./tmp/H3K4me1_coCnT/")
proj_h3kme2 = loadArchRProject("./tmp/H3K4me2_coCnT/")
proj_h3kme3 = loadArchRProject("./tmp/H3K4me3_coCnT/")

## Add TileMatrix
proj_k27 = addTileMatrix(proj_k27, tileSize = 50000, excludeChr = c("chrM", "chrY", "chrX"), force = T)
proj_h3kme1 = addTileMatrix(proj_h3kme1, tileSize = 5000, excludeChr = c("chrM", "chrY", "chrX"), force = T)
proj_h3kme2 = addTileMatrix(proj_h3kme2, tileSize = 5000, excludeChr = c("chrM", "chrY", "chrX"), force = T)
proj_h3kme3 = addTileMatrix(proj_h3kme3, tileSize = 5000, excludeChr = c("chrM", "chrY", "chrX"), force = T)

## LSI, UMAP
proj_k27 = run_cluster(proj_k27)
proj_h3kme1 = run_cluster(proj_h3kme1)
proj_h3kme2 = run_cluster(proj_h3kme2)
proj_h3kme3 = run_cluster(proj_h3kme3)

## Save project
saveArchRProject(proj_k27, outputDirectory = "./tmp/H3K27me3_coCnT_all/", load = T)
saveArchRProject(proj_h3kme1, outputDirectory = "./tmp/H3K4me1_coCnT/", load = T)
saveArchRProject(proj_h3kme2, outputDirectory = "./tmp/H3K4me2_coCnT/", load = T)
saveArchRProject(proj_h3kme3, outputDirectory = "./tmp/H3K4me3_coCnT/", load = T)

