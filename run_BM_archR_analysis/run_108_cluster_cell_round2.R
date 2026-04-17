##################################################
## Script: run_n_cluster_cell_round2.R
## Purpose: Perform clustering and UMAP dimensionality reduction on ArchR projects
## INPUT:
## - ArchR projects (directories):
## * ./tmp/H3K27me3_coCnT_final1/
## * ./tmp/H3K4me1_coCnT_final1/
## * ./tmp/H3K4me2_coCnT_final1/
## * ./tmp/H3K4me3_coCnT_final1/
## OUTPUT:
## - Updated ArchR projects with clustering and UMAP embeddings:
## * ./tmp/H3K27me3_coCnT_final1/
## * ./tmp/H3K4me1_coCnT_final1/
## * ./tmp/H3K4me2_coCnT_final1/
## * ./tmp/H3K4me3_coCnT_final1/
## - Figures of UMAP embeddings colored by clusters:
## * ./figures/108_cluster_cell_round2.pdf

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
proj_k27 = loadArchRProject("./tmp/H3K27me3_coCnT_final1/")
proj_h3kme1 = loadArchRProject("./tmp/H3K4me1_coCnT_final1/")
proj_h3kme2 = loadArchRProject("./tmp/H3K4me2_coCnT_final1/")
proj_h3kme3 = loadArchRProject("./tmp/H3K4me3_coCnT_final1/")

proj_k27 = run_cluster(proj_k27)
proj_h3kme1 = run_cluster(proj_h3kme1)
proj_h3kme2 = run_cluster(proj_h3kme2)
proj_h4kme3 = run_cluster(proj_h3kme3)


saveArchRProject(proj_k27, outputDirectory = "./tmp/H3K27me3_coCnT_final1/", load = T)
saveArchRProject(proj_h3kme1, outputDirectory = "./tmp/H3K4me1_coCnT_final1/", load = T)
saveArchRProject(proj_h3kme2, outputDirectory = "./tmp/H3K4me2_coCnT_final1/", load = T)
saveArchRProject(proj_h3kme3, outputDirectory = "./tmp/H3K4me3_coCnT_final1/", load = T)

pdf("./figures/108_cluster_cell_round2.pdf", width = 12, height = 10)
plotEmbedding(proj_k27, colorBy = "cellColData", name = "Clusters", embedding = "UMAP_Harmony") + theme(legend.position = "right")
plotEmbedding(proj_k27, colorBy = "cellColData", name = "Clusters", embedding = "UMAP") + theme(legend.position = "right")
plotEmbedding(proj_h3kme1, colorBy = "cellColData", name = "Clusters", embedding = "UMAP_Harmony") + theme(legend.position = "right")
plotEmbedding(proj_h3kme1, colorBy = "cellColData", name = "Clusters", embedding = "UMAP") + theme(legend.position = "right")
plotEmbedding(proj_h3kme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP_Harmony") + theme(legend.position = "right")
plotEmbedding(proj_h3kme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP") + theme(legend.position = "right")
plotEmbedding(proj_h3kme3, colorBy = "cellColData", name = "Clusters", embedding = "UMAP_Harmony") + theme(legend.position = "right")
plotEmbedding(proj_h3kme3, colorBy = "cellColData", name = "Clusters", embedding = "UMAP") + theme(legend.position = "right")
dev.off()
