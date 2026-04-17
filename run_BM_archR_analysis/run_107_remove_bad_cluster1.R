##################################################
## Script: run_n_remove_bad_cluster1.R
## Purpose: Remove specific clusters identified as low-quality or doublet-enriched from ArchR projects
## INPUT:
## - ArchR projects (directories):
## * ./tmp/H3K27me3_coCnT_doubletRemoved15perc/
## * ./tmp/H3K4me1_coCnT_doubletRemoved15perc/
## * ./tmp/H3K4me2_coCnT_doubletRemoved15perc/
## * ./tmp/H3K4me3_coCnT_doubletRemoved15perc/
## OUTPUT:
## - Subsetted ArchR projects with specific clusters removed (directories):
## * ./tmp/H3K27me3_coCnT_final1/
## * ./tmp/H3K4me1_coCnT_final1/
## * ./tmp/H3K4me2_coCnT_final1/
## * ./tmp/H3K4me3_coCnT_final1/
##
library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

library(ArchR)
addArchRGenome("hg38")

proj_k4me3 = loadArchRProject("./tmp/H3K4me3_coCnT_doubletRemoved15perc/")
proj_k4me2 = loadArchRProject("./tmp/H3K4me2_coCnT_doubletRemoved15perc/")
proj_k4me1 = loadArchRProject("./tmp/H3K4me1_coCnT_doubletRemoved15perc/")
proj_k27 = loadArchRProject("./tmp/H3K27me3_coCnT_doubletRemoved15perc/")


## Remove C2 in K27 and C10 in K4me3
d_k27 = proj_k27@cellColData %>% as.data.table()
cells_to_remove_k27 = d_k27[Clusters %in% c("C3", "C11")]
dim(cells_to_remove_k27)

d_k4me1 = proj_k4me1@cellColData %>% as.data.table()
cells_to_remove_k4me1 = d_k4me1[Clusters %in% c("C4")]
dim(cells_to_remove_k4me1)

d_k4me2 = proj_k4me2@cellColData %>% as.data.table()
cells_to_remove_k4me2 = d_k4me2[Clusters == "C10"]
dim(cells_to_remove_k4me2)

d_k4me3 = proj_k4me3@cellColData %>% as.data.table()
cells_to_remove_k4me3 = d_k4me3[Clusters == "C1"]
dim(cells_to_remove_k4me3)


## the cell need to remove in K27
cells_to_remove = unique(c(cells_to_remove_k27$cell_id, cells_to_remove_k4me2$cell_id, cells_to_remove_k4me3$cell_id))

subsetArchRProject(
    proj_k4me1, 
    cells = d_k4me1[!(cell_id %in% cells_to_remove), cell_id_ar],
    outputDirectory = "./tmp/H3K4me1_coCnT_final1/",
    dropCells = TRUE,
    force = TRUE
)
subsetArchRProject(
    proj_k4me2, 
    cells = d_k4me2[!(cell_id %in% cells_to_remove), cell_id_ar],
    outputDirectory = "./tmp/H3K4me2_coCnT_final1/",
    dropCells = TRUE,
    force = TRUE
)
subsetArchRProject(
    ArchRProj = proj_k4me3,
    cells = d_k4me3[!(cell_id %in% cells_to_remove), cell_id_ar],
    outputDirectory = "./tmp/H3K4me3_coCnT_final1/",
    dropCells = TRUE,
    force = TRUE
)
subsetArchRProject(
    ArchRProj = proj_k27,
    cells = d_k27[!(cell_id %in% cells_to_remove), cell_id_ar],
    outputDirectory = "./tmp/H3K27me3_coCnT_final1/",
    dropCells = TRUE,
    force = TRUE
)


