##################################################
## Script: run_n_remove_bad_cluster2.R
## Purpose: Remove specific clusters identified as low-quality or doublet-enriched from ArchR projects, only remove C4 in K4me1 (hence the corresponding cells in K27me3) 
## INPUT:
## - ArchR projects (directories):
## * ./tmp/H3K27me3_coCnT_final1/
## * ./tmp/H3K4me1_coCnT_final1/
## * ./tmp/H3K4me2_coCnT_final1/
## * ./tmp/H3K4me3_coCnT_final1/
## OUTPUT:
## - Subsetted ArchR projects with specific clusters removed (directories):
## * ./tmp/H3K27me3_coCnT_final2/
## * ./tmp/H3K4me1_coCnT_final2/
## * ./tmp/H3K4me2_coCnT_final2/
## * ./tmp/H3K4me3_coCnT_final2/
##

library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

httpgd::hgd(port = 4324)

library(ArchR)
addArchRGenome("hg39")

proj_k4me3 = loadArchRProject("./tmp/H3K4me3_coCnT_final1/")
proj_k4me2 = loadArchRProject("./tmp/H3K4me2_coCnT_final1/")
proj_k4me1 = loadArchRProject("./tmp/H3K4me1_coCnT_final1/")
proj_k27 = loadArchRProject("./tmp/H3K27me3_coCnT_final1/")


## Remove C4 in K4me1
d_k4me1 = proj_k4me1@cellColData %>% as.data.table()
cells_to_remove_k4me1 = d_k4me1[Clusters == "C4"]
dim(cells_to_remove_k4me1)

d_k27 = proj_k27@cellColData %>% as.data.table()

cells_to_remove = cells_to_remove_k4me1$cell_id %>% unique()

subsetArchRProject(
	ArchRProj = proj_k4me1,
	cells = d_k4me1[!(cell_id %in% cells_to_remove), cell_id_ar],
	outputDirectory = "./tmp/H3K4me1_coCnT_final2/",
	dropCells = TRUE,
	force = TRUE
)

subsetArchRProject(
	ArchRProj = proj_k27,
	cells = d_k27[!(cell_id %in% cells_to_remove), cell_id_ar],
	outputDirectory = "./tmp/H3K27me3_coCnT_final2/",
	dropCells = TRUE,
	force = TRUE
)

saveArchRProject(proj_k4me2, outputDirectory = "./tmp/H3K4me2_coCnT_final2/", load = T)
saveArchRProject(proj_k4me3, outputDirectory = "./tmp/H3K4me3_coCnT_final2/", load = T)


