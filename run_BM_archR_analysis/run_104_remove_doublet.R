#######################################################
### Script: run_n_remove_doublet.R
### Purpose: Remove doublets from ArchR projects for H3K27me3, H3K4me1, H3K4me2, and H3K4me3
### INPUT:
### - ArchR projects (directories):
### * ./tmp/H3K27me3_coCnT_all/
### * ./tmp/H3K4me1_coCnT/
### * ./tmp/H3K4me2_coCnT/
### * ./tmp/H3K4me3_coCnT/
### - Doublet annotation
### * ./tmp/table_H3K4me1_doublet15perc_label.tsv
### * ./tmp/table_H3K4me2_doublet15perc_label.tsv
### * ./tmp/table_H3K4me3_doublet15perc_label.tsvl
### * ./tmp/table_H3K27me3_doublet15perc_label.tsv
### OUTPUT:
### - Subsetted ArchR projects with doublets removed (directories):
### * ./tmp/H3K27me3_coCnT_doubletRemoved15perc/
### * ./tmp/H3K4me1_coCnT_doubletRemoved15perc/
### * ./tmp/H3K4me2_coCnT_doubletRemoved15perc/
### * ./tmp/H3K4me3_coCnT_doubletRemoved15perc/
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

proj_k4me3 = loadArchRProject("./tmp/H3K4me3_coCnT/")
proj_k4me2 = loadArchRProject("./tmp/H3K4me2_coCnT/")
proj_k4me1 = loadArchRProject("./tmp/H3K4me1_coCnT/")
proj_k27 = loadArchRProject("./tmp/H3K27me3_coCnT_all/")

rs_k4me3 = fread("./tmp/table_H3K4me3_doublet15perc_label.tsv")
rs_k4me2 = fread("./tmp/table_H3K4me2_doublet15perc_label.tsv")
rs_k4me1 = fread("./tmp/table_H3K4me1_doublet15perc_label.tsv")
rs_k27 = fread("./tmp/table_H3K27me3_doublet15perc_label.tsv")

rs_k4me3$target = "H3K4me3"
rs_k4me2$target = "H3K4me2"
rs_k4me1$target = "H3K4me1"
rs_k27$target = "H3K27me3"

rs_all = rbind(rs_k4me3, rs_k4me2, rs_k4me1, rs_k27)

rs_all$doublet_label %>% table
removed_cell = rs_all[doublet_label %in% c("removed", "remove")]$cell_id

rs_all_removed = rs_all[(cell_id %in% removed_cell)]
rs_all_keep = rs_all[!(cell_id %in% removed_cell)]

pdf("./figures/104_figure_doublet_removed.pdf", width = 12, height = 3)
ggplot(rs_all_keep[target == "H3K27me3"]) + aes(x = UMAP1, y = UMAP2, color = status) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_color_nejm() +
    theme_classic() 

ggplot(rs_all_keep[target == "H3K4me1"]) + aes(x = UMAP1, y = UMAP2, color = status) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_color_nejm() +
    theme_classic() 

ggplot(rs_all_keep[target == "H3K4me2"]) + aes(x = UMAP1, y = UMAP2, color = status) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_color_nejm() +
    theme_classic() 

ggplot(rs_all_keep[target == "H3K4me3"]) + aes(x = UMAP1, y = UMAP2, color = status) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_color_nejm() +
    theme_classic() 
dev.off()


d_n_remove = rs_all_removed[, .(n_remove = .N), .(target)]
d_n_keep = rs_all_keep[, .(n_keep = .N), .(target)]
merge(d_n_remove, d_n_keep, by  = "target")

d_n_remove = rs_all_removed[, .(n_remove = .N), .(target, library_id, batch)]
d_n_keep = rs_all_keep[, .(n_keep = .N), .(target, library_id, batch)]
merge(d_n_remove, d_n_keep, by  = c("library_id", "target", "batch"))

dim(rs_all)
dim(rs_all_keep)

## Subset the ArchR project to remove doublets
subsetArchRProject(
    proj_k27,
    cells = rs_all_keep[target == "H3K27me3"]$cell_id_ar,
    outputDirectory = "./tmp/H3K27me3_coCnT_doubletRemoved15perc/",
    dropCells = TRUE,
    force = TRUE
)

subsetArchRProject(
    proj_k4me1,
    cells = rs_all_keep[target == "H3K4me1"]$cell_id_ar,
    outputDirectory = "./tmp/H3K4me1_coCnT_doubletRemoved15perc/",
    dropCells = TRUE,
    force = TRUE
)

subsetArchRProject(
    proj_k4me2,
    cells = rs_all_keep[target == "H3K4me2"]$cell_id_ar,
    outputDirectory = "./tmp/H3K4me2_coCnT_doubletRemoved15perc/",
    dropCells = TRUE,
    force = TRUE
)

subsetArchRProject(
    proj_k4me3,
    cells = rs_all_keep[target == "H3K4me3"]$cell_id_ar,
    outputDirectory = "./tmp/H3K4me3_coCnT_doubletRemoved15perc/",
    dropCells = TRUE,
    force = TRUE
)


