###################################################################
### Script: run_n_label_doublet.R
### Purpose: Label doublets in ArchR projects for H3K27me3, H3K4me1, H3K4me2, and H3K4me3
### INPUT:
###  - ArchR projects (directories):
###  * ./tmp/H3K27me3_coCnT_all/
###  * ./tmp/H3K4me1_coCnT/
###  * ./tmp/H3K4me2_coCnT/
###  * ./tmp/H3K4me3_coCnT/
###
### OUTPUT:
###  - Doublet annotation 
###  * ./tmp/table_H3K4me1_doublet15perc_label.tsv
###  * ./tmp/table_H3K4me2_doublet15perc_label.tsv
###  * ./tmp/table_H3K4me3_doublet15perc_label.tsv
###  * ./tmp/table_H3K27me3_doublet15perc_label.tsv
###  - figures double in UMAP
###  * ./figures/103_figure_label_doublet.pdf
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

source("../lib/R/lib_call_doublet.R")

proj_h3k27me3 = loadArchRProject("./tmp/H3K27me3_coCnT_all/")
proj_h3k4me3 = loadArchRProject("./tmp/H3K4me3_coCnT/")
proj_h3k4me2 = loadArchRProject("./tmp/H3K4me2_coCnT/")
proj_h3k4me1 = loadArchRProject("./tmp/H3K4me1_coCnT/")

## H3k4me1
pdf("./figures/103_figure_label_doublet.pdf", width = 6, height = 5)
rs = ct_label_doublet(proj_h3k4me1, doublet_percent_cuttoff = 15, redo_LSI = F)
rs$doublet_label %>% table
write_tsv(rs, "./tmp/table_H3K4me1_doublet15perc_label.tsv")

## H3k4me2
rs = ct_label_doublet(proj_h3k4me2, doublet_percent_cuttoff = 15, redo_LSI = F)
rs$doublet_label %>% table
write_tsv(rs, "./tmp/table_H3K4me2_doublet15perc_label.tsv")

## H3k4me3
rs = ct_label_doublet(proj_h3k4me3, doublet_percent_cuttoff = 15, redo_LSI = F)
rs$doublet_label %>% table
write_tsv(rs, "./tmp/table_H3K4me3_doublet15perc_label.tsv")

## H3k27me3
rs = ct_label_doublet(proj_h3k27me3, doublet_percent_cuttoff = 15, redo_LSI = F)
rs$doublet_label %>% table
write_tsv(rs, "./tmp/table_H3K27me3_doublet15perc_label.tsv")
dev.off()

