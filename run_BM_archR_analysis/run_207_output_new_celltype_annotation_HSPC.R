library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(ggsci)

httpgd::hgd(port = 4324)

library(ArchR)
addArchRGenome("hg38")
color = ArchRPalettes$stallion
names(color) = paste0("C", (1:(length(color))))


proj_k27 = loadArchRProject("./tmp/H3K27me3_coCnT_final3/")
proj_h3k4me1 = loadArchRProject("./tmp/H3K4me1_coCnT_final3/")
proj_h3k4me2 = loadArchRProject("./tmp/H3K4me2_coCnT_final3/")
proj_h3k4me3 = loadArchRProject("./tmp/H3K4me3_coCnT_final3/")

d_proj_k27 = as.data.table(proj_k27@cellColData)
d_proj_h3kme1 = as.data.table(proj_h3k4me1@cellColData)
d_proj_h3kme2 = as.data.table(proj_h3k4me2@cellColData)
d_proj_h3kme3 = as.data.table(proj_h3k4me3@cellColData)

d_cluster = d_proj_k27[, .(cell_id, cluster_k27 = ct3)]
d_cluster = merge(d_cluster, d_proj_h3kme1[, .(cell_id, cluster_h3kme1 = Clusters)], by = "cell_id", all = T)
d_cluster = merge(d_cluster, d_proj_h3kme2[, .(cell_id, cluster_h3kme2 = Clusters)], by = "cell_id", all = T)
d_cluster = merge(d_cluster, d_proj_h3kme3[, .(cell_id, cluster_h3kme3 = Clusters)], by = "cell_id", all = T)

write_tsv(d_cluster, "./tmp/table_cell_cluster_annotation_final_HSPC.tsv")

