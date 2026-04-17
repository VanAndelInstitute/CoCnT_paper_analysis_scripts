##################################################
## Script: run_n_rna_map1.R
## Purpose: Map RNA-seq data to the clusters defined in the ArchR projects and visualize the mean expression of marker genes across clusters using a dot plot
## INPUT:
## - Seurat object containing RNA-seq data with cluster annotations:
## * ../data/triana_2021/WTA_projected.rds
## - Table of marker genes for each cluster and direction (up/down) for each dataset:
## * ./tmp/table_marker_genes_all_clusters_both_directions_round1_match.tsv
## OUTPUT:
## - Dot plot showing the mean scaled expression of marker genes across clusters:
## * ./figures/111_rna_map1.pdf
##

library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

library(Seurat)
library(ArchR)

color = ArchRPalettes$stallion
names(color) = paste0("C", 1:length(color))

so = read_rds("../data/triana_2021/WTA_projected.rds")
Idents(so) = so@meta.data$ct

so@meta.data$ct %>% table

d_maker_all_genes = fread("./tmp/table_marker_genes_all_clusters_both_directions_round1_match.tsv")

d_maker_gene = d_maker_all_genes[Direction == "Up" & marker == "K4me3" & Log2FC > 1.25, .(Cluster, gene_name = name, Log2FC, FDR)]

## choose the top 50 genes per cluster based on SD
d_maker_gene = d_maker_gene[, head(.SD, 50), by = .(Cluster)]

d_maker_gene[Cluster == "C3"]$gene_name

d_maker_gene$Cluster %<>% factor(., levels = names(color))
doup_gene = d_maker_gene[, .N, by = .(gene_name)][N > 1, gene_name]
d_maker_gene = d_maker_gene[!gene_name %in% doup_gene]
d_maker_gene = d_maker_gene[order(Cluster)]

## get genes of interest
genes_of_interest = d_maker_gene$gene_name

## get data matrix
expr_mat = GetAssayData(so, assay = "RNA", slot = "data") %>% as.matrix
dim(expr_mat)
genes_in_data = rownames(expr_mat)
expr_mat_sub = expr_mat[intersect(genes_of_interest, genes_in_data), ]
colnames(expr_mat_sub) = Idents(so)
# colnames(expr_mat_sub) = so@meta.data$ct
table(Idents(so), so@meta.data$ct)

d_plot = expr_mat_sub %>% melt %>% data.table
names(d_plot)[1:2] = c("gene_name", "cluster")
d_plot = d_plot[, .(mean_expr = mean(value)), by = .(gene_name, cluster)]

celltype_order <- c(
  # Stem & multipotent progenitors
  "HSCs & MPPs",
  "Lymphomyeloid prog",
  "Erythro-myeloid progenitors",
  
  # Erythroid / Megakaryocytic lineage
  "Early erythroid progenitor",
  "Late erythroid progenitor",
  "Aberrant erythroid",
  "Megakaryocyte progenitors",
  
  # Myeloid lineage
  "Early promyelocytes",
  "Late promyelocytes",
  "Myelocytes",
  "Eosinophil-basophil-mast cell progenitors",
  "Classical Monocytes",
  "Non-classical monocytes",
  
  # Dendritic lineage
  "Conventional dendritic cell 1",
  "Conventional dendritic cell 2",
  "Plasmacytoid dendritic cell progenitors",
  "Plasmacytoid dendritic cells",
  
  # B cell lineage
  "Pre-pro-B cells",
  "Pro-B cells",
  "Pre-B cells",
  "Small pre-B cell",
  "Immature B cells",
  "Mature naive B cells",
  "Nonswitched memory B cells",
  "Class switched memory B cells",
  "Plasma cells",
  
  # T cell lineage
  "CD4+ naive T cells",
  "CD4+ memory T cells",
  "CD69+PD-1+ memory CD4+ T cells",
  "CD8+ naive T cells",
  "CD8+ central memory T cells",
  "CD8+ effector memory T cells",
  "CD8+CD103+ tissue resident memory T cells",
  "GammaDelta T cells",
  "NK T cells",
  
  # NK lineage
  "NK cell progenitors",
  "CD56brightCD16- NK cells",
  "CD56dimCD16+ NK cells"
)

## scale by gene
d_plot[, mean_expr := scale(mean_expr, center = F), by = gene_name]
d_plot = merge(d_plot, d_maker_gene, by = "gene_name")
d_plot$gene_name = factor(d_plot$gene_name, levels = d_maker_gene$gene_name)
# d_plot$Cluster = factor(d_plot$Cluster, levels = 0:10)
# d_plot$cluster %<>% factor(levels = cluster_name[res$order])
d_plot$cluster %<>% factor(levels = rev(celltype_order))


d_plot_summary = d_plot[, .(index = mean(mean_expr)), by = .(Cluster, cluster)]

## Only average
ggplot(d_plot_summary, aes(y = cluster, x = Cluster, size = index, fill = Cluster)) +
  geom_point(shape = 21) +
  scale_size_continuous(range = c(0, 5)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cluster", y = "Reference", fill = "Mean Scaled Expression") +
  scale_fill_manual(values = color) +
  ggtitle("Mean Scaled Expression of Marker Genes Across Clusters")

ggsave(d"./figures/111_rna_map1.pdf", width = 8, height = 10)
