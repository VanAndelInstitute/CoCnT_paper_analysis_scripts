##################################################
## Script: run_n_identify_bad_cluster_round1.R
## Purpose: Compare clusters across ArchR projects for H3K27me3, H3K4me1, H3K4me2, and H3K4me3 datasets
## INPUT:
## - ArchR projects (directories):
## * ./tmp/H3K27me3_coCnT_doubletRemoved15perc/
## * ./tmp/H3K4me1_coCnT_doubletRemoved15perc/
## * ./tmp/H3K4me2_coCnT_doubletRemoved15perc/
## * ./tmp/H3K4me3_coCnT_doubletRemoved15perc/
## OUTPUT:
## - QC plots comparing clusters across datasets (PDF):
## * ./figures/106_cluster_qc_round1.pdf
## - UMAP plots with cluster labels (PDF):
## * ./figures/106_cluster_umap_labels_round1.pdf
## - Doublet percentage per cluster (PDF):
## * ./figures/106_cluster_doublet_percentage_round1.pdf
## - Violin plots of QC metrics per cluster (PDF):
## * ./figures/106_cluster_qc_violin_round1.pdf

library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(ggsci)
library(ArchR)


color = ArchRPalettes$stallion
names(color) = paste0("C", seq_along(color))


library(ArchR)
addArchRGenome("hg38")

proj_k4me3 = loadArchRProject("./tmp/H3K4me3_coCnT_doubletRemoved15perc/")
proj_k4me2 = loadArchRProject("./tmp/H3K4me2_coCnT_doubletRemoved15perc/")
proj_k4me1 = loadArchRProject("./tmp/H3K4me1_coCnT_doubletRemoved15perc/")
proj_k27 = loadArchRProject("./tmp/H3K27me3_coCnT_doubletRemoved15perc/")

plot_qc = function(proj1, proj2) {
    d1 = proj1@cellColData %>% as.data.table
    d2 = proj2@cellColData %>% as.data.table
    d_merge = merge(d1[, .(cell_id, Clusters)], d2[, .(cell_id, Clusters)], by = "cell_id", suffixes = c("_1", "_2"))
    m = table(d_merge$Clusters_1, d_merge$Clusters_2)
    plotEmbedding(proj1, embedding = "UMAP_Harmony", colorBy = "cellColData", name = "Clusters") %>% print()
    plotEmbedding(proj2, embedding = "UMAP_Harmony", colorBy = "cellColData", name = "Clusters") %>% print()
    plotEmbedding(proj1, colorBy = "cellColData", name = "Clusters") %>% print()
    plotEmbedding(proj2, colorBy = "cellColData", name = "Clusters") %>% print()
    pheatmap::pheatmap(m, scale = "row", cluster_rows = F, cluster_cols = F)
    pheatmap::pheatmap(m, scale = "column", cluster_rows = F, cluster_cols = F)
}

pdf("./figures/106_qc_plots_round1.pdf", width = 10, height = 8)

## K4me1 vs K27
plot_qc(proj_k4me1, proj_k27)

## K4me2 vs K27
plot_qc(proj_k4me2, proj_k27)

## K4me3 vs K27
plot_qc(proj_k4me3, proj_k27)

dev.off()

## Dig into C2 of K27 and C10 in K4me3
pdf("./figures/106_cluster_qc_round1.pdf", width = 6, height = 8)
# plotTSSEnrichment(proj_k27, groupBy = "Clusters")
plotFragmentSizes(proj_k27, groupBy = "Clusters")
# plotTSSEnrichment(proj_k4me1, groupBy = "Clusters")
plotFragmentSizes(proj_k4me1, groupBy = "Clusters")
# plotTSSEnrichment(proj_k4me2, groupBy = "Clusters")
plotFragmentSizes(proj_k4me2, groupBy = "Clusters")
# plotTSSEnrichment(proj_k4me3, groupBy = "Clusters")
plotFragmentSizes(proj_k4me3, groupBy = "Clusters")
dev.off()

## UMAP with cluster labels, no legend
pdf("./figures/106_cluster_umap_labels_round1.pdf", width = 6, height = 4)
plotEmbedding(proj_k4me1, embedding = "UMAP_Harmony", colorBy = "cellColData", name = "Clusters", labelAsFactor = T, labelMeans = T, size = 0.5)
plotEmbedding(proj_k4me2, embedding = "UMAP_Harmony", colorBy = "cellColData", name = "Clusters", labelAsFactor = T, labelMeans = T, size = 0.5)
plotEmbedding(proj_k4me3, embedding = "UMAP_Harmony", colorBy = "cellColData", name = "Clusters", labelAsFactor = T, labelMeans = T, size = 0.5)
plotEmbedding(proj_k27, embedding = "UMAP_Harmony", colorBy = "cellColData", name = "Clusters", labelAsFactor = T, labelMeans = T, size = 0.5)
dev.off()

## Doublet percentage
# Helper function to summarize and plot doublet/unassigned fraction and cell number per cluster
plot_cluster_summary <- function(cell_col_data, title = NULL) {
  d <- as.data.table(cell_col_data)
  
  # Summaries per cluster
  d_doublet <- d[, .(
    doublet_percent = sum(status %in% c("doublet", "unassigned"), na.rm = TRUE) / .N
  ), by = Clusters]

  
  d_cells <- d[, .(
    cell_number = .N
  ), by = Clusters]
  
  # Merge summaries
  d_plot <- merge(d_cells, d_doublet, by = "Clusters", all = TRUE)
  
  # Order clusters by cell number for readability
  d_plot[, Clusters := factor(Clusters, levels = d_plot[order(-cell_number)]$Clusters)]
  
  # Scale factor to map doublet% onto the left axis range
  sf <- max(d_plot$cell_number, na.rm = TRUE) / max(d_plot$doublet_percent, na.rm = TRUE)

  d_plot$Clusters <- factor(d_plot$Clusters, levels = paste0("C", seq_len(nrow(d_plot))))
  
  # Plot: bars = doublet% (right axis), points = cell number (left axis)
  ggplot(d_plot, aes(x = Clusters)) +
    geom_col(aes(y = doublet_percent * sf, fill = Clusters), alpha = 1) +
    geom_point(aes(y = cell_number), size = 2, color = "blue") +
    scale_y_continuous(
      name = "Cell number",
      sec.axis = sec_axis(~ . / sf, name = "Doublet fraction")
    ) +
    scale_fill_manual(values = color) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y.right = element_text(margin = margin(l = 8), color = "blue"),
    ) +
    (if (!is.null(title)) ggtitle(title) else NULL)
}


pdf("./figures/106_cluster_doublet_percentage_round1.pdf", width = 5, height = 3)
# Plot for proj_k4me1
plot_cluster_summary(proj_k4me1@cellColData, title = "proj_k4me1")

# Plot for proj_k4me2
plot_cluster_summary(proj_k4me2@cellColData, title = "proj_k4me2")

# Plot for proj_k4me3
plot_cluster_summary(proj_k4me3@cellColData, title = "proj_k4me3")

# Plot for proj_k27
plot_cluster_summary(proj_k27@cellColData, title = "proj_k27")

dev.off()


## Violin plots of QC metrics
pdf("./figures/106_cluster_qc_violin_round1.pdf", width = 8, height = 6)
plotGroups(proj_k4me1, groupBy = "Clusters", colorBy = "colData", name = "TSSEnrichment", plotAs = "violin", alpha = 0.5)
plotGroups(proj_k4me2, groupBy = "Clusters", colorBy = "colData", name = "TSSEnrichment", plotAs = "violin", alpha = 0.5)
plotGroups(proj_k4me3, groupBy = "Clusters", colorBy = "colData", name = "TSSEnrichment", plotAs = "violin", alpha = 0.5)
plotGroups(proj_k27, groupBy = "Clusters", colorBy = "colData", name = "TSSEnrichment", plotAs = "violin", alpha = 0.5)

plotGroups(proj_k4me2, groupBy = "Clusters", colorBy = "colData", name = "nFrags", plotAs = "violin", alpha = 0.5, log2Norm = T)
plotGroups(proj_k4me3, groupBy = "Clusters", colorBy = "colData", name = "nFrags", plotAs = "violin", alpha = 0.5, log2Norm = T)
plotGroups(proj_k27, groupBy = "Clusters", colorBy = "colData", name = "nFrags", plotAs = "violin", alpha = 0.5, log2Norm = T)
dev.off()
