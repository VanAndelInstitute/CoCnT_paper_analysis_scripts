##################################################
## Script: run_n_annotation_qc.R
## Purpose: Generate heatmaps to compare cluster annotations between H3K27me3 and H3K4me1, H3K4me2, H3K4me3 datasets to assess the consistency of clustering results across different histone modifications
## INPUT:
## - ArchR projects (directories):
## * ./tmp/H3K27me3_coCnT_final2/
## * ./tmp/H3K4me1_coCnT_final2/
## * ./tmp/H3K4me2_coCnT_final2/
## * ./tmp/H3K4me3_coCnT_final2/
## OUTPUT:
## - Heatmaps showing cluster correlations between H3K27me3 and each of the H3K4me1, H3K4me2, H3K4me3 datasets
## * ./figures/118_cluster_correlation_heatmap.pdf
##

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
color = ArchRPalettes$stallion
names(color) = paste0("C", (1:(length(color))))

proj_list = list(
  H3K27me3 = "./tmp/H3K27me3_coCnT_final2/",
  H3K4me1 = "./tmp/H3K4me1_coCnT_final2/",
  H3K4me2 = "./tmp/H3K4me2_coCnT_final2/",
  H3K4me3 = "./tmp/H3K4me3_coCnT_final2/"
)

## Heatmap shows how the cluster corrlate with each other annotations
proj_k27 = loadArchRProject(proj_list$H3K27me3)
proj_h3k4me1 = loadArchRProject(proj_list$H3K4me1)
proj_h3k4me2 = loadArchRProject(proj_list$H3K4me2)
proj_h3k4me3 = loadArchRProject(proj_list$H3K4me3)

plot_qc = function(proj1, proj2, main) {
    d1 = proj1@cellColData %>% as.data.table
    d2 = proj2@cellColData %>% as.data.table
    d_merge = merge(d1[, .(cell_id, Clusters)], d2[, .(cell_id, Clusters)], by = "cell_id", suffixes = c("_1", "_2"))
    m = table(d_merge$Clusters_1, d_merge$Clusters_2)
    ## order the cluster by name, C1, C2, instead of C1, C10, C2
    m = m[order(as.numeric(gsub("C", "", rownames(m)))), order(as.numeric(gsub("C", "", colnames(m))))]

    # pheatmap::pheatmap(m, scale = "row", cluster_rows = F, cluster_cols = F)
    m  = t(t(m) / colSums(m))
    pheatmap::pheatmap(m, scale = 'none', cluster_rows = F, cluster_cols = F, col = colorRampPalette(c("white", "Navy"))(100), main = main)
}

pdf("./figures/118_cluster_correlation_heatmap.pdf", width = 4, height = 3)

## K4me1 vs K27
plot_qc(proj_h3k4me1, proj_k27, main = "K4me1 vs K27")

## K4me2 vs K27
plot_qc(proj_h3k4me2, proj_k27, main = "K4me2 vs K27")

## K4me3 vs K27
plot_qc(proj_h3k4me3, proj_k27, main = "K4me3 vs K27")

dev.off()

## Remove later {{{
# d_plot_list = lapply(names(proj_list), function(target_x) {
#   proj = loadArchRProject(proj_list[[target_x]])
#   umap = getEmbedding(proj, embedding = "UMAP_Harmony", returnDF = TRUE)
#   names(umap) = c("UMAP1", "UMAP2")
#   d_meta = as.data.table(proj@cellColData)
#   d_plot = cbind(d_meta, umap)
#   d_plot[, target := target_x]
#   d_plot
# })
#
# d_plot = rbindlist(d_plot_list)
#
# ## Annotate antibody information
# d_plot[, antibody_target := paste0(antibody_target, "_", antibody_species)]
# d_plot[,
#   co_antibody_target := paste0(co_antibody_target, "_", co_antibody_species)
# ]
# d_plot[,
#   antibody_combination := paste0(antibody_target, "_", co_antibody_target)
# ]
# d_plot[
#   library_type %in% c("CoCnt1", "CoCnt2"),
#   antibody_target := paste0("z_", antibody_combination)
# ]
# table(d_plot$antibody_target)
#
# ## Color name
# annot_ct2 = c(
#   "C1" = "NK cells1",
#   "C2" = "T cells",
#   "C3" = "Mature B cells1",
#   "C4" = "Plasma cells",
#   "C5" = "Mature B cells2",
#   "C6" = "Memory B cells",
#   "C7" = "Aberrant erythroid",
#   "C8" = "Late erythroid progenitors",
#   "C9" = "NK cells2",
#   "C10" = "HSPC",
#   "C11" = "Plasmacytoid dendritic cells",
#   "C12" = "Myelocyte/Classical Monocytes1",
#   "C13" = "Myelocyte/Classical Monocytes2",
#   "C14" = "Myelocyte/Classical Monocytes3",
#   "C15" = "Non-classical Monocytes",
#   "C16" = "Pre-Pro-B cells",
#   "C17" = "Immature B cells",
#   "C18" = "Early Erythroid progenitors",
#   "C19" = "EoBaso progenitors"
# )
#
#
# ## Arrange the order of cell types
# ct_by_diff <- c(
#   "HSPC",
#   "EoBaso progenitors",
#   "Plasmacytoid dendritic cells",
#   "Myelocyte/Classical Monocytes1",
#   "Myelocyte/Classical Monocytes2",
#   "Myelocyte/Classical Monocytes3",
#   "Non-classical Monocytes",
#   "Erythroid progenitors1",
#   "Erythroid progenitors2",
#   "Erythroid progenitors3",
#   "NK cells1",
#   "NK cells2",
#   "T cells",
#   "Pre-Pro-B cells",
#   "Immature B cells",
#   "Mature B cells1",
#   "Mature B cells2",
#   "Memory B cells",
#   "Plasma cells"
# )
#
#
# ## Update the cell type annotation
# d_plot[, cell_type := factor(ct2, levels = ct_by_diff)]
#
#
# ## Update the color scheme
# names(color) = annot_ct2[names(color)]
#
# ## color antibody species
# ggplot(d_plot, aes(x = UMAP1, y = UMAP2, color = antibody_species)) +
#   geom_point(size = 0.5, alpha = 0.8) +
#   theme_classic() +
#   facet_wrap(~target) +
#   scale_color_nejm() +
#   theme(
#     legend.position = "right",
#     legend.title = element_blank()
#   )
#
# ## nFrags per cell by antibody: Move to the bed file QC
# ggplot(d_plot, aes(x = antibody_target, y = log10(nFrags))) +
#   geom_boxplot() +
#   theme_classic() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#   )
#
# ## Proportion of cell types per antibody target
# ggplot(
#   na.omit(d_plot[, .N, .(antibody_target, cell_type)]),
#   aes(x = antibody_target, y = N, fill = cell_type)
# ) +
#   geom_bar(stat = "identity", position = "fill") +
#   theme_classic() +
#   scale_fill_manual("fill", values = color) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   ) +
#   labs(x = "Antibody Target", y = "Proportion of Cell Types")
#
# ggplot(d_plot[!is.na(cell_type)], aes(x = cell_type, y = log10(nFrags))) +
#   geom_boxplot() +
#   facet_wrap(~antibody_target) +
#   theme_classic() +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1)
#   )
## }}}
