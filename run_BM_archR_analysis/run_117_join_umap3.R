##################################################
## Script: run_n_join_umap3.R
## Purpose: Join UMAP embeddings from multiple ArchR projects, visualize them, and save updated ArchR projects with combined cluster annotations
## INPUT:
## - ArchR projects (directories):
## * ./tmp/H3K27me3_coCnT_final2/
## * ./tmp/H3K4me1_coCnT_final2/
## * ./tmp/H3K4me2_coCnT_final2/
## * ./tmp/H3K4me3_coCnT_final2/
## OUTPUT:
## - Annotation table with combined cluster annotations for all projects:
## * ./tmp/table_cell_cluster_annotation_final_round2.tsv
## - Combined UMAP embeddings and cluster annotations for all projects
## * ./figures/117_figure_UMAP_coCnt_final_round2_cluster_k27.pdf
## * ./figures/117_figure_UMAP_coCnt_final_round2_cluster_k27_line.pdf
## * ./figures/117_figure_UMAP_coCnt_final_round2_cluster_k27_no_line.pdf
## * ./figures/117_figure_UMAP_coCnt_final_round2_cluster_h3k4me1.pdf
## * ./figures/117_figure_UMAP_coCnt_final_round2_cluster_h3k4me2.pdf
## * ./figures/117_figure_UMAP_coCnt_final_round2_cluster_h3k4me3.pdf
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

color <- c(
    C1  = "#8C8178",
    C2  = "#A3623A",
    C3  = "#B4009E",
    C4  = "#7D0033",
    C5  = "#89288F",
    C6  = "#B6B0FF",
    C7  = "#499195",
    C8  = "#208A42",
    C9  = "#4A5B5D",
    C10 = "#D51F26",
    C11 = "#FFBBED",
    C12 = "#8B3E1D",
    C13 = "#F47D2B",
    C14 = "#E09E5A",
    C15 = "#BDA687",
    C16 = "#6C2B85",
    C17 = "#8250C4",
    C18 = "#14532D",
    C19 = "#8AD4EB"
)

proj_k27 = loadArchRProject("./tmp/H3K27me3_coCnT_final2/")
proj_h3k4me1 = loadArchRProject("./tmp/H3K4me1_coCnT_final2/")
proj_h3k4me2 = loadArchRProject("./tmp/H3K4me2_coCnT_final2/")
proj_h3k4me3 = loadArchRProject("./tmp/H3K4me3_coCnT_final2/")

d_proj_k27 = as.data.table(proj_k27@cellColData)
d_proj_h3kme1 = as.data.table(proj_h3k4me1@cellColData)
d_proj_h3kme2 = as.data.table(proj_h3k4me2@cellColData)
d_proj_h3kme3 = as.data.table(proj_h3k4me3@cellColData)

d_cluster = d_proj_k27[, .(cell_id, cluster_k27 = Clusters)]
d_cluster = merge(d_cluster, d_proj_h3kme1[, .(cell_id, cluster_h3kme1 = Clusters)], by = "cell_id", all = T)
d_cluster = merge(d_cluster, d_proj_h3kme2[, .(cell_id, cluster_h3kme2 = Clusters)], by = "cell_id", all = T)
d_cluster = merge(d_cluster, d_proj_h3kme3[, .(cell_id, cluster_h3kme3 = Clusters)], by = "cell_id", all = T)


# Level 2
annot_ct2 = c(
    "C1" = "NK cells1",
    "C2" = "T cells",
    "C3" = "Mature B cells1",
    "C4" = "Plasma cells",
    "C5" = "Mature B cells2",
    "C6" = "Memory B cells",
    "C7" = "Erythroid progenitors3",
    "C8" = "Erythroid progenitors2",
    "C9" = "NK cells2",
    "C10" = "HSPC",
    "C11" = "Plasmacytoid dendritic cells",
    "C12" = "Myelocyte/Classical Monocytes1",
    "C13" = "Myelocyte/Classical Monocytes2",
    "C14" = "Myelocyte/Classical Monocytes3",
    "C15" = "Non-classical Monocytes",
    "C16" = "Pre-Pro-B cells",
    "C17" = "Immature B cells",
    "C18" = "Erythroid progenitors1",
    "C19" = "EoBaso progenitors"
)


color_ct2 = color
names(color_ct2) = annot_ct2

ordered_ct <- c(
    "HSPC",
    "EoBaso progenitors",
    "Plasmacytoid dendritic cells",
    "Myelocyte/Classical Monocytes1",
    "Myelocyte/Classical Monocytes2",
    "Myelocyte/Classical Monocytes3",
    "Non-classical Monocytes",
    "Erythroid progenitors1",
    "Erythroid progenitors2",
    "Erythroid progenitors3",
    "NK cells1",
    "NK cells2",
    "T cells",
    "Pre-Pro-B cells",
    "Immature B cells",
    "Mature B cells1",
    "Mature B cells2",
    "Memory B cells",
    "Plasma cells"
)

annot_ct2 %>% as.character

d_cluster$ct2 = revalue(d_cluster$cluster_k27, replace = annot_ct2)

## save annotation table
write_tsv(d_cluster, "./tmp/table_cell_cluster_annotation_final_round2.tsv")


umap_k27 = as.data.table( getEmbedding( ArchRProj = proj_k27, embedding = "UMAP_Harmony"))
umap_h3kme1 = as.data.table( getEmbedding( ArchRProj = proj_h3k4me1, embedding = "UMAP_Harmony"))
umap_h3kme2 = as.data.table( getEmbedding( ArchRProj = proj_h3k4me2, embedding = "UMAP_Harmony"))
umap_h3kme3 = as.data.table( getEmbedding( ArchRProj = proj_h3k4me3, embedding = "UMAP_Harmony"))

names(umap_k27) = c("UMAP1", "UMAP2")
names(umap_h3kme1) = c("UMAP1", "UMAP2")
names(umap_h3kme2) = c("UMAP1", "UMAP2")
names(umap_h3kme3) = c("UMAP1", "UMAP2")

scale2one = function(x) {
    x / (max(x) - min(x)) / 1.5
}

## Plot triple UMAP
u = 0.85

umap_k27[, UMAP1 := scale2one(UMAP1)]
umap_k27[, UMAP2 := scale2one(UMAP2)]
umap_h3kme1[, UMAP1 := scale2one(UMAP1) - 0.8 * u]
umap_h3kme1[, UMAP2 := scale2one(UMAP2) - 0.8 * u]
umap_h3kme2[, UMAP1 := scale2one(UMAP1) + 0 * u]
umap_h3kme2[, UMAP2 := scale2one(UMAP2) + 1.1 * u]
umap_h3kme3[, UMAP1 := scale2one(UMAP1) + 0.8 * u]
umap_h3kme3[, UMAP2 := scale2one(UMAP2) - 0.8 * u]

umap_k27[, cell_id := d_proj_k27$cell_id]
umap_h3kme1[, cell_id := d_proj_h3kme1$cell_id]
umap_h3kme2[, cell_id := d_proj_h3kme2$cell_id]
umap_h3kme3[, cell_id := d_proj_h3kme3$cell_id]

umap_k27[, mark := "H3K27me3"]
umap_h3kme1[, mark := "H3K4me1"]
umap_h3kme2[, mark := "H3K4me2"]
umap_h3kme3[, mark := "H3K4me3"]


d_plot = rbind(umap_k27, umap_h3kme1, umap_h3kme2, umap_h3kme3)
d_plot = merge(d_plot, d_cluster, by = "cell_id")

color_ct2 = color
names(color_ct2) = annot_ct2
d_plot$ct2 %<>% factor(levels = ordered_ct)


d_plot$cluster_k27 %<>% factor(levels = c(paste0("C", 1:length(unique(d_plot$cluster_k27)))))
g = ggplot(d_plot[, .(UMAP1, UMAP2, ct2, cell_id)] %>% na.omit, aes(x = UMAP1, y = UMAP2, color = ct2)) +
    geom_point(size = 0.50, alpha = 0.3, stroke = 0) +
    scale_color_manual(values = color_ct2) +
    guides(color = guide_legend(override.aes = list(size=5, alpha = 0.8))) +
    theme_classic() +
    theme(
	legend.position = "right",
	axis.text = element_blank(),
	axis.ticks = element_blank(),
	axis.title = element_blank()
    )
g
ggsave("./figures/117_figure_UMAP_coCnt_final_round2_cluster_k27.pdf", g, width = 9, height = 7)

g_line = ggplot(d_plot[, .(UMAP1, UMAP2, ct2, cell_id)] %>% na.omit, aes(x = UMAP1, y = UMAP2, color = ct2)) +
    geom_line(aes(group = cell_id), alpha = 0.01, size = 0.001) +
    geom_point(size = 0.5, alpha = 0.3, stroke = 0) +
    scale_color_manual(values = color_ct2) +
    guides(color = guide_legend(override.aes = list(size=5, alpha = 0.8))) +
    theme_classic() +
    theme(
	legend.position = "right",
	axis.text = element_blank(),
	axis.ticks = element_blank(),
	axis.title = element_blank()
    )

g_line
ggsave("./figures/117_figure_UMAP_coCnt_final_round2_cluster_k27_line.pdf", g_line, width = 9, height = 7)

cdolor_ct2 = color
names(color_ct2) = annot_ct2
d_plot$ct2 %<>% factor(levels = ordered_ct)

ggplot(d_plot[, .(UMAP1, UMAP2, ct2, cell_id)] %>% na.omit, aes(x = UMAP1, y = UMAP2, color = ct2)) +
    # geom_line(aes(group = cell_id), alpha = 0.03, size = 0.05) +
    geom_point(size = 0.3, alpha = 0.8) +
    scale_color_manual(values = color_ct2) +
    theme_classic() +
    guides(color = guide_legend(override.aes = list(size=10))) +
    theme(
	legend.position = "right",
	axis.text = element_blank(),
	axis.ticks = element_blank(),
	axis.title = element_blank()
    )
ggsave("./figures/117_figure_UMAP_coCnt_final_round2_cluster_k27_no_line.pdf", width = 9, height = 7)


## Plot dual UMAP
# d_plot$cluster_h3kme1 %<>% factor(levels = c(paste0("C", 1:length(unique(d_plot$cluster_h3kme1)))), NA)
ggplot(d_plot[rev(order(cluster_h3kme1))], aes(x = UMAP1, y = UMAP2, color = cluster_h3kme1)) +
    # geom_line(aes(group = cell_id), alpha = 0.03, size = 0.05) +
    geom_point(size = 0.05, alpha = 0.3) +
    scale_color_manual(values = color) +
    theme_classic() +
    theme(
	legend.position = "right",
	axis.text = element_blank(),
	axis.ticks = element_blank(),
	axis.title = element_blank()
    )
ggsave("./figures/117_figure_UMAP_coCnt_final_round2_cluster_h3kme1.pdf", width = 9, height = 7)


ggplot(d_plot[rev(order(cluster_h3kme2))], aes(x = UMAP1, y = UMAP2, color = cluster_h3kme2)) +
    # geom_line(aes(group = cell_id), alpha = 0.03, size = 0.05) +
    geom_point(size = 0.05, alpha = 0.3) +
    scale_color_manual(values = color) +
    theme_classic() +
    theme(
	legend.position = "right",
	axis.text = element_blank(),
	axis.ticks = element_blank(),
	axis.title = element_blank()
    )
ggsave("./figures/117_figure_UMAP_coCnt_final_round2_cluster_h3kme2.pdf", width = 9, height = 7)

ggplot(d_plot[rev(order(cluster_h3kme3))], aes(x = UMAP1, y = UMAP2, color = cluster_h3kme3)) +
    # geom_line(aes(group = cell_id), alpha = 0.03, size = 0.05) +
    geom_point(size = 0.05, alpha = 0.3) +
    scale_color_manual(values = color) +
    theme_classic() +
    theme(
	legend.position = "right",
	axis.text = element_blank(),
	axis.ticks = element_blank(),
	axis.title = element_blank()
    )
ggsave("./figures/117_figure_UMAP_coCnt_final_round2_cluster_h3kme3.pdf", width = 9, height = 7)

## Update achR proj
v_cluster_plus = d_cluster$ct2[match(proj_k27@cellColData$cell_id, d_cluster$cell_id)]
proj_k27 = addCellColData(ArchRProj = proj_k27, data = v_cluster_plus, name = "ct2", cells = proj_k27@cellColData$cell_id_ar, force = T)
saveArchRProject(ArchRProj = proj_k27, outputDirectory = "./tmp/H3K27me3_coCnT_final2/", load = F)

v_cluster_plus = d_cluster$ct2[match(proj_h3k4me1@cellColData$cell_id, d_cluster$cell_id)]
proj_h3k4me1 = addCellColData(ArchRProj = proj_h3k4me1, data = v_cluster_plus, name = "ct2", cells = proj_h3k4me1@cellColData$cell_id_ar, force = T)
saveArchRProject(ArchRProj = proj_h3k4me1, outputDirectory = "./tmp/H3K4me1_coCnT_final2/", load = F)

v_cluster_plus = d_cluster$ct2[match(proj_h3k4me2@cellColData$cell_id, d_cluster$cell_id)]
proj_h3k4me2 = addCellColData(ArchRProj = proj_h3k4me2, data = v_cluster_plus, name = "ct2", cells = proj_h3k4me2@cellColData$cell_id_ar, force = T)
saveArchRProject(ArchRProj = proj_h3k4me2, outputDirectory = "./tmp/H3K4me2_coCnT_final2/", load = F)

v_cluster_plus = d_cluster$ct2[match(proj_h3k4me3@cellColData$cell_id, d_cluster$cell_id)]
proj_h3k4me3 = addCellColData(ArchRProj = proj_h3k4me3, data = v_cluster_plus, name = "ct2", cells = proj_h3k4me3@cellColData$cell_id_ar, force = T)
saveArchRProject(ArchRProj = proj_h3k4me3, outputDirectory = "./tmp/H3K4me3_coCnT_final2/", load = F)

