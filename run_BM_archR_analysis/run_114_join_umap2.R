##################################################
## Script: run_n_join_umap2.R
## Purpose: Join UMAP embeddings from multiple ArchR projects and visualize them together
## INPUT:
## - ArchR projects (directories):
## * ./tmp/H3K27me3_coCnT_final2/
## * ./tmp/H3K4me1_coCnT_final2/
## * ./tmp/H3K4me2_coCnT_final2/
## * ./tmp/H3K4me3_coCnT_final2/
## OUTPUT:
## - Combined UMAP embeddings and cluster annotations for all projects
## * ./figures/114_join_umap2.pdf
## - Updated ArchR projects with combined cluster annotations:
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
library(ggsci)

library(ArchR)
addArchRGenome("hg38")
color <- ArchRPalettes$stallion
names(color) <- paste0("C", (1:(length(color))))

proj_k27 <- loadArchRProject("./tmp/H3K27me3_coCnT_final2/")
proj_h3k4me1 <- loadArchRProject("./tmp/H3K4me1_coCnT_final2/")
proj_h3k4me2 <- loadArchRProject("./tmp/H3K4me2_coCnT_final2/")
proj_h3k4me3 <- loadArchRProject("./tmp/H3K4me3_coCnT_final2/")


d_proj_k27 <- as.data.table(proj_k27@cellColData)
d_proj_h3kme1 <- as.data.table(proj_h3k4me1@cellColData)
d_proj_h3kme2 <- as.data.table(proj_h3k4me2@cellColData)
d_proj_h3kme3 <- as.data.table(proj_h3k4me3@cellColData)

d_cluster <- d_proj_k27[, .(cell_id, cluster_k27 = Clusters)]
d_cluster <- merge(d_cluster, d_proj_h3kme1[, .(cell_id, cluster_h3kme1 = Clusters)], by = "cell_id", all = T)
d_cluster <- merge(d_cluster, d_proj_h3kme2[, .(cell_id, cluster_h3kme2 = Clusters)], by = "cell_id", all = T)
d_cluster <- merge(d_cluster, d_proj_h3kme3[, .(cell_id, cluster_h3kme3 = Clusters)], by = "cell_id", all = T)

umap_k27 <- as.data.table(getEmbedding(ArchRProj = proj_k27, embedding = "UMAP_Harmony"))
umap_h3kme1 <- as.data.table(getEmbedding(ArchRProj = proj_h3k4me1, embedding = "UMAP_Harmony"))
umap_h3kme2 <- as.data.table(getEmbedding(ArchRProj = proj_h3k4me2, embedding = "UMAP_Harmony"))
umap_h3kme3 <- as.data.table(getEmbedding(ArchRProj = proj_h3k4me3, embedding = "UMAP_Harmony"))

names(umap_k27) <- c("UMAP1", "UMAP2")
names(umap_h3kme1) <- c("UMAP1", "UMAP2")
names(umap_h3kme2) <- c("UMAP1", "UMAP2")
names(umap_h3kme3) <- c("UMAP1", "UMAP2")

scale2one <- function(x) {
    x / (max(x) - min(x)) / 1.5
}

## Plot triple UMAP

umap_k27[, UMAP1 := scale2one(UMAP1)]
umap_k27[, UMAP2 := scale2one(UMAP2)]
umap_h3kme1[, UMAP1 := scale2one(UMAP1) - 1]
umap_h3kme1[, UMAP2 := scale2one(UMAP2) - 1]
umap_h3kme2[, UMAP1 := scale2one(UMAP1) + 0]
umap_h3kme2[, UMAP2 := scale2one(UMAP2) + 1]
umap_h3kme3[, UMAP1 := scale2one(UMAP1) + 1]
umap_h3kme3[, UMAP2 := scale2one(UMAP2) - 1]

umap_k27[, cell_id := d_proj_k27$cell_id]
umap_h3kme1[, cell_id := d_proj_h3kme1$cell_id]
umap_h3kme2[, cell_id := d_proj_h3kme2$cell_id]
umap_h3kme3[, cell_id := d_proj_h3kme3$cell_id]

umap_k27[, mark := "H3K27me3"]
umap_h3kme1[, mark := "H3K4me1"]
umap_h3kme2[, mark := "H3K4me2"]
umap_h3kme3[, mark := "H3K4me3"]


d_plot <- rbind(umap_k27, umap_h3kme1, umap_h3kme2, umap_h3kme3)
d_plot <- merge(d_plot, d_cluster, by = "cell_id")

pdf("./figures/114_join_umap2.pdf", width = 12, height = 10)
ggplot(d_plot[, .(UMAP1, UMAP2, cluster_k27)] %>% na.omit(), aes(x = UMAP1, y = UMAP2, color = cluster_k27)) +
    # geom_line(aes(group = cell_id), alpha = 0.03, size = 0.05) +
    geom_point(size = 0.3, alpha = 0.8) +
    scale_color_manual(values = color) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    theme_classic() +
    theme(
        legend.position = "right",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()
    )


## Plot dual UMAP
# d_plot$cluster_h3kme1 %<>% factor(levels = c(paste0("C", 1:length(unique(d_plot$cluster_h3kme1)))), NA)
ggplot(d_plot[rev(order(cluster_h3kme1))], aes(x = UMAP1, y = UMAP2, color = cluster_h3kme1)) +
    # geom_line(aes(group = cell_id), alpha = 0.03, size = 0.05) +
    geom_point(size = 0.05, alpha = 0.2) +
    scale_color_manual(values = color) +
    theme_classic() +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    theme(
        legend.position = "right",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()
    )


ggplot(d_plot[rev(order(cluster_h3kme2))], aes(x = UMAP1, y = UMAP2, color = cluster_h3kme2)) +
    # geom_line(aes(group = cell_id), alpha = 0.03, size = 0.05) +
    geom_point(size = 0.05, alpha = 0.3) +
    scale_color_manual(values = color) +
    theme_classic() +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    theme(
        legend.position = "right",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()
    )

ggplot(d_plot[rev(order(cluster_h3kme3))], aes(x = UMAP1, y = UMAP2, color = cluster_h3kme3)) +
    # geom_line(aes(group = cell_id), alpha = 0.03, size = 0.05) +
    geom_point(size = 0.05, alpha = 0.3) +
    scale_color_manual(values = color) +
    theme_classic() +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    theme(
        legend.position = "right",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()
    )
dev.off()

## Update achR proj
v_cluster_plus <- d_cluster$cluster_k27[match(proj_k27@cellColData$cell_id, d_cluster$cell_id)]
proj_k27 <- addCellColData(ArchRProj = proj_k27, data = v_cluster_plus, name = "Clusters_plus", cells = proj_k27@cellColData$cell_id_ar, force = T)
saveArchRProject(ArchRProj = proj_k27, outputDirectory = "./tmp/H3K27me3_coCnT_final2/", load = F)

v_cluster_plus <- d_cluster$cluster_k27[match(proj_h3k4me1@cellColData$cell_id, d_cluster$cell_id)]
proj_h3k4me1 <- addCellColData(ArchRProj = proj_h3k4me1, data = v_cluster_plus, name = "Clusters_plus", cells = proj_h3k4me1@cellColData$cell_id_ar, force = T)
saveArchRProject(ArchRProj = proj_h3k4me1, outputDirectory = "./tmp/H3K4me1_coCnT_final2/", load = F)

v_cluster_plus <- d_cluster$cluster_k27[match(proj_h3k4me2@cellColData$cell_id, d_cluster$cell_id)]
proj_h3k4me2 <- addCellColData(ArchRProj = proj_h3k4me2, data = v_cluster_plus, name = "Clusters_plus", cells = proj_h3k4me2@cellColData$cell_id_ar, force = T)
saveArchRProject(ArchRProj = proj_h3k4me2, outputDirectory = "./tmp/H3K4me2_coCnT_final2/", load = F)

v_cluster_plus <- d_cluster$cluster_k27[match(proj_h3k4me3@cellColData$cell_id, d_cluster$cell_id)]
proj_h3k4me3 <- addCellColData(ArchRProj = proj_h3k4me3, data = v_cluster_plus, name = "Clusters_plus", cells = proj_h3k4me3@cellColData$cell_id_ar, force = T)
saveArchRProject(ArchRProj = proj_h3k4me3, outputDirectory = "./tmp/H3K4me3_coCnT_final2/", load = F)

