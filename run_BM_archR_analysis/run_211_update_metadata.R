library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

library(ArchR)
color = ArchRPalettes$stallion


proj = loadArchRProject("./tmp/H3K27me3_coCnT_final3/")

d_meta = as.data.table(proj@cellColData)
d_umap = as.data.table(getEmbedding(proj, embedding = "UMAP_Harmony"))
names(d_umap) = c("UMAP_1", "UMAP_2")

d_meta = cbind(d_meta, d_umap)

write_tsv(d_umap, "./tmp/table_umap_h3k4me27_final3.tsv")
## TODO: Use the ArchR object to update the metadata and umap tables, instead of reading and writing tsv files.

d_meta = fread("../app_trajectory_shiny/data/cell_metadata.tsv")

write_tsv(d_meta, "./tmp/table_metadata_h3k27me3.tsv")

d_umap = fread("../app_trajectory_shiny/data/table_umap_h3k27me3_final3.tsv")

write_tsv(d_umap, "./tmp/table_umap_h3k4me27_final3.tsv")


d_plot = cbind(
	d_meta[, .(cell_id, ct3)],
	d_umap[, .(UMAP1, UMAP2)]
)

names(color) = unique(d_meta$ct3)[-1]

pdf("./figures/figure_umap_h3k27me3_final3.pdf", width = 6, height = 5)
ggplot(d_plot, aes(x = -UMAP1, y = UMAP2, color = ct3)) +
	geom_point(size = 0.5) +
	scale_color_manual(values = color) +
	theme_void() +
	theme(legend.position = "none")
dev.off()

# B cell clusters
b_cell_clusters = c(
    "HSC&MPP",
    "Pre-Pro-B cells",
    "Immature B cells",
    "Mature B cells1",
    "Mature B cells2",
    "Memory B cells",
    "Plasma cells"
)
# Monocyte clusters
m_cell_clusters = c(
    "HSC&MPP",
    "GMP",
    "Myelocyte/Classical Monocytes1",
    "Myelocyte/Classical Monocytes2",
    "Myelocyte/Classical Monocytes3"
)
# Erythroid clusters
e_cell_clusters = c(
    "HSC&MPP",
    "MEP",
    "Erythroid progenitors1",
    "Erythroid progenitors2",
    "Erythroid progenitors3"
)



## Add trajectories to the ArchR object
proj_h3k27 <- addTrajectory(
    ArchRProj = proj_h3k27, 
    name = "Bcell_Trajectory", 
    groupBy = "ct3",
    trajectory = b_cell_clusters, 
    embedding = "UMAP_Harmony", 
    force = TRUE
)

proj_h3k27 <- addTrajectory(
    ArchRProj = proj_h3k27, 
    name = "Monocyte_Trajectory", 
    groupBy = "ct3",
    trajectory = m_cell_clusters, 
    embedding = "UMAP_Harmony", 
    force = TRUE
)

proj_h3k27 <- addTrajectory(
    ArchRProj = proj_h3k27, 
    name = "Erythroid_Trajectory", 
    groupBy = "ct3",
    trajectory = e_cell_clusters, 
    embedding = "UMAP_Harmony", 
    force = TRUE
)

Cairo::CairoPDF("./tmp/figure_H3K27me3_Trajectory.pdf", width = 6, height = 5)
plotTrajectory(proj_h3k27, trajectory = "Bcell_Trajectory", colorBy = "cellColData", name = "Bcell_Trajectory", embedding = "UMAP_Harmony")
plotTrajectory(proj_h3k27, trajectory = "Monocyte_Trajectory", colorBy = "cellColData", name = "Monocyte_Trajectory", embedding = "UMAP_Harmony")
plotTrajectory(proj_h3k27, trajectory = "Erythroid_Trajectory", colorBy = "cellColData", name = "Erythroid_Trajectory", embedding = "UMAP_Harmony")
dev.off()


saveArchRProject(proj_h3k27, outputDirectory = "./tmp/H3K27me3_coCnT_final3/", load = T)

d_meta = as.data.table(proj_h3k27@cellColData)

write_tsv(d_meta, "./tmp/table_metadata_h3k27me3_final3_with_trajectories.tsv")

umap_h3k27 = as.data.table( getEmbedding( ArchRProj = proj_h3k27, embedding = "UMAP"))
names(umap_h3k27) = c("UMAP1", "UMAP2")
# umap_h3k27$cell_id = d_meta$cell_id

write_tsv(umap_h3k27, "./tmp/table_umap_h3k27me3_final3.tsv")


