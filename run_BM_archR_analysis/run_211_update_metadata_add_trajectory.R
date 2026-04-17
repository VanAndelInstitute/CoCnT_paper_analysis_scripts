library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

library(ArchR)
color = ArchRPalettes$stallion

## Update Metadata and UMAP
proj_h3k27 = loadArchRProject("./tmp/H3K27me3_coCnT_final3/")

## Add trajectories to the ArchR object
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

## Output the trajectory plots
Cairo::CairoPDF("./figures/211_H3K27me3_Trajectory.pdf", width = 6, height = 5)
plotTrajectory(proj_h3k27, trajectory = "Bcell_Trajectory", colorBy = "cellColData", name = "Bcell_Trajectory", embedding = "UMAP_Harmony")
plotTrajectory(proj_h3k27, trajectory = "Monocyte_Trajectory", colorBy = "cellColData", name = "Monocyte_Trajectory", embedding = "UMAP_Harmony")
plotTrajectory(proj_h3k27, trajectory = "Erythroid_Trajectory", colorBy = "cellColData", name = "Erythroid_Trajectory", embedding = "UMAP_Harmony")
dev.off()

## Update the ArchR project
saveArchRProject(proj_h3k27, outputDirectory = "./tmp/H3K27me3_coCnT_final3/", load = T)

## Output the metadata and umap tables

d_meta = as.data.table(proj_h3k27@cellColData)

write_tsv(d_meta, "./tmp/table_metadata_h3k4me27_final3_HSPC_with_trajectories.tsv")

umap_h3k27 = as.data.table( getEmbedding( ArchRProj = proj_h3k27, embedding = "UMAP"))
names(umap_h3k27) = c("UMAP1", "UMAP2")

write_tsv(umap_h3k27, "./tmp/table_umap_h3k4me27_final3_HSPC.tsv")


