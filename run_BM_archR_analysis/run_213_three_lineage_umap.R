## script: run_213_three_lineage_umap.R
## purpose: plot the UMAP with three lineage trajectories
## input: 
## - Updated metadata table with new cell type annotation and trajectory scores
## * ./tmp/table_metadata_h3k4me27_final3_HSPC.tsv
## output:
## - UMAP plot with the trajectory scores
## * ./figures/213_figure_umap_trajectory_h3k27me3_final3.pdf
##
library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

d_plot = fread("./tmp/table_metadata_h3k4me27_final3_HSPC.tsv")

table(d_meta$ct3)

d_plot = d_plot[, .(cell_id, Bcell_Trajectory, Monocyte_Trajectory, Erythroid_Trajectory, ct3, UMAP1, UMAP2)][
    ct3 %in% c(
	"HSC&MPP",
	"Pre-Pro-B cells",
	"Immature B cells",
	"Mature B cells1",
	"Mature B cells2",
	"Memory B cells",
	"Plasma cells",
	"GMP",
	"Myelocyte/Classical Monocytes1",
	"Myelocyte/Classical Monocytes2",
	"Myelocyte/Classical Monocytes3",
	"MEP",
	"Erythroid progenitors1",
	"Erythroid progenitors2",
	"Erythroid progenitors3"
)]

pal_B <- colorRampPalette(c("white", "navy"))(100)
pal_M <- colorRampPalette(c("white", "darkgreen"))(100)
pal_E <- colorRampPalette(c("white", "red"))(100)

d_plot2 <- rbindlist(list(
  d_plot[ct3 == "HSC&MPP", .(UMAP1, UMAP2, color = "grey90")],
  d_plot[ct3 != "HSC&MPP" & !is.na(Bcell_Trajectory),
         .(UMAP1, UMAP2, color = pal_B[pmin(100, pmax(1, round(Bcell_Trajectory)))])],
  d_plot[ct3 != "HSC&MPP" & !is.na(Monocyte_Trajectory),
         .(UMAP1, UMAP2, color = pal_M[pmin(100, pmax(1, round(Monocyte_Trajectory)))])],
  d_plot[ct3 != "HSC&MPP" & !is.na(Erythroid_Trajectory),
         .(UMAP1, UMAP2, color = pal_E[pmin(100, pmax(1, round(Erythroid_Trajectory)))])],
  d_plot[ct3 != "HSC&MPP" & is.na(Bcell_Trajectory) & is.na(Monocyte_Trajectory) & is.na(Erythroid_Trajectory),
		 .(UMAP1, UMAP2, color = "grey60")]
))

p = ggplot(d_plot2, aes(x = -UMAP1, y = UMAP2, color = color)) +
    geom_point(size = 0.5) +
    scale_color_identity() +
    theme_void() +
    theme(legend.position = "none")

ggsave("./figures/213_figure_umap_trajectory_h3k27me3_final3.pdf", p, width = 6, height = 5)


