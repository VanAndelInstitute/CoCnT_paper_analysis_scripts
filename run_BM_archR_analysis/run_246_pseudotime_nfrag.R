library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

httpgd::hgd(port = 4323)

d_meta_k27 <- fread("./tmp/table_metadata_h3k27me3_final3_with_trajectories.tsv")
d_meta_k4me1 <- fread("./tmp/table_metadata_h3k4me1_final3.tsv")
d_meta_k4me2 <- fread("./tmp/table_metadata_h3k4me2_final3.tsv")
d_meta_k4me3 <- fread("./tmp/table_metadata_h3k4me3_final3.tsv")
d_meta_k4me1_cooc <- fread("./tmp/table_metadata_h3k4me1_final3_with_cooc.tsv")
d_meta_k4me2_cooc <- fread("./tmp/table_metadata_h3k4me2_final3_with_cooc.tsv")
d_meta_k4me3_cooc <- fread("./tmp/table_metadata_h3k4me3_final3_with_cooc.tsv")

d_meta_l <- list(
    k27 = d_meta_k27,
    k4me1 = d_meta_k4me1,
    k4me2 = d_meta_k4me2,
    k4me3 = d_meta_k4me3,
    k4me1_cooc = d_meta_k4me1_cooc,
    k4me2_cooc = d_meta_k4me2_cooc,
    k4me3_cooc = d_meta_k4me3_cooc
)


d_meta <- lapply(1:length(d_meta_l), function(i) {
    x <- d_meta_l[[i]][, .(cell_id, nFrags)]
    x$marker <- names(d_meta_l)[i]
    x
}) %>% rbindlist()

d_trajectory <- d_meta_k27[, .(cell_id, b_trajectory = Bcell_Trajectory, e_trajectory = Erythroid_Trajectory, m_trajectory = Monocyte_Trajectory, ct3)]

# m_list <- readRDS("./tmp/list_matrix_imputed_gene_score.rds")

d_meta <- merge(d_meta, d_trajectory, by = "cell_id", all.x = TRUE)

ct_by_diff <- c(
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
)


d_meta$ct3 %<>% factor(levels = ct_by_diff)


## nFrags vs trajectory {{{

## Three lineages
lineages_v <- c("e_trajectory", "b_trajectory", "m_trajectory")

## times 7 markers
markers_v <- c("k27", "k4me1", "k4me2", "k4me3", "k4me1_cooc", "k4me2_cooc", "k4me3_cooc")

pdf("./figures/246_nfrag_trajectory_trajectory_joint_lineage.pdf", width = 6, height = 4)
d_plot = d_meta[marker == "k27", .(cell_id, nFrags, b_trajectory, e_trajectory, m_trajectory, ct3)]
d_plot_long <- melt(d_plot, id.vars = c("cell_id", "nFrags", "ct3"), measure.vars = c("b_trajectory", "e_trajectory", "m_trajectory"), variable.name = "lineage", value.name = "trajectory")
ggplot(d_plot_long, aes(x = trajectory, y = nFrags, color = lineage)) +
	geom_point(alpha = 0.1) +
	geom_smooth(method = "loess", se = FALSE, span = 0.2) +
	theme_classic() +
	scale_y_log10() +
	labs(title = "k27 - Trajectory", x = "Trajectory")
	# facet_wrap(~lineage, scales = "free_x")
dev.off()


pdf("./figures/242_nfrag_trajectory.pdf", width = 6, height = 4)

for (i in 1:length(lineages_v)) {
    for (j in 1:length(markers_v)) {
        x_lineages <- lineages_v[i]
        x_markers <- markers_v[j]

        print(paste0("Plotting ", x_markers, " - ", x_lineages))

        d_plot <- d_meta[marker == x_markers, .(cell_id, nFrags, trajectory = get(x_lineages), ct3)] %>% na.omit()

        p1 <- ggplot(d_plot, aes(x = trajectory, y = nFrags, color = ct3)) +
            geom_point() +
            theme_classic() +
            scale_y_log10() +
            labs(title = paste0(x_markers, " - ", x_lineages), x = x_lineages)

        p2 <- ggplot(d_plot, aes(x = ct3, y = nFrags, color = trajectory)) +
            geom_violin() +
            theme_classic() +
            scale_y_log10() +
            labs(title = paste0(x_markers, " - ", x_lineages), x = "Cell Type")

        print(p1)
        print(p2)
    }
}
dev.off()

## }}}
