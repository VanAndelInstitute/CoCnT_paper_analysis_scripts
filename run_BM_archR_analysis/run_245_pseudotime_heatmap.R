## script: run_245_pseudotime_heatmap.R
## purpose: Generate heatmaps of gene activity along pseudotime trajectories for specific gene sets categorized
## by chromatin state transitions, and compare activation times between different transition types.
## input:
## - Imputed gene score matrices for different histone marks
## * ./tmp/list_matrix_imputed_gene_score.rds
## - Metadata with pseudotime trajectories for each cell
## * ./tmp/table_metadata_h3k4me27_final3_HSPC.tsv
## - Gene categorization based on chromatin state transitions
## * ./tmp/table_gene_chromatin_status_trajectory_consensus_remove_dup_long.tsv
## output:
## - Heatmaps of gene activity along pseudotime trajectories for genes categorized by chromatin state
## transitions (Un->Active and Bivalent->Active) for each lineage
## * ./figures/245_pseudotime_heatmap_activation_order_50_me2.pdf
## - Comparison of activation times between Un->Active and Bivalent->Active gene sets for each lineage
## * ./figures/245_activation_time_comparison_50_me2.pdf
## 
library(data.table)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(readr)
library(stringr)
library(pheatmap)

m_list <- readRDS("./tmp/list_matrix_imputed_gene_score.rds")
d_meta <- fread("./tmp/table_metadata_h3k4me27_final3_HSPC.tsv")
gene_category <- fread("./tmp/table_gene_chromatin_status_trajectory_consensus_remove_dup_long.tsv")

d <- fread("./tmp/table_gene_chromatin_status_consensus_by_celltype.tsv")


marker_colors <- c(
    H3K27me3 = "#3B8FC4",
    H3K4me1 = "#FFD000",
    H3K4me2 = "#F39C12",
    H3K4me3 = "#E74C3C",
    H3K4me1_cooc = "#56B870",
    H3K4me2_cooc = "#D16BA5",
    H3K4me3_cooc = "#9B59B6"
)

trajectory_map <- list(
    "B" = list(name = "Bcell_Trajectory"),
    "M" = list(name = "Monocyte_Trajectory"),
    "E" = list(name = "Erythroid_Trajectory")
)


available_markers <- intersect(names(m_list), names(marker_colors))
available_genes <- rownames(m_list[[1]])

## Function to process the heatmap data and plot {{{
process_heatmap <- function(mat) {
    # Smoothing with loess
    m_fit <- t(apply(mat, 1, function(x) {
        fit <- loess(x ~ seq_along(x))
        predict(fit)
    }))
    m_fit[!is.finite(m_fit)] <- NA
    keep_finite <- apply(m_fit, 1, function(x) any(is.finite(x)))
    m_fit <- m_fit[keep_finite, , drop = FALSE]
    m_fit[!is.finite(m_fit)] <- 0
    m_fit[m_fit < 0] <- 0

    m_fit2 <- m_fit

    # Scale each gene to [0,1]
    m_scaled <- t(apply(m_fit2, 1, function(x) {
        rng <- max(x) - min(x)
        if (rng == 0) rep(0, length(x)) else (x - min(x)) / rng
    }))

    # Activation time: first time >= 0.95
    activation_time <- apply(m_scaled, 1, function(x) {
        idx <- which(x >= 0.50)
        if (length(idx) == 0) Inf else idx[1]
    })

    activation_time <- activation_time / ncol(m_scaled) * 100

    # Keep only genes with a valid activation event
    keep_act <- is.finite(activation_time) & activation_time > 1
    m_act <- m_scaled[keep_act, , drop = FALSE]
    activation_time_act <- activation_time[keep_act]

    # Optional tie-breakers
    peak_time <- apply(m_act, 1, which.max)
    activation_strength <- apply(m_act, 1, function(x) max(x) - mean(x))

    ord <- order(activation_time_act, peak_time, -activation_strength)
    m <- m_act[ord, , drop = FALSE]

    # png(filename = "./figures/245_pseudotime_heatmap_example.png", width = 800, height = 1200)
    pheatmap(
        m,
        show_rownames = FALSE,
        show_colnames = FALSE,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("white", "red"))(100)
    )
    # dev.off()
    return(activation_time_act)
}

## }}}

x_marker <- "H3K4me2"

pdf("./figures/245_pseudotime_heatmap_activation_order_50_me2.pdf", width = 8, height = 12)

lapply(1:length(trajectory_map), function(i) {
    cat("Processing lineage:", names(trajectory_map)[i], "\n")
    lineage_x <- names(trajectory_map)[i]

    group_x <- paste0(lineage_x, "_consensus")
    x_trajectory <- trajectory_map[[lineage_x]]$name

    m_x <- m_list[[x_marker]]

    order_cell_d <- d_meta[order(get(x_trajectory)), .(cell_id, get(x_trajectory))] %>% na.omit()
    order_cell_d <- order_cell_d[cell_id %in% colnames(m_x)]


    m_x_cell_order <- m_x[, order_cell_d$cell_id]

    m_x_cell_order
    dim(m_x_cell_order)

    gene_v_u2a <- gene_category[trans %in% c("Un->Active", "Un->Active->Un") & group == group_x, gene]
    gene_v_b2a <- gene_category[trans %in% c("Bivalent->Active", "Bivalent->Active->Bivalent") & group == group_x, gene]

    m_x_cell_order_u2a <- m_x_cell_order[gene_v_u2a, ]
    m_x_cell_order_b2a <- m_x_cell_order[gene_v_b2a, ]

    dim(m_x_cell_order_u2a)
    dim(m_x_cell_order_b2a)

    ## Process heatmaps for each
    activation_time_u2a <- process_heatmap(mat = m_x_cell_order_u2a)
    activation_time_b2a <- process_heatmap(mat = m_x_cell_order_b2a)
    rbind(
        data.table(gene = names(activation_time_u2a), activation_time = activation_time_u2a, type = "U2A", lineage = lineage_x),
        data.table(gene = names(activation_time_b2a), activation_time = activation_time_b2a, type = "B2A", lineage = lineage_x)
    )
}) -> activation_time_list

dev.off()


activation_time_list %>% rbindlist() -> activation_time_df_all

#########################
#  Compare U2A and B2A  #
#########################

## normalize the pseudotime index to [0,100] for better interpretability by lineage
# activation_time_df_all[, max_time := max(activation_time[is.finite(activation_time)], na.rm = TRUE), by = lineage]
# activation_time_df_all[, activation_time := activation_time / max_time * 100]

pdf("./figures/245_activation_time_comparison_50_me2.pdf", width = 6, height = 4)
ggplot(activation_time_df_all[activation_time > 1], aes(x = type, y = activation_time, fill = type)) +
    geom_violin(outlier.shape = NA) +
    # geom_boxplot(width = 0.1, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
    facet_wrap(~lineage) +
    labs(x = "Transition Type", y = "Activation Time (normalized)") +
    # stat_compare_means(method = "wilcox.test", label = "p.format") +
    theme_classic() +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("U2A" = "#56B870", "B2A" = "#D16BA5"))
dev.off()

compare_means(activation_time ~ type, data = activation_time_df_all[activation_time > 1], method = "wilcox.test", group.by = "lineage", p.adjust.method = "BH")



