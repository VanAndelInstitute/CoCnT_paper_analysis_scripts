# =============================================================================
# Script: run_222_chromatin_status_dynamics_v2.R
#
# Input:
#   - table_gene_marker_average_score_by_celltype.tsv.gz: marker gene scores by cell type (compressed).
#
# Output:
#   - A PDF file ('./figures/223_figure_chromatin_status_dynamics_v2.pdf') visualizing the chromatin status dynamics for each cell type.
# =============================================================================

library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(ggsci)


# d_cnt = fread("./tmp/table_gene_marker_average_score_by_celltype_ct2_merge.tsv.gz")
d_cnt = fread("./tmp/table_gene_marker_average_score_by_celltype.tsv.gz")
d_plot = d_cnt

## Get gene list with CpG islands
d_cpg_gene = readLines("./tmp/table_gene_cpg_island_promoter_hg19.txt")

d_plot[, is_cpg := gene %in% d_cpg_gene]

## Categorize the genes into 
## 1. Bivalent
## 2. Active
## 3. Repressive
## 4. Others


plot_chromatin_status <- function(file_path, plot_title = NULL) {
    d_wide <- fread(file_path)
    d_wide$cat <- factor(d_wide$cat, levels = c("Active", "Repressive", "Bivalent", "Un", "Other"))

    # Remove "Other" category genes for HSPC
    d_plot1 <- d_wide[cell_type == "HSPC", .N, by = cat][cat != "Other"][order(cat)]
    d_plot1[, fraction := N / sum(N)]

    v_cat <- unique(d_plot1$cat)

    d_plot_n <- lapply(seq_along(v_cat), function(i) {
        cat_i <- v_cat[i]
        v_gene <- unique(d_wide[cell_type == "HSPC" & cat == cat_i]$gene)
        d_plot_i <- d_wide[gene %in% v_gene, .N, by = .(cell_type, cat)]
        d_plot_i <- d_plot_i[cat != "Other"]
        d_plot_i[, fraction := N / sum(N), by = cell_type]
        d_plot_i[, belong_cat := cat_i]
        d_plot_i
    }) %>% rbindlist()

    colors_cat <- c(
        "Active" = "#21854e",
        "Repressive" = "#bc3b2a",
        "Bivalent" = "#e18727",
        "Un" = "#1872b5"
    )

    ct_by_diff <- c(
        "HSPC",
        "Pre-Pro-B cells",
        "Immature B cells",
        "Mature B cells1",
        "Mature B cells2",
        "Memory B cells",
        "Plasma cells",
        "Myelocyte/Classical Monocytes1",
        "Myelocyte/Classical Monocytes2",
        "Myelocyte/Classical Monocytes3",
        "Erythroid progenitors1",
        "Erythroid progenitors2",
        "Erythroid progenitors3"
    )

    plot(c(1, 15), y = c(0, 1), type = "n", xlab = "Cell Type", ylab = "Fraction of genes", xaxt = "n", main = plot_title)

    # Plot HSPC bar
    x <- 1
    y1 <- 0
    for (i in seq_len(nrow(d_plot1))) {
        y2 <- y1 + d_plot1[i]$fraction
        rect(xleft = x - 0.4, ybottom = y1, xright = x + 0.4, ytop = y2, col = colors_cat[d_plot1[i]$cat], border = NA)
        y1 <- y2
    }

    # Plot each cell type bar
    x <- x + 1
    for (j in 2:length(ct_by_diff)) {
        cell_type_now <- ct_by_diff[j]
        d_plot_i <- d_plot_n[cell_type == cell_type_now]
        d_plot_i <- d_plot_i[order(belong_cat, cat)]
        y1 <- 0
        for (i in seq_len(nrow(d_plot_i))) {
            y2 <- y1 + d_plot_i[i]$fraction * (d_plot1[cat == d_plot_i[i]$belong_cat]$fraction)
            rect(xleft = x - 0.4, ybottom = y1, xright = x + 0.4, ytop = y2, col = colors_cat[d_plot_i[i]$cat], border = NA)
            y1 <- y2
        }
        x <- x + 1
    }
}

# Plot for each chromatin status file
Cairo::CairoPDF("./figures/223_figure_chromatin_status_dynamics_v2.pdf", width = 10, height = 6)
plot_chromatin_status("./tmp/table_gene_chromatin_status_by_celltype_me1.tsv", plot_title = "ME1")
plot_chromatin_status("./tmp/table_gene_chromatin_status_by_celltype_me2.tsv", plot_title = "ME2")
plot_chromatin_status("./tmp/table_gene_chromatin_status_by_celltype_me3.tsv", plot_title = "ME3")
plot_chromatin_status("./tmp/table_gene_chromatin_status_consensus_by_celltype.tsv", plot_title = "Consensus")
dev.off()
