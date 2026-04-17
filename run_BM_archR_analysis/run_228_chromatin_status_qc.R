library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

d_cnt = fread("./tmp/table_gene_marker_average_score_by_celltype.tsv.gz")

d_state_consensus = fread("./tmp/table_gene_chromatin_status_consensus_by_celltype_wide.tsv")
d_state_consensus_long = fread("./tmp/table_gene_chromatin_status_consensus_by_celltype.tsv")

httpgd::hgd(port = 4322)



## Bivalent gene number per cell type {{{
cell_type_order = c(
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


pdf("./figures/227_bivalent_gene_number_per_cell_type.pdf", width = 10, height = 5)
d_plot = d_state_consensus_long[consensus_count > 1]
d_plot = d_plot[, .N, by = .(cell_type, cat)]
d_plot$cell_type = factor(d_plot$cell_type, levels = cell_type_order)

ggplot(d_plot[cat != "Other"], aes(x = cell_type, y = N, fill = cat)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
    scale_fill_manual(values = c(
	    "Active" = "#21854e",
	    "Repressive" = "#bc3b2a",
	    "Bivalent" = "#e18727",
	    "Un" = "#1872b5"
	    )) +
    facet_wrap(~ cat, scales = "free_y", nrow = 1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(d_plot[cat != "Other"], aes(x = cell_type, y = N, color = cat)) +
    # geom_bar(stat = "identity", position = "dodge") +
    geom_point(size = 2) +
    geom_line(group = 1) +
    scale_color_manual(values = c(
	    "Active" = "#21854e",
	    "Repressive" = "#bc3b2a",
	    "Bivalent" = "#e18727",
	    "Un" = "#1872b5"
	    )) +
    theme_classic() +
    facet_wrap(~ cat, scales = "free_y", ncol = 1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


dev.off()
## }}}

## Scatter plot of H3K4me2 and H3K27me3 and H3K4me2-cooc {{{
pdf("./figures/222_scatter_plot_threshold.pdf")
d_plot_w = d_plot[cell_type %in% c("HSPC")] %>% 
    dcast(gene + cell_type + is_cpg ~ marker, value.var = "gene_score")

ggplot(d_plot_w, aes(x = H3K27me3 + 0.01, y = H3K4me2 + 0.1, color = H3K4me2_cooc + 0.1)) +
    geom_point(size = 0.7, alpha = 0.5) +
    geom_density2d(color = "black", alpha = 0.3) +
    scale_color_gradient(low = "grey", high = "orange", transform = "log10") +
    theme_classic() +
    geom_vline(xintercept = cutoff["H3K27me3"]) +
    geom_hline(yintercept = cutoff["H3K4me2"]) +
    scale_y_log10() + scale_x_log10() +
    ggtitle("HSPC: H3K4me2 vs H3K27me3 colored by H3K4me2 co-occurrence score")
dev.off()

## }}}

## Gene score per status {{{
state_v = d_state_consensus[, .(gene, HSPC)]
d_plot = merge(d_cnt[cell_type == "HSPC"], state_v, by = "gene")
d_plot = d_plot[HSPC %in% c("Active", "Bivalent", "Repressive", "Un")]
d_plot$HSPC = factor(d_plot$HSPC, levels = c("Active", "Repressive", "Bivalent", "Un"))

pdf("./figures/228_gene_score_by_chromatin_status.pdf")
comparisons = list(
	c("Active", "Bivalent"),
	c("Repressive", "Bivalent"),
	c("Bivalent", "Un")
)


ggplot(d_plot, aes(x = HSPC, y = gene_score, fill = HSPC)) +
    geom_boxplot(outliers = F, position="dodge") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~ marker, scales = "free") +
    stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = comparisons) +
    scale_fill_manual(values = c(
	    "Active" = "#21854e",
	    "Repressive" = "#bc3b2a",
	    "Bivalent" = "#e18727",
	    "Un" = "#1872b5"
	    )) +
    ggtitle("Gene score by chromatin status in HSPC") +
    xlab("Chromatin status") + ylab("Gene score")


ggplot(d_plot[marker %in% c("H3K4me2", "H3K27me3", "H3K4me2_cooc")], aes(x = HSPC, y = gene_score, fill = HSPC)) +
    geom_boxplot(outliers = F, position="dodge") +
    theme_classic() +
    facet_wrap(~ marker, scales = "free_y", nrow = 3) +
    scale_fill_manual(values = c(
	    "Active" = "#21854e",
	    "Repressive" = "#bc3b2a",
	    "Bivalent" = "#e18727",
	    "Un" = "#1872b5"
	    )) +
    ggtitle("Gene score by chromatin status in HSPC") +
    xlab("Chromatin status") + ylab("Gene score")

## do the statistical test for comparisons
stat_test = compare_means(gene_score ~ HSPC, data = d_plot[marker %in% c("H3K4me2", "H3K27me3", "H3K4me2_cooc")], group.by = "marker", method = "wilcox.test", p.adjust.method = "BH", comparisons = comparisons)
stat_test

dev.off()

## }}}
