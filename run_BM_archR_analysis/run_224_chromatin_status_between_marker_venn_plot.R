## script: run_224_chromatin_status_between_marker_venn_plot.R
## purpose: Compare the bivalent genes defined by different marks, and do gene set enrichment for the genes defined by multiple marks vs single mark
## input:
## - chromatin status by cell type for me1, me2, me3
## * ./tmp/table_gene_chromatin_status_wide_me1.tsv
## * ./tmp/table_gene_chromatin_status_wide_me2.tsv
## * ./tmp/table_gene_chromatin_status_wide_me3.tsv
## output:
## - Venn plot for bivalent genes defined by different marks in each cell type
## * ./figures/224_figure_chromatin_status_bivalent_venn.pdf
## - Gene set enrichment for genes defined by multiple marks vs single mark in HSPC
## * ./figures/224_figure_chromatin_status_bivalent_venn_hspc_enrichment.pdf
##
library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

d_me1 = fread("./tmp/table_gene_chromatin_status_wide_me1.tsv") %>% melt(id.vars = c("gene"), variable.name = "cell_type", value.name = "cat")
d_me2 = fread("./tmp/table_gene_chromatin_status_wide_me2.tsv") %>% melt(id.vars = c("gene"), variable.name = "cell_type", value.name = "cat")
d_me3 = fread("./tmp/table_gene_chromatin_status_wide_me3.tsv") %>% melt(id.vars = c("gene"), variable.name = "cell_type", value.name = "cat")

cell_type_v = d_me1$cell_type %>% unique()

library(eulerr)
library(grid)

pdf("./figures/224_figure_chromatin_status_bivalent_venn.pdf", width = 6, height = 6)
for (cell_type_x in cell_type_v) {
    print(cell_type_x)
    l_plot = list(
	me1_bi = d_me1[cell_type == cell_type_x & cat == "Bivalent", gene],
	me2_bi = d_me2[cell_type == cell_type_x & cat == "Bivalent", gene],
	me3_bi = d_me3[cell_type == cell_type_x & cat == "Bivalent", gene]
    )

    fit <- euler(l_plot)
    p <- plot(
	fit,
	fills = c("red", "green", "blue"),
	alpha = 0.5,
	labels = list(cex = 2),
	quantities = list(cex = 2),
	main = paste0("Bivalent genes in ", cell_type_x)
    )
    grid::grid.newpage()
    grid::grid.draw(p)
}
dev.off()

d_hspc = list(
    me1_bi_hspc = d_me1[cell_type == "HSPC" & cat == "Bivalent", .(gene, mark = "me1")],
    me2_bi_hspc = d_me2[cell_type == "HSPC" & cat == "Bivalent", .(gene, mark = "me2")],
    me3_bi_hspc = d_me3[cell_type == "HSPC" & cat == "Bivalent", .(gene, mark = "me3")]
) %>% rbindlist


v_gene_by_multi = d_hspc[, .N, by = gene][N > 1, gene]
v_gene_by_1 = d_hspc[, .N, by = gene][N == 1, gene]

## Do gene set enrichment
library(clusterProfiler)
library(org.Hs.eg.db)
eg = bitr(v_gene_by_multi, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ego_multi = enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = T)
ego_multi_df = as.data.table(ego_multi)
ego_multi_df[1:20, .(Description, GeneRatio, qvalue)]

eg = bitr(v_gene_by_1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ego_1 = enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = T)
ego_1_df = as.data.table(ego_1)
ego_1_df[1:20, .(Description, GeneRatio, qvalue)]

d_plot = rbind(
    ego_multi_df[, .(Description, FoldEnrichment, qvalue, group = "Multi-mark")],
    ego_1_df[, .(Description, FoldEnrichment = - FoldEnrichment, qvalue, group = "Single-mark")]
)

ggplot(d_plot, aes(x = FoldEnrichment, y = -log10(qvalue), color = group)) +
	geom_point() +
	theme_bw() +
	lims(x = c(-9, 9)) +
	## Add label for top 10 terms in each group
	geom_text(data = d_plot[qvalue < 0.05][order(qvalue)][, head(.SD, 20), by = group], aes(label = Description), position = position_jitter(width = 0.2, height = 0), size = 3) +
	theme_classic()
ggsave("./figures/224_figure_chromatin_status_bivalent_venn_hspc_enrichment.pdf", width = 6, height = 6)

## Negative control Mature B
d_me1$cell_type %>% unique
d_hspc = list(
    me1_bi_hspc = d_me1[cell_type == "Mature B cells1" & cat == "Active", .(gene, mark = "me1")],
    me2_bi_hspc = d_me2[cell_type == "Mature B cells1" & cat == "Active", .(gene, mark = "me2")],
    me3_bi_hspc = d_me3[cell_type == "Mature B cells1" & cat == "Active", .(gene, mark = "me3")]
) %>% rbindlist


v_gene_by_multi = d_hspc[, .N, by = gene][N > 1, gene]
v_gene_by_1 = d_hspc[, .N, by = gene][N == 1, gene]

## Do gene set enrichment
library(clusterProfiler)
library(org.Hs.eg.db)
eg = bitr(v_gene_by_multi, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ego_multi = enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = T)
ego_multi_df = as.data.table(ego_multi)
ego_multi_df[1:20, .(Description, GeneRatio, qvalue)]

eg = bitr(v_gene_by_1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ego_1 = enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = T)
ego_1_df = as.data.table(ego_1)
ego_1_df[1:20, .(Description, GeneRatio, qvalue)]



# d_me1 = fread("./tmp/table_chromatin_status_transitions_M_me1.tsv")
# d_me2 = fread("./tmp/table_chromatin_status_transitions_M_me2.tsv")
# d_me3 = fread("./tmp/table_chromatin_status_transitions_M_me3.tsv")
#
#
#
# library(eulerr)
# library(grid)
#
# pdf("./tmp/figure_chromatin_status_bivalent_venn_M.pdf", width = 6, height = 6)
# for (cell_type_x in cell_type_v) {
#     print(cell_type_x)
#     l_plot = list(
# 	me1_bi = d_me1[cell_type == cell_type_x & cat == "Bivalent", gene],
# 	me2_bi = d_me2[cell_type == cell_type_x & cat == "Bivalent", gene],
# 	me3_bi = d_me3[cell_type == cell_type_x & cat == "Bivalent", gene]
#     )
#
#
#     fit <- euler(l_plot, shape = "ellipse")
#     p <- plot(
# 	fit,
# 	fills = c("red", "green", "blue"),
# 	alpha = 0.5,
# 	labels = list(cex = 2),
# 	quantities = list(cex = 2),
# 	main = paste0("Bivalent genes in ", cell_type_x)
#     )
#     grid::grid.newpage()
#     grid::grid.draw(p)
# }
# dev.off()
#
# d_me1 = fread("./tmp/table_chromatin_status_transitions_E_me1.tsv")
# d_me2 = fread("./tmp/table_chromatin_status_transitions_E_me2.tsv")
# d_me3 = fread("./tmp/table_chromatin_status_transitions_E_me3.tsv")
#
# cell_type_v = d_me1$cell_type %>% unique()
#
#
# library(eulerr)
# library(grid)
#
# pdf("./tmp/figure_chromatin_status_bivalent_venn_E.pdf", width = 6, height = 6)
# for (cell_type_x in cell_type_v) {
#     print(cell_type_x)
#     l_plot = list(
# 	me1_bi = d_me1[cell_type == cell_type_x & cat == "Bivalent", gene],
# 	me2_bi = d_me2[cell_type == cell_type_x & cat == "Bivalent", gene],
# 	me3_bi = d_me3[cell_type == cell_type_x & cat == "Bivalent", gene]
#     )
#
#
#     fit <- euler(l_plot, shape = "ellipse")
#     p <- plot(
# 	fit,
# 	fills = c("red", "green", "blue"),
# 	alpha = 0.5,
# 	labels = list(cex = 2),
# 	quantities = list(cex = 2),
# 	main = paste0("Bivalent genes in ", cell_type_x)
#     )
#     grid::grid.newpage()
#     grid::grid.draw(p)
# }
# dev.off()

