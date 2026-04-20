## script: run_227_x2active.R
## purpose: Compare the chromatin status in HSPC with the active state in other cell type 
## input:
## - Chromatin status consensus table
## * ./tmp/table_gene_chromatin_status_consensus_by_celltype_wide.tsv
## - HSPC repressive score
## * ./tmp/table_gene_marker_average_score_by_celltype.tsv.gz
## output:
## - Proportion of genes with at least 1 active state in other cell type for each
## chromatin status category in HSPC
## * ./figures/227_figure_chromatin_status_x2active_proportion_active_once.pdf
## * ./figures/227_figure_chromatin_status_x2active_proportion_active.pdf
## * ./figures/227_figure_chromatin_status_x2active_proportion_n_celltype_active.pdf
## - HSPC repressive score for genes with Repressive->Active transition and Repressive without Active transition
## * ./figures/227_figure_chromatin_status_x2active_repressive_score.pdf
## - Proportion of genes with at least 1 repressive state in other cell type for each chromatin status category in HSPC
## * ./figures/227_figure_chromatin_status_x2repressive_proportion_active.pdf
##
library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(ggsci)

## Proportion of each category endup with at least once active state {{{
d = fread("./tmp/table_gene_chromatin_status_consensus_by_celltype_wide.tsv")
d_status_m = d[, c(-1, -5), with = T]
active_count = apply(d_status_m, 1, function(x) (sum(x == "Active", na.rm = T)))

lapply(1:11, function(i) {
    res = d[, .(hspc_cat = HSPC, active_n = active_count)] %>%
	.[, .(n = sum(active_n >= i) / .N), by = hspc_cat]
    res$active_n = i
    res
}) %>% rbindlist() -> d_active_count

d_plot = d_active_count[active_n == 1, .(hspc_cat, n)][hspc_cat %in% c("Active", "Repressive", "Un", "Bivalent")]
d_plot$hspc_cat = factor(d_plot$hspc_cat, levels = c("Active", "Repressive", "Bivalent", "Un"))
p = ggplot(d_plot, aes(x = hspc_cat, y = n, fill = hspc_cat)) +
	geom_bar(stat = "identity", position = "dodge") +
	theme_classic() +
	scale_fill_manual(values = c(
	    "Active" = "#21854e",
	    "Repressive" = "#bc3b2a",
	    "Bivalent" = "#e18727",
	    "Un" = "#1872b5"
	)) +
	labs(x = "Chromatin status in HSPC", y = "Proportion of genes with at least 1 active state in other cell type") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("./figures/227_figure_chromatin_status_x2active_proportion_active_once.pdf", p, width = 3, height = 3)


p = ggplot(d_active_count[hspc_cat %in% c("Active", "Repressive", "Un", "Bivalent")], aes(x = active_n, y = n, color = hspc_cat)) +
    geom_point() +
    geom_line() +
    theme_classic() +
    scale_color_manual(values = c(
	    "Active" = "#21854e",
	    "Repressive" = "#bc3b2a",
	    "Bivalent" = "#e18727",
	    "Un" = "#1872b5"
    )) +
    labs(x = "More than n cell type with active state", y = "Proportion of genes") +
    ## add x == 1 axis tick
    scale_x_continuous(breaks = 1:11)


ggsave("./figures/227_figure_chromatin_status_x2active_proportion_active.pdf", p , width = 5, height = 3)



lapply(1:11, function(i) {
    res = d[, .(hspc_cat = HSPC, active_n = active_count)] %>%
	.[, .(n = sum(active_n == i) / .N), by = hspc_cat]
    res$active_n = i
    res
}) %>% rbindlist() -> d_active_count


p = ggplot(d_active_count[hspc_cat %in% c("Active", "Repressive", "Un", "Bivalent")], aes(x = active_n, y = n, color = hspc_cat)) +
    geom_point() +
    geom_line() +
    theme_classic() +
    scale_color_manual(values = c(
	    "Active" = "#21854e",
	    "Repressive" = "#bc3b2a",
	    "Bivalent" = "#e18727",
	    "Un" = "#1872b5"
    )) +
    labs(x = "active state cell type number", y = "Proportion of genes") +
    ## add x == 1 axis tick
    scale_x_continuous(breaks = 1:11)


ggsave("./figures/227_figure_chromatin_status_x2active_proportion_n_celltype_active.pdf", p, width = 5, height = 3)



## }}}

## Does the gene with Repressive->Active has lower Repressive gene score in HSPC {{{

# the gene has repressive genes that has ever been active once

d = fread("./tmp/table_gene_chromatin_status_consensus_by_celltype_wide.tsv")
d_status_m = d[, c(-1, -5), with = T]
active_count = apply(d_status_m, 1, function(x) (sum(x == "Active", na.rm = T)))

d_plot = d[, .(gene = d$gene, hspc_cat = HSPC, active_n = active_count)]

# get Hspc repressive score
d_score = fread("./tmp/table_gene_marker_average_score_by_celltype.tsv.gz")

d_plot = merge(
    d_plot,
    d_score[cell_type == "HSPC" & marker == "H3K27me3", .(gene, gene_score)],
    by = "gene",
    all.x = T
)

pdf("./figures/227_figure_chromatin_status_x2active_repressive_score.pdf", width = 5, height = 4)
ggplot(d_plot[hspc_cat %in% c("Bivalent", "Un", "Repressive")]) + aes(x = hspc_cat, y = gene_score + 0.1, fill = active_n >= 1) +
    geom_boxplot(outliers = F) +
    scale_fill_manual(values = c("FALSE" = "#e18727", "TRUE" = "#21854e"), name = "Ever active in other cell type") +
    geom_hline(yintercept = 0.205, linetype = "dashed") +
    theme_classic() +
	labs(x = "Chromatin status in HSPC", y = "H3K27me3 score in HSPC") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
	scale_y_log10()

ggplot(d_plot[hspc_cat %in% c("Bivalent", "Un", "Repressive")]) + aes(x = hspc_cat, y = gene_score + 0.1, fill = active_n >= 1) +
    geom_boxplot(outliers = F) +
    scale_fill_manual(values = c("FALSE" = "#e18727", "TRUE" = "#21854e"), name = "Ever active in other cell type") +
    geom_hline(yintercept = 0.205, linetype = "dashed") +
    theme_classic() +
	labs(x = "Chromatin status in HSPC", y = "H3K27me3 score in HSPC") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

dev.off()

## }}}

## Will bivalent prevent repression {{{

d = fread("./tmp/table_gene_chromatin_status_consensus_by_celltype_wide.tsv")
d_status_m = d[, c(-1, -5), with = T]
active_count = apply(d_status_m, 1, function(x) (sum(x == "Repressive", na.rm = T)))

lapply(1:11, function(i) {
    res = d[, .(hspc_cat = HSPC, active_n = active_count)] %>%
	.[, .(n = sum(active_n >= i) / .N), by = hspc_cat]
    res$active_n = i
    res
}) %>% rbindlist() -> d_active_count

# ggplot(d_active_count[hspc_cat %in% c("Active", "Repressive", "Un", "Bivalent")], aes(x = active_n, y = n, color = hspc_cat)) +
#     geom_point() +
#     geom_line() +
#     theme_classic() +
#     scale_color_manual(values = c(
# 	    "Active" = "#21854e",
# 	    "Repressive" = "#bc3b2a",
# 	    "Bivalent" = "#e18727",
# 	    "Un" = "#1872b5"
#     )) +
#     labs(x = "> Number of cell types with repressive state", y = "Proportion of genes")

p = ggplot(d_active_count[hspc_cat %in% c("Active", "Repressive", "Un", "Bivalent") & active_n >= 1], 
    aes(x = hspc_cat, y = n, fill = hspc_cat)) +
	geom_bar(stat = "identity", position = "dodge") +
	theme_classic() +
	scale_fill_manual(values = c(
	    "Active" = "#21854e",
	    "Repressive" = "#bc3b2a",
	    "Bivalent" = "#e18727",
	    "Un" = "#1872b5"
	)) +
	labs(x = "> Number of cell types with repressive state", y = "Proportion of genes")

ggsave("./figures/227_figure_chromatin_status_x2repressive_proportion_active.pdf", p, width = 5, height = 3)


d = fread("./tmp/table_gene_chromatin_status_consensus_by_celltype_wide.tsv")
d_status_m = d[, c(-1, -5), with = T]
d_status_m_B = d_status_m[, c(4, 5, 6, 7, 11, 12), with = F]
d_status_m_M = d_status_m[, c(1:3), with = F]
d_status_m_E = d_status_m[, c(8:10), with = F]

d_status_m_l = list(
    B = d_status_m_B,
    M = d_status_m_M,
    E = d_status_m_E
)

lapply(d_status_m_l, function(x) {
	apply(x, 1, function(x) (sum(x == "Repressive", na.rm = T)))
}) -> repressive_count_list
names(repressive_count_list) = names(d_status_m_l)

lapply(1:length(repressive_count_list), function(i) {
    active_count = repressive_count_list[[i]]
    res = d[, .(hspc_cat = HSPC, active_n = active_count)] %>%
	.[, .(n = sum(active_n >= 1) / .N), by = hspc_cat]
    res$lineage = names(repressive_count_list)[i]
    res
}) %>% rbindlist() -> d_active_count


p = ggplot(d_active_count[hspc_cat %in% c("Active", "Un", "Bivalent")],
    aes(x = lineage, y = n, fill = hspc_cat)) +
    geom_bar(stat = "identity", position = "dodge") +
	theme_classic() +
	scale_fill_manual(values = c(
	    "Active" = "#21854e",
	    "Repressive" = "#bc3b2a",
	    "Bivalent" = "#e18727",
	    "Un" = "#1872b5"
	)) +
	labs(x = "> Number of cell types with repressive state", y = "Proportion of genes")

ggsave("./figures/227_figure_chromatin_status_x2repressive_proportion_active.pdf", p, width = 5, height = 3)


## }}}

############
#  Backup  #
############

## Trajectory {{{
# d_u2a = fread("./tmp/raw_table_B_cell_Avg_u2a_genes_missing_7.tsv")
# d_b2a = fread("./tmp/raw_table_B_cell_Avg_b2a_genes_missing_7.tsv")
#
# d_plot = rbind(
# 	d_u2a[, .(gene_score, pseudotime, group = "Un->Active")],
# 	d_b2a[, .(gene_score, pseudotime, group = "Bivalent->Active")]
# )
#
# p = ggplot(d_plot, aes(x = pseudotime, y = gene_score, color = group)) +
#     geom_point(shape = ".", alpha = 0.3) +
#     geom_smooth(method = "loess", se = F) +
#     theme_classic() +
#     scale_color_manual(values = c("Un->Active" = "#1872b5", "Bivalent->Active" = "#e18727")) +
#     labs(x = "Pseudotime", y = "Gene score")
# ggsave("./figures/227_figure_chromatin_status_x2active_trajectory.pdf",p, width = 5, height = 3)
## }}}

### {{{
# v_bi2active = readLines("../data/All_Biv_Active_genes.txt")
# v_un2active = readLines("../data/All_Un_Active_genes.txt")
#
# library(clusterProfiler)
# library(org.Hs.eg.db)
# eg = bitr(v_bi2active, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# ego_bi2active = enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = T)
# ego_bi2active_df = as.data.table(ego_bi2active)
# ego_bi2active_df[1:20, .(Description, GeneRatio, qvalue)]
#
# eg = bitr(v_un2active, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# ego_un2active = enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = T)
# ego_un2active_df = as.data.table(ego_un2active)
# ego_un2active_df[1:20, .(Description, GeneRatio, qvalue)]
#
# d_plot = rbind(
#     ego_bi2active_df[, .(Description, FoldEnrichment, qvalue, group = "Bivalent->Active")],
#     ego_un2active_df[, .(Description, FoldEnrichment = -FoldEnrichment, qvalue, group = "Un->Active")]
# )
#
# library(ggrepel)
#
# g = ggplot(d_plot, aes(x = FoldEnrichment, y = -log10(qvalue), color = group)) +
# 	geom_point() +
# 	theme_bw() +
# 	# lims(x = c(-9, 9)) +
# 	## Add label for top 10 terms in each group
# 	geom_label_repel(data = d_plot[qvalue < 0.05][order(qvalue)][, head(.SD, 20), by = group], aes(label = Description), size = 1) +
# 	theme_classic()
#
# ggsave("./figures/227_figure_chromatin_status_x2active_enrichment.pdf", g, width = 6, height = 6)
#
#
#
# ## KEGG enrichment
# eg = bitr(v_bi2active, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# ego_bi2active_kegg = enrichKEGG(eg$ENTREZID, organism = "hsa", pAdjustMethod = "BH", qvalueCutoff = 0.05)
# ego_bi2active_kegg_df = as.data.table(ego_bi2active_kegg)
# ego_bi2active_kegg_df[1:20, .(Description, GeneRatio, qvalue)]
#
# eg = bitr(v_un2active, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# ego_un2active_kegg = enrichKEGG(eg$ENTREZID, organism = "hsa", pAdjustMethod = "BH", qvalueCutoff = 0.05)
# ego_un2active_kegg_df = as.data.table(ego_un2active_kegg)
# ego_un2active_kegg_df[1:20, .(Description, GeneRatio, qvalue)]
#
# d_plot_kegg = rbind(
# 	ego_bi2active_kegg_df[, .(Description, FoldEnrichment, qvalue, group = "Bivalent->Active")],
# 	ego_un2active_kegg_df[, .(Description, FoldEnrichment = -FoldEnrichment, qvalue, group = "Un->Active")]
# )
#
# g_kegg = ggplot(d_plot_kegg, aes(x = FoldEnrichment, y = -log10(qvalue), color = group)) +
# 	geom_point() +
# 	theme_bw() +
# 	# lims(x = c(-9, 9)) +
# 	## Add label for top 10 terms in each group
# 	geom_label_repel(data = d_plot_kegg[qvalue < 0.05][order(qvalue)][, head(.SD, 20), by = group], aes(label = Description), size = 1) +
# 	theme_classic()
#
# ggsave("./figures/227_figure_chromatin_status_x2active_kegg_enrichment.pdf", g_kegg, width = 6, height = 6)
#
# ## GSEA C2 enrichment
# library(clusterProfiler)
# library(org.Hs.eg.db)
#
# term2gene = msigdbr(species = "Homo sapiens", category = "C2") %>%
#     dplyr::select(gs_name, ncbi_gene) %>%
#     dplyr::distinct()
#
# eg = bitr(v_bi2active, fromType = "SYMBOL", toType   = "ENTREZID", OrgDb    = org.Hs.eg.db)
# ego_bi2active_c2 = enricher( gene = eg$ENTREZID, TERM2GENE = term2gene, pAdjustMethod = "BH", qvalueCutoff  = 0.05)
# ego_bi2active_c2_df = as.data.table(ego_bi2active_c2)
# ego_bi2active_c2_df[1:20, .(Description, qvalue)]
#
# eg = bitr(v_un2active, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# ego_un2active_c2 = enricher( gene = eg$ENTREZID, TERM2GENE = term2gene, pAdjustMethod = "BH", qvalueCutoff  = 0.05)
# ego_un2active_c2_df = as.data.table(ego_un2active_c2)
# ego_un2active_c2_df[1:20, .(Description, qvalue)]
#
# d_plot_c2 = rbind(
#     ego_bi2active_c2_df[, .(Description, FoldEnrichment, qvalue, group = "Bivalent->Active")],
#     ego_un2active_c2_df[, .(Description, FoldEnrichment = -FoldEnrichment, qvalue, group = "Un->Active")]
# )
#
# g_c2 = ggplot(d_plot_c2, aes(x = FoldEnrichment, y = -log10(qvalue), color = group)) +
# 	geom_point() +
# 	theme_bw() +
# 	# lims(x = c(-9, 9)) +
# 	## Add label for top 10 terms in each group
# 	geom_label_repel(data = d_plot_c2[qvalue < 0.05][order(qvalue)][, head(.SD, 20), by = group], aes(label = Description), size = 1) +
# 	theme_classic()
#
# ggsave("./figures/227_figure_chromatin_status_x2active_c2_enrichment.pdf", g_c2, width = 6, height = 6)
#
# ## GSEA C8 enrichment
# library(clusterProfiler)
# library(org.Hs.eg.db)
#
# term2gene = msigdbr(species = "Homo sapiens", category = "C8") %>%
#     dplyr::select(gs_name, ncbi_gene) %>%
#     dplyr::distinct()
#
# eg = bitr(v_bi2active, fromType = "SYMBOL", toType   = "ENTREZID", OrgDb    = org.Hs.eg.db)
# ego_bi2active_c8 = enricher( gene = eg$ENTREZID, TERM2GENE = term2gene, pAdjustMethod = "BH", qvalueCutoff  = 0.05)
# ego_bi2active_c8_df = as.data.table(ego_bi2active_c8)
# ego_bi2active_c8_df[1:20, .(Description, qvalue)]
#
# eg = bitr(v_un2active, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# ego_un2active_c8 = enricher( gene = eg$ENTREZID, TERM2GENE = term2gene, pAdjustMethod = "BH", qvalueCutoff  = 0.05)
# ego_un2active_c8_df = as.data.table(ego_un2active_c8)
# ego_un2active_c8_df[1:20, .(Description, qvalue)]
#
# d_plot_c8 = rbind(
#     ego_bi2active_c8_df[, .(Description, FoldEnrichment, qvalue, group = "Bivalent->Active")],
#     ego_un2active_c8_df[, .(Description, FoldEnrichment = -FoldEnrichment, qvalue, group = "Un->Active")]
# )
#
# g_c8 = ggplot(d_plot_c8, aes(x = FoldEnrichment, y = -log10(qvalue), color = group)) +
# 	geom_point() +
# 	theme_bw() +
# 	# lims(x = c(-9, 9)) +
# 	## Add label for top 10 terms in each group
# 	geom_label_repel(data = d_plot_c8[qvalue < 0.05][order(qvalue)][, head(.SD, 20), by = group], aes(label = Description), size = 1) +
# 	theme_classic()
#
# ggsave("./figures/227_figure_chromatin_status_x2active_c8_enrichment.pdf", g_c8, width = 6, height = 6)
#
# ## CPG island enrichment
# library(ggsci)
# d_cpg_gene = readLines("./tmp/table_gene_cpg_island_promoter_hg19.txt")
#
# d_plot_cpg = data.table(
# 	gene = c(v_bi2active, v_un2active),
# 	group = c(rep("Bivalent->Active", length(v_bi2active)), rep("Un->Active", length(v_un2active)))
# )
#
# d_plot_cpg[, is_cpg := ifelse(gene %in% d_cpg_gene, "With CpG island", "Without CpG island")]
#
# g = ggplot(d_plot_cpg, aes(x = group, fill = is_cpg)) +
# 	geom_bar(position = "fill") +
# 	scale_fill_nejm(name = "CpG island in promoter") +
# 	theme_classic() +
# 	theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
#
# ggsave("./figures/227_figure_chromatin_status_x2active_cpg_enrichment.pdf", g, width = 4, height = 4)
#
#
### }}}
