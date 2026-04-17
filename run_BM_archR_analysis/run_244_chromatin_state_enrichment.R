library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

d_state = fread("./tmp/table_gene_chromatin_status_consensus_by_celltype_wide.tsv") 

d = d_state[, .(gene, HSPC)]

## do the Gene ontology enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)

## convert gene symbol to entrez id
gene_symbol = d$gene
gene_entrez = bitr(gene_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_entrez = gene_entrez[!duplicated(gene_entrez$ENTREZID), ]

d = merge(d, gene_entrez, by.x = "gene", by.y = "SYMBOL")

## run the enrichment analysis

## Bivalent
x = d[HSPC == "Bivalent", ]$ENTREZID
ego = enrichGO(gene = x, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 1,qvalueCutoff = 1, readable = TRUE, maxGSSize = 1000, minGSSize = 50)
ego_result = as.data.frame(ego)
ego_result_b = ego_result[order(ego_result$p.adjust), ]

## Active
x = d[HSPC == "Active", ]$ENTREZID
ego = enrichGO(gene = x, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE, maxGSSize = 1000, minGSSize = 50)
ego_result = as.data.frame(ego)
ego_result_a = ego_result[order(ego_result$p.adjust), ]

## Repressive
x = d[HSPC == "Repressive", ]$ENTREZID
ego = enrichGO(gene = x, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE, maxGSSize = 1000, minGSSize = 50)
ego_result = as.data.frame(ego)
ego_result_r = ego_result[order(ego_result$p.adjust), ]

## Unmarked
x = d[HSPC == "Un", ]$ENTREZID
ego = enrichGO(gene = x, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE, maxGSSize = 1000, minGSSize = 50)
ego_result = as.data.frame(ego)
ego_result_u = ego_result[order(ego_result$p.adjust), ]


## Put together
ego_result_b$HSPC_state <- "Bivalent"
ego_result_a$HSPC_state <- "Active"
ego_result_r$HSPC_state <- "Repressive"
ego_result_u$HSPC_state <- "Unmarked"

ego_all <- rbind(
    ego_result_b,
    ego_result_a,
    ego_result_r,
    ego_result_u
) %>% data.table

## choose all the terms under the development process
library(GO.db)

# dev_go <- c("GO:0032502")
dev_go = "GO:0030097"

# get offspring terms
dev_children <- unique(c(
    dev_go,
    unlist(as.list(GOBPOFFSPRING[dev_go]))
))

ego_all_dev <- ego_all[ID %in% dev_children, ]
ego_all_dev$HSPC_state %>% table

library(dplyr)
library(ggplot2)

plot_df <- as.data.table(ego_all_dev) %>%
    group_by(HSPC_state) %>%
    slice_min(order_by = pvalue, n = 45, with_ties = FALSE) %>%
    ungroup()

# httpgd::hgd(port = 4322)
colors_cat <- c(
    "Active" = "#21854e",
    "Repressive" = "#bc3b2a",
    "Bivalent" = "#e18727",
    "Un" = "#1872b5"
)

write_tsv(plot_df, "./tmp/table_HSPC_state_GO_enrichment.tsv")

plot_df = fread("./tmp/table_HSPC_state_GO_enrichment.tsv")

httpgd::hgd(port = 4322)

pdf("./figures/244_figure_chromatin_state_enrichment.pdf", width = 10, height = 7)

ggplot(plot_df, aes(x = FoldEnrichment, y = -log10(pvalue))) + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
    facet_wrap(~HSPC_state, ncol = 1) +
    theme_classic() +
    geom_point()


ggplot(plot_df[HSPC_state %in% c("Bivalent", "Unmarked")], aes(x = FoldEnrichment, y = -log10(pvalue))) + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
    facet_wrap(~HSPC_state, ncol = 1) +
    theme_classic() +
    geom_point()

term_order = plot_df[HSPC_state == "Bivalent"][order(FoldEnrichment), Description]

plot_df$Description = factor(plot_df$Description, levels = term_order)

ggplot(plot_df, aes(x = FoldEnrichment, y = Description)) +
    geom_point(aes(size = -log10(pvalue), color = HSPC_state)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey") +
    facet_wrap(~HSPC_state, nrow = 1) +
    scale_color_manual(values = colors_cat) +
    theme_classic() +
    ylab("") +
    xlab("Fold Enrichment") 
dev.off()
