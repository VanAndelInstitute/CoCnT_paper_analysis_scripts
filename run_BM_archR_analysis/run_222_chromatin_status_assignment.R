# =============================================================================
# Script: run_222_chromatin_status_assignment.R
#
# Input files:
#   - ./tmp/table_gene_marker_average_score_by_celltype.tsv.gz
#       (Table containing gene marker average scores by cell type.)
#   - ./tmp/table_gene_cpg_island_promoter_hg19.txt
#       (List of genes with CpG island promoter annotation.)
#
# Output file:
#   - ./tmp/table_gene_chromatin_status_by_celltype_me1.tsv
#   - ./tmp/table_gene_chromatin_status_by_celltype_me2.tsv
#   - ./tmp/table_gene_chromatin_status_by_celltype_me3.tsv
#   - ./tmp/table_gene_chromatin_status_wide_me1.tsv
#   - ./tmp/table_gene_chromatin_status_wide_me2.tsv
#   - ./tmp/table_gene_chromatin_status_wide_me3.tsv
#
# Description:
#   - This script reads in marker gene average scores and the list of genes with CpG islands.
#   - It annotates genes with CpG island promoter status.
#   - It assigns each gene to a chromatin status (Bivalent, Active, Repressive, Others)
#     for each cell type and writes the result to the output file for downstream analysis.
# =============================================================================


library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(ggsci)

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

cutoff = c(
    H3K27me3 = 0.21,
    H3K4me1 = 0.305,
    H3K4me2 = 0.139,
    H3K4me3 = 0.181,
    H3K4me1_cooc = 0.447,
    H3K4me2_cooc = 0.314,
    H3K4me3_cooc = 0.455
)

ct_v_B = c(
    "HSPC",
    "Pre-Pro-B cells",
    "Immature B cells",
    "Mature B cells1",
    "Mature B cells2",
    "Memory B cells",
    "Plasma cells"
)

ct_v_E = c(
    "HSPC",
    "Erythroid progenitors1",
    "Erythroid progenitors2",
    "Erythroid progenitors3"
)

ct_v_M = c(
    "HSPC",
    "Myelocyte/Classical Monocytes1",
    "Myelocyte/Classical Monocytes2",
    "Myelocyte/Classical Monocytes3"
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

## only include the choosen cell types
d_plot = d_plot[cell_type %in% ct_by_diff]

assign_chromatin_status <- function(d_plot, histone_marks, cooc_mark, cutoff, output_prefix) {
    d_sub <- d_plot[marker %in% c(histone_marks, "H3K27me3", cooc_mark), ]
    d_wide <- dcast(d_sub, gene + cell_type + is_cpg ~ marker, value.var = "gene_score")

    # Bivalent: co-occurrence AND both marks high
    d_wide[get(cooc_mark) > cutoff[cooc_mark] & get(histone_marks) > cutoff[histone_marks] & H3K27me3 > cutoff["H3K27me3"], cat := "Bivalent"]
    # Active: histone mark high, H3K27me3 low
    d_wide[get(histone_marks) > cutoff[histone_marks] & H3K27me3 < cutoff["H3K27me3"], cat := "Active"]
    # Repressive: histone mark low, H3K27me3 high
    d_wide[get(histone_marks) < cutoff[histone_marks] & H3K27me3 > cutoff["H3K27me3"], cat := "Repressive"]
    # Un: both low
    d_wide[get(histone_marks) < cutoff[histone_marks] & H3K27me3 < cutoff["H3K27me3"], cat := "Un"]
    # Other: any not classified above
    d_wide[is.na(cat), cat := "Other"]

    unique_other_genes <- unique(d_wide[cat == "Other", gene])

    write_tsv(d_wide, sprintf("./tmp/table_gene_chromatin_status_by_celltype_%s.tsv", output_prefix))

    d_wide <- fread(sprintf("./tmp/table_gene_chromatin_status_by_celltype_%s.tsv", output_prefix))

    # Each row: gene; columns: cell types
    d_wide_w <- dcast(d_wide, gene ~ cell_type, value.var = "cat")
    setcolorder(d_wide_w, c("gene", ct_by_diff))

    write_tsv(d_wide_w, sprintf("./tmp/table_gene_chromatin_status_wide_%s.tsv", output_prefix))
}

## Run for ME1, ME2, ME3
assign_chromatin_status(d_plot, "H3K4me1", "H3K4me1_cooc", cutoff, "me1")
assign_chromatin_status(d_plot, "H3K4me2", "H3K4me2_cooc", cutoff, "me2")
assign_chromatin_status(d_plot, "H3K4me3", "H3K4me3_cooc", cutoff, "me3")


### {{{
# plot_cat_gene_exp = function(d_wide) {
#     d_wide$cat = factor(d_wide$cat, levels = c("Active", "Repressive", "Bivalent", "Un", "Other"))

#     d_wide$cell_type = factor(d_wide$cell_type, levels = ct_by_diff)

#     d_wide = d_wide %<>% na.omit

#     # g1 = ggplot(d_wide, aes(x = cat, y = rna_expression + 0.0001, fill = cell_type)) +
# 	    # geom_boxplot(outliers = F) +
# 	    # scale_fill_nejm() +
# 	    # scale_y_log10() +
# 	    # theme_classic() 

#     d_wide_count = d_wide[, .N, by = .(cell_type, cat)]
#     g2 = ggplot(d_wide_count, aes(x = cell_type, y = N, fill = cat)) +
# 	    geom_bar(stat = "identity", position = "fill") +
# 	    scale_fill_manual(
# 		values = c(
# 		    "Active" = "#21854e",
# 		    "Repressive" = "#bc3b2a",
# 		    "Bivalent" = "#e18727",
# 		    "Un" = "#1872b5",
# 		    "Other" = "grey"
# 		)) +
# 	    theme_classic() +
# 	    theme(axis.text.x = element_text(angle = 45, hjust = 1))

#     print(g2)
# }
# plot_cat_gene_exp(d_wide[cell_type %in% ct_v_B & gene %in% hspc_gene_list[cat == "Bivalent"]$gene])
# plot_cat_gene_exp(d_wide[is_cpg == F])
#
# ## Dynamics of chromatin status along the Erythroid differentiation trajectory
# # NOTE: change here
# ct_set = ct_v_E
# d_sub = d_wide[cell_type %in% ct_set]
# d_sub$cell_type = factor(d_sub$cell_type, levels = ct_set)
#
# ## 1. sort by gene, then by cell_type
# setorder(d_sub, gene, cell_type)
#
# ## 2. for each gene, get the state in the *next* cell type
# d_sub[, `:=`(
#   next_cat  = shift(cat, type = "lead"),
#   next_cell = shift(as.character(cell_type), type = "lead")
# ), by = gene]
#
# ## keep only real transitions (within the trajectory)
# transitions <- d_sub[!is.na(next_cat)]
#
# write_tsv(transitions, "./tmp/table_chromatin_status_transitions_E_me1.tsv")
#
# ## 3a. pooled transition counts & probabilities (all steps together)
# trans_probs <- transitions[
#   , .N, by = .(from = cat, to = next_cat)
# ][
#   , total_from := sum(N), by = from
# ][
#   , prob := N / total_from
# ][order(from, to)]
#
# ## 3b. if you want probabilities *per step* (HSPC→Ery1, Ery1→Ery2, …):
# transitions[, step := paste0(cell_type, "->", next_cell)]
#
# trans_probs_by_step <- transitions[
#   , .N, by = .(step, from = cat, to = next_cat)
# ][
#   , total_from := sum(N), by = .(step, from)
# ][
#   , prob := N / total_from
# ][order(step, from, to)]
#
# trans_probs_by_step[prob > 0.10 & from != "Other" & to != "Other"]
# trans_probs_by_step[from == "Bivalent" & to == "Active"][order(to, step)]
# trans_probs_by_step[from == "Un" & to == "Active"][order(to, step)]
# trans_probs_by_step[from == "Un" & to == "Repressive"][order(to, step)]
#
# trans_probs_by_step$step %<>% factor(levels = ct_set[-length(ct_set)] %>% 
#   paste0("->", ct_set[-1]))
#
# ggplot(trans_probs_by_step, aes(x = step, y = prob, group = interaction(from, to), color = to)) +
#   geom_line() +
#   geom_point() +
#   scale_color_manual(
# 	values = c(
# 	  "Active" = "#21854e",
# 	  "Repressive" = "#bc3b2a",
# 	  "Bivalent" = "#e18727",
# 	  "Un" = "#1872b5"
# 	)
#   ) +
#   facet_wrap(~from) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   ggtitle("Chromatin status transition probabilities along B cell differentiation trajectory")
#
# ## NOTE: Change here
# write_tsv(trans_probs_by_step, "./tmp/table_chromatin_status_transition_probabilities_by_step_E_me3.tsv")
#
# ## Get gene list
# gene_lists_by_step <- transitions[
#   , .(genes = list(unique(gene))),       # collect genes as a list column
#   by = .(step, from = cat, to = next_cat)
# ][order(step, from, to)]
#
# gene_list1 = c("MECOM", "PRDM16", "HOXB5", "PAX5", "EBF1", "CEBPD", "CEBPA", "CEBPE", "FOSB", "JUNB", "ZFPM1", "GATA1")
#
# gene_lists_by_step
#
# ## NOTE: Change here
# write_rds(gene_lists_by_step, "./tmp/rds_gene_lists_by_step_M_me3.rds")
#
# out = gene_lists_by_step[, .(gene = unlist(genes)), by =.(step, from, to)]#[gene %in% gene_list1]
# ## NOTE: Change here
# write_tsv(out, "./tmp/table_selected_genes_chromatin_status_transitions_M_all_me3.tsv")
#
#
## }}}
