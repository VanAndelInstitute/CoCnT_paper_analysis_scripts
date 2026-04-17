################################################################################
# Script: run_n_rna_pseudo_bulk_match_cocnt_cluster.R
# Purpose: Map Triana 2021 scRNA-seq cell-type annotations to coCnT-style merged
#          cell classes and generate RNA pseudo-bulk expression matrices for
#          both merged and original annotations
#
# INPUT:
#   - Seurat object:
#     * ../data/triana_2021/WTA_projected.bk.rds
#
# OUTPUT:
#   - RNA pseudo-bulk tables:
#     * ./tmp/table_triana_2021_rna_pseudo_bulk_by_ct2_merge.tsv
#     * ./tmp/table_triana_2021_rna_pseudo_bulk_by_ct.tsv
################################################################################

library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

library(Seurat)

so = read_rds("../../data/triana_2021/WTA_projected.bk.rds")

Idents(so) = so@meta.data$ct

x = so@meta.data$ct %>% table
selected = names(x[x > 10])
y = x[x > 10]
y / sum(y)
so_sub = subset(so, idents = selected)

d_ann = so_sub@meta.data
d_ann$ct %>% table

map_rules <- data.table(
    pattern = c(
	# B lineage (stage-specific)
	"(?i)plasma",                                  # Plasma cells
	"(?i)immature\\s*b",                           # Immature B
	"(?i)class\\s*switched|non\\s*switched|memory\\s*b",  # Memory B
	"(?i)mature\\s*naive\\s*b|mature\\s*b",        # Mature/naive B
	"(?i)pre[- ]?pro[- ]?b|pre[- ]?b|pro[- ]?b|small\\s*pre\\s*b", # Pre/Pro-B

	# T / NK
	"(?i)gamma\\s*delta|cd4\\+|cd8\\+|t\\s*cells?|t\\s*cell", # T cell family
	"(?i)^nk\\b|nk\\s*cells?|nk\\s*t|cd56",         # NK & NKT (group to NK Cells)
	"(?i)nk\\s*progenitor",                         # NK progenitors -> NK Cells

	# DCs
	"(?i)Plasmacytoid\\s*dendritic",                # pDC (and pDC progenitors)

	# Monocytes / Granulocytes
	"(?i)non[- ]?classical\\s*monocyte",            # Non-classical Monocytes
	"(?i)classical\\s*monocyte",                    # Classical Monocytes
	"(?i)myelocyte|promyelocyte",                   # Myelocyte/Promyelocyte
	"(?i)myelocyte/classical\\s*monocytes?",        # mixed label in coCnT

	# Erythroid / Meg / Eo-Baso-Mast
	"(?i)\\s*erythroid",                       	    # Early erythroid prog
	"(?i)megakaryocyte\\s*progenitor",              # Meg progenitors
	"(?i)eosinophil|basophil|eobaso",          	    # Eo/Baso/Mast progenitors

	# Stem/progenitor umbrella
	"(?i)hsc|mpp|hspc",                             # HSPC/MPP

	# catch-all NA/empty
	"^\\s*$"
	),
    label = c(
	"Plasma cells",
	"Immature B cells",
	"Mature B cells",
	"Mature B cells",
	"Pre/Pro-B cells",

	"T cells",
	"NK Cells",
	"NK Cells",

	"pDC",

	"Non-classical Monocytes",
	"Myelocyte/Classical Monocytes",
	"Myelocyte/Classical Monocytes",
	"Myelocyte/Classical Monocytes",

	"Erythroid progenitors",
	"HSPC",
	"EoBaso progenitors",

	"HSPC",

	"Other"
    )
)

## helper to apply mapping
map_ann <- function(x, rules = map_rules) {
    y <- rep(NA_character_, length(x))
    x0 <- ifelse(is.na(x), "", x)
    for (i in seq_len(nrow(rules))) {
	hit <- str_detect(x0, rules$pattern[i])
	y[is.na(y) & hit] <- rules$label[i]
    }
    # Anything still NA after rules → "Other"
    y[is.na(y)] <- "Other"
    y
}

## Apply mapping
d_ann$ct2_merge = map_ann(as.character(d_ann$ct))
unique(d_ann$ct)
unique(d_ann$ct2_merge)

so_sub@meta.data = d_ann

# Aggregate to pseudo-bulk by coCut&Tag cluster
so_pb = AverageExpression(so_sub, group.by = "ct2_merge", assays = "RNA", slot = "data")
d_rna_bulk = data.table(as.data.frame(so_pb$RNA), keep.rownames = "gene")
write_tsv(d_rna_bulk, "./tmp/table_triana_2021_rna_pseudo_bulk_by_ct2_merge.tsv")

# Aggregate to pseudo-bulk by original cell type
so_pb2 = AverageExpression(so_sub, group.by = "ct", assays = "RNA", slot = "data")
d_rna_bulk2 = data.table(as.data.frame(so_pb2$RNA), keep.rownames = "gene")
write_tsv(d_rna_bulk2, "./tmp/table_triana_2021_rna_pseudo_bulk_by_ct.tsv")
