library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)


gene_category = fread("./tmp/table_gene_chromatin_status_trajectory_consensus_remove_dup_long.tsv")
d_cg_gene = readLines("./tmp/table_gene_cpg_island_promoter_hg19.txt")

## get the cpg island promoter genes percentage for each categories

gene_category[, is_cpg_promoter := ifelse(gene %in% d_cg_gene, "cpg_promoter", "non_cpg_promoter")]
gene_category[, lineage := substr(group, 1, 1)]


pdf("./figures/238_active2repressive_cpg_percentage.pdf", width = 8, height = 6)

transitions <- c("Active->Un->Repressive", "Active->Bivalent->Repressive", "Active->Repressive")
# marks <- c("me2", "me1", "me3")

# for (mark in marks) {
  cpg_percent <- gene_category[,
    # grepl("\\+", marker),
    .(cpg_promoter_percentage = sum(is_cpg_promoter == "cpg_promoter") / .N * 100),
    by = .(lineage, trans)
  ][trans %in% transitions] # & lineage != "M"]

  print(
    ggplot(cpg_percent, aes(x = lineage, y = cpg_promoter_percentage, fill = trans)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Transition", y = "Percentage of CpG Island Promoter Genes") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
  )
# }

dev.off()
