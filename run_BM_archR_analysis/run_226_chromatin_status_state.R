library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

f_v = dir("./tmp/", "table_gene_chromatin_status_by_celltype") %>% grep("123", ., invert = TRUE, value = T) 
d = lapply(f_v, function(f) {
    d = fread(paste0("./tmp/", f))
    names(d) = c("gene", "celltype", "is_cpg", "K27", "K4", "cooc", "cat")
    d[, mark := str_extract(f, "me[123]")]
    return(d)
}) %>% rbindlist()

# plot(1)
# ggplot(d) + aes(x = K27, y = K4, color = cat) + geom_point() + facet_wrap(~mark) +
    # scale_y_log10() + scale_x_log10() + theme_bw() + theme(panel.grid = element_blank()) 
# d[K27 > 10 & K4 > 10, gene] %>% unique

d_consensus = d[, .(consensus_count = .N, consensus = paste(sort(mark), collapse = "+")), by = .(gene, celltype, cat)]
d_consensus

d_consensus$cell_type = d_consensus$celltype

write_tsv(d_consensus, "./tmp/table_gene_chromatin_status_consensus_by_celltype.tsv")

d_consensus_w = dcast(d_consensus[consensus_count > 1], gene ~ celltype, value.var = "cat")

write_tsv(d_consensus_w, "./tmp/table_gene_chromatin_status_consensus_by_celltype_wide.tsv")


