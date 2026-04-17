library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(GenomicRanges)
httpgd::hgd(port = 4322)


## gene body annotation
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
gr_gene <- genes(txdb)   # GRanges: seqnames, start, end, strand, gene_id
gene_body_ann <- as.data.frame(gr_gene)[, c("seqnames","start","end","strand","gene_id")] %>% data.table
gene_body_ann$gene_name <- mapIds(org.Hs.eg.db, keys = gene_body$gene_id, column = "SYMBOL", keytype = "ENTREZID")
gene_body_ann = gene_body_ann[seqnames %in% paste0("chr", c(1:22, "X", "Y"))]

## celltype order {{{
cell_type_order = c(
    "HSPC",
    "Pre-Pro-B cells",
    "Immature B cells",
    "Mature B cells1",
    "Mature B cells2",
    "Memory B cells",
    "Plasma cells"
    # "Myelocyte/Classical Monocytes1",
    # "Myelocyte/Classical Monocytes2",
    # "Myelocyte/Classical Monocytes3",
    # "Erythroid progenitors1",
    # "Erythroid progenitors2",
    # "Erythroid progenitors3"
)


## }}}

## color {{{
annot_ct2 = c(
    "C1" = "NK cells1",
    "C2" = "T cells",
    "C3" = "Mature B cells1",
    "C4" = "Plasma cells",
    "C5" = "Mature B cells2",
    "C6" = "Memory B cells",
    "C7" = "Erythroid progenitors3",
    "C8" = "Erythroid progenitors2",
    "C9" = "NK cells2",
    "C10" = "HSPC",
    "C11" = "Plasmacytoid dendritic cells",
    "C12" = "Myelocyte/Classical Monocytes1",
    "C13" = "Myelocyte/Classical Monocytes2",
    "C14" = "Myelocyte/Classical Monocytes3",
    "C15" = "Non-classical Monocytes",
    "C16" = "Pre-Pro-B cells",
    "C17" = "Immature B cells",
    "C18" = "Erythroid progenitors1",
    "C19" = "EoBaso progenitors"
)

color <- c(
    C1  = "#8C8178",
    C2  = "#A3623A",
    C3  = "#B4009E",
    C4  = "#7D0033",
    C5  = "#89288F",
    C6  = "#B6B0FF",
    C7  = "#499195",
    C8  = "#208A42",
    C9  = "#4A5B5D",
    C10 = "#D51F26",
    C11 = "#FFBBED",
    C12 = "#8B3E1D",
    C13 = "#F47D2B",
    C14 = "#E09E5A",
    C15 = "#BDA687",
    C16 = "#6C2B85",
    C17 = "#8250C4",
    C18 = "#14532D",
    C19 = "#8AD4EB"
)

color_ct2 = color
names(color_ct2) = annot_ct2

## }}}

d_meta <- fread("./tmp/table_metadata_h3k27me3_final3_with_trajectories.tsv")

d_meta$ct2 %<>% factor(levels = cell_type_order)

# ggplot(d_meta) +
#     geom_point(aes(x = Bcell_Trajectory, y = PromoterRatio)) 
# ggsave("figures/test_promoter_ratio_trajectory.pdf", width = 4, height = 3)

sample_sheet = rbind(
    fread("../../data/sample_sheet_PR001798.tsv"),
    fread("../../data/sample_sheet_PR001799.tsv"),
    fread("../../data/sample_sheet_PR001855.tsv"),
    fread("../../data/sample_sheet_PR001856.tsv")
    )

d_bed = fread("../run_BM_bed_profile/tmp/table_good_cells_fragments.tsv")

gene_category = fread("./tmp/table_shared_genes_between_markers_in_trajectory.tsv")
gene_category_bi = gene_category[grepl("\\+", marker)]

# gene_ann = fread("../../data/hg38_TSS_list.bed")
# names(gene_ann) = c("seq", "start", "end", "gene_name", "strand")

d_cell_id_sample = d_bed[, .(cell_id, sample_name)] %>% unique

## B lineage, H3K27me3, CoCnT {{{

cell_type_order = c(
    "HSPC",
    "Pre-Pro-B cells",
    "Immature B cells",
    "Mature B cells1",
    "Mature B cells2",
    "Memory B cells",
    "Plasma cells"
    # "Myelocyte/Classical Monocytes1",
    # "Myelocyte/Classical Monocytes2",
    # "Myelocyte/Classical Monocytes3",
    # "Erythroid progenitors1",
    # "Erythroid progenitors2",
    # "Erythroid progenitors3"
)

antibody_target_x = "H3K27me3"
sample_type_x = "CoCnT"

sample_name_v = sample_sheet[antibody_target %in% antibody_target_x & library_type %in% sample_type_x, sample_name]

d_bed_sub = d_bed[sample_name %in% sample_name_v]

transition_v = c("whole_genome", "Active->Bivalent->Repressive", "Active->Repressive", "Active->Un->Repressive")

dm = expand.grid(
    cell_type_x = cell_type_order,
    transition_x = transition_v,
    lienage_x = c("B")
    )

# dm = dm[-which(dm$lienage == "E" & !(grepl("Ery", dm$cell_type_x))), ]
# dm = dm[-which(dm$lienage == "B" & (grepl("Ery", dm$cell_type_x))), ]


lapply(1:nrow(dm), function(i) {
    print(i)
    cell_type_x = dm$cell_type_x[i]
    transition_x = dm$transition_x[i]
    lineage_x = dm$lienage_x[i]

    # marker_x = "me2"

    if (transition_x == "whole_genome") {
	gene_v = gene_body_ann$gene_name
	gene_ann_sub = gene_body_ann[gene_name %in% gene_v]
    } else {
	gene_v = gene_category_bi[trans %in% transition_x & lineage == lineage_x, gene]
	gene_ann_sub = gene_body_ann[gene_name %in% gene_v]
    }

    cell_id_v = d_meta[ct2 %in% cell_type_x, cell_id]
    d_bed_sub2 = d_bed_sub[cell_id %in% cell_id_v]

    gr_promoter = GRanges(
	seqnames = gene_ann_sub$seq,
	ranges = IRanges(
	    start = ifelse(gene_ann_sub$strand == "+", gene_ann_sub$start - 1000, gene_ann_sub$end - 1000),
	    end = ifelse(gene_ann_sub$strand == "+", gene_ann_sub$start + 1000, gene_ann_sub$end + 1000)
	    ),
	strand = gene_ann_sub$strand,
	gene_name = gene_ann_sub$gene_name
    )

    gr_body = GRanges(
	seqnames = gene_ann_sub$seq,
	ranges = IRanges(
	    start = gene_ann_sub$start,
	    end = gene_ann_sub$end
	    ),
	strand = gene_ann_sub$strand,
	gene_name = gene_ann_sub$gene_name
    )

    gr_fragments = GRanges(
	seqnames = d_bed_sub2$chr,
	ranges = IRanges(
	    start = d_bed_sub2$start,
	    end = d_bed_sub2$end
	    ),
	cell_id = d_bed_sub2$cell_id
    )

    hits_prom <- countOverlaps(gr_promoter, gr_fragments, type = "any")
    hits_body <- countOverlaps(gr_body, gr_fragments, type = "any")

    data.table(
	cell_type = cell_type_x,
	transition = transition_x,
	hits_promoter = hits_prom,
	hits_body = hits_body,
	percentage_promoter = hits_prom / (hits_prom + hits_body + 1) * 100,
	gene = gene_ann_sub$gene_name
	)
    }) %>% rbindlist -> d_result

write_tsv(d_result, "./tmp/table_promoter_body_hits_percentage_h3k27me3_B.tsv")

d_result2 = d_result

d_cell_count = d_meta[, .(cell_n = .N), by = ct2]
d_result2 = merge(d_result2, d_cell_count, by.x = "cell_type", by.y = "ct2")
d_result2 = d_result2[(hits_promoter + hits_body) > (cell_n / 20)]
d_result2[, .N, by = .(cell_type, transition)] %>% print
d_result2$transition = factor(d_result2$transition, levels = transition_v)

p = ggplot(d_result2, aes(x = cell_type, y = percentage_promoter, fill = transition)) +
    geom_boxplot(outliers = F) +
    labs(x = "Cell Type", y = "Mean Percentage of Promoter Reads") +
    theme_classic() + ggtitle(paste0("Mean Percentage of Promoter Reads in Different Transitions (", antibody_target_x, ", ", sample_type_x, ")")) +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1))
p
ggsave(paste0("figures/239_promoter_reads_percent_boxplot_", antibody_target_x, "_", sample_type_x, ".pdf"), p, width = 8, height = 6)

d_plot = d_result2[, .(median = median(percentage_promoter), upper =  quantile(percentage_promoter, 0.75), 
    lower = quantile(percentage_promoter, 0.25)), by = .(cell_type, transition)]

ggplot(d_plot, aes(x = cell_type, y = median, color = transition, group = transition)) +
	geom_point(position = position_dodge(0.9)) +
	geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.9), width = 0.25) +
	geom_line(position = position_dodge(0.9)) +
	facet_wrap(~transition, ncol = 1) +
	labs(x = "Cell Type", y = "Mean Percentage of Promoter Reads") +
	theme_classic() + ggtitle(paste0("Mean Percentage of Promoter Reads in Different Transitions (", antibody_target_x, ", ", sample_type_x, ")")) +
	theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(paste0("figures/239_promoter_reads_percent_lineplot_", antibody_target_x, "_", sample_type_x, ".pdf"), width = 6, height = 6)

## }}}

## E lineage, H3K27me3, CoCnT {{{

cell_type_order = c(
    "HSPC",
    # "Pre-Pro-B cells",
    # "Immature B cells",
    # "Mature B cells1",
    # "Mature B cells2",
    # "Memory B cells",
    # "Plasma cells"
    # "Myelocyte/Classical Monocytes1",
    # "Myelocyte/Classical Monocytes2",
    # "Myelocyte/Classical Monocytes3",
    "Erythroid progenitors1",
    "Erythroid progenitors2",
    "Erythroid progenitors3"
)

antibody_target_x = "H3K27me3"
sample_type_x = "CoCnT"

sample_name_v = sample_sheet[antibody_target %in% antibody_target_x & library_type %in% sample_type_x, sample_name]

d_bed_sub = d_bed[sample_name %in% sample_name_v]

transition_v = c("whole_genome", "Active->Bivalent->Repressive", "Active->Repressive", "Active->Un->Repressive")

dm = expand.grid(
    cell_type_x = cell_type_order,
    transition_x = transition_v,
    lienage_x = c("E")
    )


lapply(1:nrow(dm), function(i) {
    print(i)
    cell_type_x = dm$cell_type_x[i]
    transition_x = dm$transition_x[i]
    lineage_x = dm$lienage_x[i]

    if (transition_x == "whole_genome") {
	gene_v = gene_body_ann$gene_name
	gene_ann_sub = gene_body_ann[gene_name %in% gene_v]
    } else {
	gene_v = gene_category_bi[trans %in% transition_x & lineage == lineage_x, gene]
	gene_ann_sub = gene_body_ann[gene_name %in% gene_v]
    }

    cell_id_v = d_meta[ct2 %in% cell_type_x, cell_id]
    d_bed_sub2 = d_bed_sub[cell_id %in% cell_id_v]

    gr_promoter = GRanges(
	seqnames = gene_ann_sub$seq,
	ranges = IRanges(
	    start = ifelse(gene_ann_sub$strand == "+", gene_ann_sub$start - 1000, gene_ann_sub$end - 1000),
	    end = ifelse(gene_ann_sub$strand == "+", gene_ann_sub$start + 1000, gene_ann_sub$end + 1000)
	    ),
	strand = gene_ann_sub$strand,
	gene_name = gene_ann_sub$gene_name
    )

    gr_body = GRanges(
	seqnames = gene_ann_sub$seq,
	ranges = IRanges(
	    start = gene_ann_sub$start,
	    end = gene_ann_sub$end
	    ),
	strand = gene_ann_sub$strand,
	gene_name = gene_ann_sub$gene_name
    )

    gr_fragments = GRanges(
	seqnames = d_bed_sub2$chr,
	ranges = IRanges(
	    start = d_bed_sub2$start,
	    end = d_bed_sub2$end
	    ),
	cell_id = d_bed_sub2$cell_id
    )

    hits_prom <- countOverlaps(gr_promoter, gr_fragments, type = "any")
    hits_body <- countOverlaps(gr_body, gr_fragments, type = "any")

    data.table(
	cell_type = cell_type_x,
	transition = transition_x,
	hits_promoter = hits_prom,
	hits_body = hits_body,
	percentage_promoter = hits_prom / (hits_prom + hits_body + 1) * 100,
	gene = gene_ann_sub$gene_name
	)
    }) %>% rbindlist -> d_result

write_tsv(d_result, "./tmp/table_promoter_body_hits_percentage_h3k27me3_E.tsv")

d_result2 = d_result

d_cell_count = d_meta[, .(cell_n = .N), by = ct2]
d_result2 = merge(d_result2, d_cell_count, by.x = "cell_type", by.y = "ct2")
d_result2 = d_result2[(hits_promoter + hits_body) > (cell_n / 20)]
d_result2[, .N, by = .(cell_type, transition)] %>% print
d_result2$transition = factor(d_result2$transition, levels = transition_v)

p = ggplot(d_result2, aes(x = cell_type, y = percentage_promoter, fill = transition)) +
    geom_boxplot(outliers = F) +
    labs(x = "Cell Type", y = "Mean Percentage of Promoter Reads") +
    theme_classic() + ggtitle(paste0("Mean Percentage of Promoter Reads in Different Transitions (", antibody_target_x, ", ", sample_type_x, ")")) +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1))
p
ggsave(paste0("figures/239_promoter_reads_percent_boxplot_", antibody_target_x, "_", sample_type_x, "_E.pdf"), p, width = 8, height = 6)

d_plot = d_result2[, .(median = median(percentage_promoter), upper =  quantile(percentage_promoter, 0.75), 
    lower = quantile(percentage_promoter, 0.25)), by = .(cell_type, transition)]

ggplot(d_plot, aes(x = cell_type, y = median, color = transition, group = transition)) +
	geom_point(position = position_dodge(0.9)) +
	geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.9), width = 0.25) +
	geom_line(position = position_dodge(0.9)) +
	facet_wrap(~transition, ncol = 1) +
	labs(x = "Cell Type", y = "Mean Percentage of Promoter Reads") +
	theme_classic() + ggtitle(paste0("Mean Percentage of Promoter Reads in Different Transitions (", antibody_target_x, ", ", sample_type_x, ")")) +
	theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(paste0("figures/239_promoter_reads_percent_lineplot_", antibody_target_x, "_", sample_type_x, "_E.pdf"), width = 6d, height = 6)

## }}}

## B lineage, H3K4me2, CoCnT {{{
antibody_target_x = "H3K4me2"
sample_type_x = "CoCnT"

sample_name_v = sample_sheet[antibody_target %in% antibody_target_x & library_type %in% sample_type_x, sample_name]

d_bed_sub = d_bed[sample_name %in% sample_name_v]

transition_v = c("whole_genome", "Active->Bivalent->Repressive", "Active->Repressive", "Active->Un->Repressive")

dm = expand.grid(
    cell_type_x = cell_type_order,
    transition_x = transition_v,
    lienage_x = c("B")
    )

# dm = dm[-which(dm$lienage == "E" & !(grepl("Ery", dm$cell_type_x))), ]
# dm = dm[-which(dm$lienage == "B" & (grepl("Ery", dm$cell_type_x))), ]


lapply(1:nrow(dm), function(i) {
    print(i)
    cell_type_x = dm$cell_type_x[i]
    transition_x = dm$transition_x[i]
    lineage_x = dm$lienage_x[i]

    marker_x = "me2"

    if (transition_x == "whole_genome") {
	gene_v = gene_body_ann$gene_name
	gene_ann_sub = gene_body_ann[gene_name %in% gene_v]
    } else {
	gene_v = gene_category[trans %in% transition_x & grepl(marker_x, marker) & lineage == lineage_x, gene]
	gene_ann_sub = gene_body_ann[gene_name %in% gene_v]
    }

    cell_id_v = d_meta[ct2 %in% cell_type_x, cell_id]
    d_bed_sub2 = d_bed_sub[cell_id %in% cell_id_v]

    gr_promoter = GRanges(
	seqnames = gene_ann_sub$seq,
	ranges = IRanges(
	    start = ifelse(gene_ann_sub$strand == "+", gene_ann_sub$start - 1000, gene_ann_sub$end - 1000),
	    end = ifelse(gene_ann_sub$strand == "+", gene_ann_sub$start + 1000, gene_ann_sub$end + 1000)
	    ),
	strand = gene_ann_sub$strand,
	gene_name = gene_ann_sub$gene_name
    )

    gr_body = GRanges(
	seqnames = gene_ann_sub$seq,
	ranges = IRanges(
	    start = gene_ann_sub$start,
	    end = gene_ann_sub$end
	    ),
	strand = gene_ann_sub$strand,
	gene_name = gene_ann_sub$gene_name
    )

    gr_fragments = GRanges(
	seqnames = d_bed_sub2$chr,
	ranges = IRanges(
	    start = d_bed_sub2$start,
	    end = d_bed_sub2$end
	    ),
	cell_id = d_bed_sub2$cell_id
    )

    hits_prom <- countOverlaps(gr_promoter, gr_fragments, type = "any")
    hits_body <- countOverlaps(gr_body, gr_fragments, type = "any")

    data.table(
	cell_type = cell_type_x,
	transition = transition_x,
	hits_promoter = hits_prom,
	hits_body = hits_body,
	percentage_promoter = hits_prom / (hits_prom + hits_body + 1) * 100,
	gene = gene_ann_sub$gene_name
	)
    }) %>% rbindlist -> d_result

write_tsv(d_result, "./tmp/table_promoter_body_hits_percentage_h3k4me2_B.tsv")

d_result2 = d_result

d_cell_count = d_meta[, .(cell_n = .N), by = ct2]
d_result2 = merge(d_result2, d_cell_count, by.x = "cell_type", by.y = "ct2")
d_result2 = d_result2[(hits_promoter + hits_body) > (cell_n / 20)]
d_result2[, .N, by = .(cell_type, transition)] %>% print
d_result2$transition = factor(d_result2$transition, levels = transition_v)

p = ggplot(d_result2, aes(x = cell_type, y = percentage_promoter, fill = transition)) +
    geom_boxplot(outliers = F) +
    labs(x = "Cell Type", y = "Mean Percentage of Promoter Reads") +
    theme_classic() + ggtitle(paste0("Mean Percentage of Promoter Reads in Different Transitions (", antibody_target_x, ", ", sample_type_x, ")")) +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1))
p
# ggsave(paste0("figures/239_promoter_reads_percent_boxplot_", antibody_target_x, "_", sample_type_x, ".pdf"), p, width = 8, height = 6)

d_plot = d_result2[, .(median = median(percentage_promoter), upper =  quantile(percentage_promoter, 0.75), 
    lower = quantile(percentage_promoter, 0.25)), by = .(cell_type, transition)]

ggplot(d_plot, aes(x = cell_type, y = median, color = transition, group = transition)) +
	geom_point(position = position_dodge(0.9)) +
	geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.9), width = 0.25) +
	geom_line(position = position_dodge(0.9)) +
	facet_wrap(~transition, ncol = 1) +
	labs(x = "Cell Type", y = "Mean Percentage of Promoter Reads") +
	theme_classic() + ggtitle(paste0("Mean Percentage of Promoter Reads in Different Transitions (", antibody_target_x, ", ", sample_type_x, ")")) +
	theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1))

# ggsave(paste0("figures/239_promoter_reads_percent_lineplot_", antibody_target_x, "_", sample_type_x, ".pdf"), width = 8, height = 6)

## }}}

## divide each gene into 100 bins and calculate the reads per bin

## B lineage, H3K27me3, CoCnT {{{

cell_type_order = c(
    "HSPC",
    "Pre-Pro-B cells",
    "Immature B cells",
    "Mature B cells1",
    "Mature B cells2",
    "Memory B cells",
    "Plasma cells"
    # "Myelocyte/Classical Monocytes1",
    # "Myelocyte/Classical Monocytes2",
    # "Myelocyte/Classical Monocytes3",
    # "Erythroid progenitors1",
    # "Erythroid progenitors2",
    # "Erythroid progenitors3"
)



antibody_target_x = "H3K27me3"
sample_type_x = "CoCnT"

sample_name_v = sample_sheet[antibody_target %in% antibody_target_x & library_type %in% sample_type_x, sample_name]

d_bed_sub = d_bed[sample_name %in% sample_name_v]

transition_v = c("whole_genome", "Active->Bivalent->Repressive", "Active->Repressive", "Active->Un->Repressive")

dm = expand.grid(
    cell_type_x = cell_type_order,
    transition_x = transition_v,
    lienage_x = c("B")
    )

# dm = dm[-which(dm$lienage == "E" & !(grepl("Ery", dm$cell_type_x))), ]
# dm = dm[-which(dm$lienage == "B" & (grepl("Ery", dm$cell_type_x))), ]

library(parallel)

lapply(1:nrow(dm), function(i) {
    print(i)
    cell_type_x = dm$cell_type_x[i]
    transition_x = dm$transition_x[i]
    lineage_x = dm$lienage_x[i]

    marker_x = "me2"

    if (transition_x == "whole_genome") {
	gene_v = gene_body_ann$gene_name
	gene_ann_sub = gene_body_ann[gene_name %in% gene_v]
    } else {
	gene_v = gene_category[trans %in% transition_x & grepl(marker_x, marker) & lineage == lineage_x, gene]
	gene_ann_sub = gene_body_ann[gene_name %in% gene_v]
    }

    cell_id_v = d_meta[ct2 %in% cell_type_x, cell_id]
    d_bed_sub2 = d_bed_sub[cell_id %in% cell_id_v]

    gr_fragments = GRanges(
	seqnames = d_bed_sub2$chr,
	ranges = IRanges(
	    start = d_bed_sub2$start,
	    end = d_bed_sub2$end
	    ),
	cell_id = d_bed_sub2$cell_id
    )

    mclapply(1:20, function(j) {
	gr_bin = GRanges(
	    seqnames = gene_ann_sub$seq,
	    ranges = IRanges(
		start = ifelse(gene_ann_sub$strand == "+", gene_ann_sub$start + (j-1) * (gene_ann_sub$end - gene_ann_sub$start) / 20, gene_ann_sub$end - j * (gene_ann_sub$end - gene_ann_sub$start) / 20),
		end = ifelse(gene_ann_sub$strand == "+", gene_ann_sub$start + j * (gene_ann_sub$end - gene_ann_sub$start) / 20, gene_ann_sub$end - (j-1) * (gene_ann_sub$end - gene_ann_sub$start) / 20)
		),
	    strand = gene_ann_sub$strand,
	    gene_name = gene_ann_sub$gene_name,
	    bin = j
	)

	hits_bin <- countOverlaps(gr_bin, gr_fragments, type = "any")
	data.table(
	    cell_type = cell_type_x,
	    transition = transition_x,
	    hits_bin = hits_bin,
	    gene = gene_ann_sub$gene_name,
	    gene_length = gene_ann_sub$end - gene_ann_sub$start,
	    bin = j
	    )
    }, mc.cores = 10) %>% rbindlist -> d_bin_hits
    d_bin_hits
}) %>% rbindlist -> d_result

write_tsv(d_result, "./tmp/table_gene_body_bin_hits_h3k27me3_B.tsv")

d_result[gene_length > 20000 & bin  == 1, .N, by = .(cell_type, transition, bin)] %>% print

d_result2 = d_result[gene_length > 20000]

d_plot = d_result2[, .(fpk = sum(hits_bin / gene_length * 1e3 * 20) / .N), by = .(cell_type, transition, bin)]

ggplot(d_plot, aes(x = bin, y = fpk, color = transition)) +
	geom_line() +
	facet_wrap(~cell_type, ncol = 1, scales = "free_y") +
	labs(x = "Gene Body Bin", y = "Mean FPK") +
	theme_classic() + ggtitle(paste0("Mean FPK Across Gene Body Bins (", antibody_target_x, ", ", sample_type_x, ")")) +
	theme(plot.title = element_text(hjust = 0.5))

ggplot(d_plot, aes(x = bin, y = fpk, color = cell_type)) +
	geom_line() +
	facet_wrap(~transition, ncol = 1, scales = "free_y") +
	labs(x = "Gene Body Bin", y = "Mean FPK") +
	theme_classic() + ggtitle(paste0("Mean FPK Across Gene Body Bins (", antibody_target_x, ", ", sample_type_x, ")")) +
	theme(plot.title = element_text(hjust = 0.5))


d_plot[, sum_fpk := sum(fpk), by = .(cell_type, transition)] %>% print
d_plot[, fpk_ratio := fpk / sum_fpk]

ggplot(d_plot, aes(x = -bin, y = fpk_ratio, color = transition)) +
	geom_line() +
	geom_point() +
	facet_wrap(~cell_type, ncol = 1) +
	labs(x = "Gene Body Bin", y = "FPK Ratio") +
	theme_classic() + ggtitle(paste0("FPK Ratio Across Gene Body Bins (", antibody_target_x, ", ", sample_type_x, ")")) +
	theme(plot.title = element_text(hjust = 0.5)) +
	coord_flip()

ggsave(paste0("figures/239_gene_body_bin_fpk_ratio_lineplot_", antibody_target_x, "_", sample_type_x, ".pdf"), width = 8, height = 10)

## }}}

## E lineage, H3K27me3, CoCnT {{{

cell_type_order = c(
    "HSPC",
    # "Pre-Pro-B cells",
    # "Immature B cells",
    # "Mature B cells1",
    # "Mature B cells2",
    # "Memory B cells",
    # "Plasma cells"
    # "Myelocyte/Classical Monocytes1",
    # "Myelocyte/Classical Monocytes2",
    # "Myelocyte/Classical Monocytes3",
    "Erythroid progenitors1",
    "Erythroid progenitors2",
    "Erythroid progenitors3"
)



antibody_target_x = "H3K27me3"
sample_type_x = "CoCnT"

sample_name_v = sample_sheet[antibody_target %in% antibody_target_x & library_type %in% sample_type_x, sample_name]

d_bed_sub = d_bed[sample_name %in% sample_name_v]

transition_v = c("whole_genome", "Active->Bivalent->Repressive", "Active->Repressive", "Active->Un->Repressive")

dm = expand.grid(
    cell_type_x = cell_type_order,
    transition_x = transition_v,
    lienage_x = c("E")
    )

# dm = dm[-which(dm$lienage == "E" & !(grepl("Ery", dm$cell_type_x))), ]
# dm = dm[-which(dm$lienage == "B" & (grepl("Ery", dm$cell_type_x))), ]


lapply(1:nrow(dm), function(i) {
    print(i)
    cell_type_x = dm$cell_type_x[i]
    transition_x = dm$transition_x[i]
    lineage_x = dm$lienage_x[i]

    marker_x = "me2"

    if (transition_x == "whole_genome") {
	gene_v = gene_body_ann$gene_name
	gene_ann_sub = gene_body_ann[gene_name %in% gene_v]
    } else {
	gene_v = gene_category[trans %in% transition_x & grepl(marker_x, marker) & lineage == lineage_x, gene]
	gene_ann_sub = gene_body_ann[gene_name %in% gene_v]
    }

    cell_id_v = d_meta[ct2 %in% cell_type_x, cell_id]
    d_bed_sub2 = d_bed_sub[cell_id %in% cell_id_v]

    gr_fragments = GRanges(
	seqnames = d_bed_sub2$chr,
	ranges = IRanges(
	    start = d_bed_sub2$start,
	    end = d_bed_sub2$end
	    ),
	cell_id = d_bed_sub2$cell_id
    )

    mclapply(1:20, function(j) {
	gr_bin = GRanges(
	    seqnames = gene_ann_sub$seq,
	    ranges = IRanges(
		start = ifelse(gene_ann_sub$strand == "+", gene_ann_sub$start + (j-1) * (gene_ann_sub$end - gene_ann_sub$start) / 20, gene_ann_sub$end - j * (gene_ann_sub$end - gene_ann_sub$start) / 20),
		end = ifelse(gene_ann_sub$strand == "+", gene_ann_sub$start + j * (gene_ann_sub$end - gene_ann_sub$start) / 20, gene_ann_sub$end - (j-1) * (gene_ann_sub$end - gene_ann_sub$start) / 20)
		),
	    strand = gene_ann_sub$strand,
	    gene_name = gene_ann_sub$gene_name,
	    bin = j
	)

	hits_bin <- countOverlaps(gr_bin, gr_fragments, type = "any")
	data.table(
	    cell_type = cell_type_x,
	    transition = transition_x,
	    hits_bin = hits_bin,
	    gene = gene_ann_sub$gene_name,
	    gene_length = gene_ann_sub$end - gene_ann_sub$start,
	    bin = j
	    )
    }, mc.cores = 10) %>% rbindlist -> d_bin_hits
    d_bin_hits
}) %>% rbindlist -> d_result

write_tsv(d_result, "./tmp/table_gene_body_bin_hits_h3k27me3_E.tsv")

d_result[gene_length > 20000 & bin  == 1, .N, by = .(cell_type, transition, bin)] %>% print

d_result2 = d_result[gene_length > 20000]

d_plot = d_result2[, .(fpk = sum(hits_bin / gene_length * 1e3 * 20) / .N), by = .(cell_type, transition, bin)]

ggplot(d_plot, aes(x = bin, y = fpk, color = transition)) +
	geom_line() +
	facet_wrap(~cell_type, ncol = 1, scales = "free_y") +
	labs(x = "Gene Body Bin", y = "Mean FPK") +
	theme_classic() + ggtitle(paste0("Mean FPK Across Gene Body Bins (", antibody_target_x, ", ", sample_type_x, ")")) +
	theme(plot.title = element_text(hjust = 0.5))

ggplot(d_plot, aes(x = bin, y = fpk, color = cell_type)) +
	geom_line() +
	facet_wrap(~transition, ncol = 1, scales = "free_y") +
	labs(x = "Gene Body Bin", y = "Mean FPK") +
	theme_classic() + ggtitle(paste0("Mean FPK Across Gene Body Bins (", antibody_target_x, ", ", sample_type_x, ")")) +
	theme(plot.title = element_text(hjust = 0.5))


d_plot[, sum_fpk := sum(fpk), by = .(cell_type, transition)] %>% print
d_plot[, fpk_ratio := fpk / sum_fpk]

ggplot(d_plot, aes(x = -bin, y = fpk_ratio, color = transition)) +
	geom_line() +
	geom_point() +
	facet_wrap(~cell_type, ncol = 1) +
	labs(x = "Gene Body Bin", y = "FPK Ratio") +
	theme_classic() + ggtitle(paste0("FPK Ratio Across Gene Body Bins (", antibody_target_x, ", ", sample_type_x, ")")) +
	theme(plot.title = element_text(hjust = 0.5)) +
	coord_flip()

ggsave(paste0("figures/239_gene_body_bin_fpk_ratio_lineplot_", antibody_target_x, "_", sample_type_x, "_E.pdf"), width = 8, height = 10)

## }}}


