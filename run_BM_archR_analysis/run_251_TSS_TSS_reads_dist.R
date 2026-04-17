library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(parallel)
library(httpgd)

httpgd::hgd(port = 4322)

# =========================================================
# 1. Gene annotation
# =========================================================
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
gr_gene <- genes(txdb)

gene_body_ann <- as.data.frame(gr_gene)[, c("seqnames", "start", "end", "strand", "gene_id")] %>%
    data.table()

gene_body_ann[, gene_name := mapIds(
    org.Hs.eg.db,
    keys = gene_id,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
)]

gene_body_ann <- gene_body_ann[
    seqnames %in% paste0("chr", c(1:22, "X", "Y")) &
        !is.na(gene_name)
]

# gene length
gene_body_ann[, gene_length := end - start + 1L]

# =========================================================
# 2. Metadata and input tables
# =========================================================
cell_type_order <- c(
    "HSPC",
    "Pre-Pro-B cells",
    "Immature B cells",
    "Mature B cells1",
    "Mature B cells2",
    "Memory B cells",
    "Plasma cells"
)

d_meta <- fread("./tmp/table_metadata_h3k27me3_final3_with_trajectories.tsv")
d_meta[, ct2 := factor(ct2, levels = cell_type_order)]

sample_sheet <- rbind(
    fread("../../data/sample_sheet_PR001798.tsv"),
    fread("../../data/sample_sheet_PR001799.tsv"),
    fread("../../data/sample_sheet_PR001855.tsv"),
    fread("../../data/sample_sheet_PR001856.tsv")
)

d_bed <- fread("../run_BM_bed_profile/tmp/table_good_cells_fragments.tsv")

d_bed_clean = d_bed[end > start]

gene_category <- fread("./tmp/table_shared_genes_between_markers_in_trajectory.tsv")
gene_category_bi <- gene_category[grepl("\\+", marker)]

## merge (Un->Active and Un->Active->Un into "Un->Active", Bivalent->Active and Bivalent->Active->Bivalent into "Bivalent->Active")
gene_category_bi[trans == "Un->Active->Un", trans := "Un->Active"]
gene_category_bi[trans == "Bivalent->Active->Bivalent", trans := "Bivalent->Active"]


# =========================================================
# 3. Helper: build metagene bins
# =========================================================
make_meta_bins <- function(gene_dt, n_body = 20, flank_bp = 10000, flank_bins = 20) {
    res_list <- list()
    up_bin_size <- flank_bp / flank_bins

    make_iranges_safe <- function(start_vec, end_vec) {
        s <- as.integer(floor(pmin(start_vec, end_vec)))
        e <- as.integer(ceiling(pmax(start_vec, end_vec)))
        s[s < 1L] <- 1L
        e[e < 1L] <- 1L
        IRanges(start = s, end = e)
    }

    # ---------- upstream ----------
    for (j in seq_len(flank_bins)) {
        start_vec <- ifelse(
            gene_dt$strand == "+",
            gene_dt$start - flank_bp + (j - 1) * up_bin_size,
            gene_dt$end + (flank_bins - j) * up_bin_size + 1
        )
        end_vec <- ifelse(
            gene_dt$strand == "+",
            gene_dt$start - flank_bp + j * up_bin_size - 1,
            gene_dt$end + (flank_bins - j + 1) * up_bin_size
        )

        gr_up <- GRanges(
            seqnames = gene_dt$seqnames,
            ranges = make_iranges_safe(start_vec, end_vec),
            strand = gene_dt$strand,
            gene_name = gene_dt$gene_name,
            region = "upstream",
            bin = j,
            bin_label = paste0("U", j)
        )
        res_list[[length(res_list) + 1]] <- gr_up
    }

    # ---------- gene body ----------
    for (j in seq_len(n_body)) {
        start_vec <- ifelse(
            gene_dt$strand == "+",
            gene_dt$start + (j - 1) * gene_dt$gene_length / n_body,
            gene_dt$end - j * gene_dt$gene_length / n_body + 1
        )
        end_vec <- ifelse(
            gene_dt$strand == "+",
            gene_dt$start + j * gene_dt$gene_length / n_body - 1,
            gene_dt$end - (j - 1) * gene_dt$gene_length / n_body
        )

        gr_body <- GRanges(
            seqnames = gene_dt$seqnames,
            ranges = make_iranges_safe(start_vec, end_vec),
            strand = gene_dt$strand,
            gene_name = gene_dt$gene_name,
            region = "body",
            bin = j,
            bin_label = paste0("B", j)
        )
        res_list[[length(res_list) + 1]] <- gr_body
    }

    # ---------- downstream ----------
    for (j in seq_len(flank_bins)) {
        start_vec <- ifelse(
            gene_dt$strand == "+",
            gene_dt$end + (j - 1) * up_bin_size + 1,
            gene_dt$start - j * up_bin_size
        )
        end_vec <- ifelse(
            gene_dt$strand == "+",
            gene_dt$end + j * up_bin_size,
            gene_dt$start - (j - 1) * up_bin_size - 1
        )

        gr_down <- GRanges(
            seqnames = gene_dt$seqnames,
            ranges = make_iranges_safe(start_vec, end_vec),
            strand = gene_dt$strand,
            gene_name = gene_dt$gene_name,
            region = "downstream",
            bin = j,
            bin_label = paste0("D", j)
        )
        res_list[[length(res_list) + 1]] <- gr_down
    }

    gr_all <- do.call(c, res_list)

    bin_levels <- c(
        paste0("U", seq_len(flank_bins)),
        paste0("B", seq_len(n_body)),
        paste0("D", seq_len(flank_bins))
    )
    gr_all$bin_global <- factor(gr_all$bin_label, levels = bin_levels)

    gr_all
}


# =========================================================
# 4. Main function to calculate metagene profile
# =========================================================
calc_metagene_profile <- function(
  d_bed,
  d_meta,
  sample_sheet,
  gene_body_ann,
  gene_category,
  cell_type_order,
  antibody_target_x = "H3K4me2",
  sample_type_x = "CoCnT",
  lineage_x = "B",
  # marker_x = "me2",
  # transition_v = c("whole_genome", "Un->Active", "Un->Active->Un", "Bivalent->Active", "Bivalent->Active->Bivalent"),
  transition_v = c("Un->Active", "Bivalent->Active"),
  min_gene_length = 20000,
  n_body = 20,
  flank_bp = 10000,
  flank_bins = 2,
  mc.cores = 10
) {
    sample_name_v <- sample_sheet[
        antibody_target %in% antibody_target_x & library_type %in% sample_type_x,
        sample_name
    ]

    d_bed_sub <- d_bed[sample_name %in% sample_name_v]

    dm <- expand.grid(
        cell_type_x = cell_type_order,
        transition_x = transition_v,
        lineage_x = lineage_x,
        stringsAsFactors = FALSE
    )

    res_all <- lapply(seq_len(nrow(dm)), function(i) {
        message("Processing ", i, "/", nrow(dm))

        cell_type_x <- dm$cell_type_x[i]
        transition_x <- dm$transition_x[i]
        lineage_x_i <- dm$lineage_x[i]

	print(paste0("cell_type_x: ", cell_type_x, "; transition_x: ", transition_x, "; lineage_x_i: ", lineage_x_i))

        if (transition_x == "whole_genome") {
            gene_ann_sub <- gene_body_ann[gene_length > min_gene_length]
        } else {
            gene_v <- gene_category[
                # trans %in% transition_x & grepl(marker_x, marker) & lineage == lineage_x_i,
                trans %in% transition_x & lineage == lineage_x_i,
                gene
            ] %>% unique()

            gene_ann_sub <- gene_body_ann[
                gene_name %in% gene_v & gene_length > min_gene_length
            ]
        }

        if (nrow(gene_ann_sub) == 0) {
            return(NULL)
        }

        cell_id_v <- d_meta[ct2 %in% cell_type_x, cell_id]
        d_bed_sub2 <- d_bed_sub[cell_id %in% cell_id_v]
	d_bed_sub2 <- d_bed_sub2[chr %in% paste0("chr", c(1:22, "X", "Y"))]


        if (nrow(d_bed_sub2) == 0) {
            return(NULL)
        }

        n_cells <- uniqueN(d_bed_sub2$cell_id)

        gr_fragments <- GRanges(
            seqnames = d_bed_sub2$chr,
            ranges = IRanges(start = d_bed_sub2$start, end = d_bed_sub2$end),
            cell_id = d_bed_sub2$cell_id
        )

        gr_bins <- make_meta_bins(
            gene_dt = gene_ann_sub,
            n_body = n_body,
            flank_bp = flank_bp,
            flank_bins = flank_bins
        )

        hits_bin <- countOverlaps(gr_bins, gr_fragments, type = "any")

        d_out <- data.table(
            cell_type = cell_type_x,
            transition = transition_x,
            lineage = lineage_x_i,
            gene = gr_bins$gene_name,
            region = gr_bins$region,
            bin = as.integer(gr_bins$bin),
            bin_label = as.character(gr_bins$bin_label),
            bin_global = as.character(gr_bins$bin_global),
            hits_bin = hits_bin,
            n_cells = n_cells
        )

        # per-bin width for normalization
        d_out[, bin_width := fifelse(
            region == "body",
            gene_ann_sub[match(gene, gene_name), gene_length] / n_body,
            flank_bp / flank_bins
        )]

	d_out[, signal_per_kb_1000cells := hits_bin / (bin_width / 1000) / n_cells * 1000]

        d_out
    }) %>% rbindlist(fill = TRUE)

    # order bins
    bin_levels <- c(
        paste0("U", seq_len(flank_bins)),
        paste0("B", seq_len(n_body)),
        paste0("D", seq_len(flank_bins))
    )
    d_out <- copy(res_all)
    d_out[, bin_global := factor(bin_global, levels = bin_levels)]

    # normalize as mean fragments per kb per gene
    d_plot <- d_out[
        ,
        .(signal = mean(signal_per_kb_1000cells, na.rm = TRUE)),
        by = .(cell_type, transition, lineage, bin_global)
    ]

    d_plot[, bin_index := as.integer(bin_global)]

    list(raw = d_out, plot = d_plot)
}

# =========================================================
# 5. Run for B lineage
# =========================================================

## H3K4me3: Un->Active and Bivalent->Active transitions {{{
res_B <- calc_metagene_profile(
    d_bed = d_bed_clean,
    d_meta = d_meta,
    sample_sheet = sample_sheet,
    gene_body_ann = gene_body_ann,
    gene_category = gene_category_bi,
    cell_type_order = cell_type_order[-3],
    transition_v = c("Un->Active", "Bivalent->Active"),
    antibody_target_x = "H3K4me3",
    sample_type_x = "CoCnT",
    lineage_x = "B",
    # marker_x = "me2",
    min_gene_length = 5000,
    n_body = 20,
    flank_bp = 10000,
    flank_bins = 20,
    mc.cores = 10
)

d_plot_B <- res_B$plot

d_plot_B <- d_plot_B[transition %in% c("Un->Active", "Bivalent->Active")]

# custom x-axis labels
x_breaks <- c(1, 10, 20, 30, 40, 50, 60)
x_labels <- c("-10 kb", "-5 kb", "TSS (0%)", "50%", "TTS (100%)", "+5 kb", "+10 kb")

d_plot_B$bin_global %>% table()

d_plot_B$cell_type %<>% factor(., levels = cell_type_order)

p_B <- ggplot(d_plot_B, aes(x = bin_index, y = signal, color = transition, group = transition)) +
    geom_line(linewidth = 1) +
    facet_wrap(~cell_type, ncol = 1) +
    scale_x_continuous(
        breaks = x_breaks,
        labels = x_labels,
        expand = c(0.01, 0.01)
    ) +
    labs(
        x = "Gene coordinate",
        y = "Fragment density (per kb per 1000 cells)",
        title = "Metagene profile across gene body and flanking regions (B lineage)"
    ) +
    scale_color_manual(values = c("Un->Active" = "#1972b5", "Bivalent->Active" = "#e18726")) +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    geom_vline(xintercept = c(20, 40), linetype = 2, color = "grey50")

print(p_B)

ggsave(
    "figures/251_metagene_profile_H3K27me3_CoCnT_B.pdf",
    p_B,
    width = 8, height = 10
)
## }}}

## H3K4me2: Un->Active and Bivalent->Active transitions {{{
res_B <- calc_metagene_profile(
    d_bed = d_bed,
    d_meta = d_meta,
    sample_sheet = sample_sheet,
    gene_body_ann = gene_body_ann,
    gene_category = gene_category_bi,
    cell_type_order = cell_type_order,
    antibody_target_x = "H3K4me2",
    sample_type_x = "CoCnT",
    lineage_x = "B",
    # marker_x = "me2",
    min_gene_length = 5000,
    n_body = 20,
    flank_bp = 10000,
    flank_bins = 20,
    mc.cores = 10
)

d_plot_B <- res_B$plot

d_plot_B <- d_plot_B[transition %in% c("Un->Active", "Bivalent->Active")]

# custom x-axis labels
x_breaks <- c(1, 10, 20, 30, 40, 50, 60)
x_labels <- c("-10 kb", "-5 kb", "TSS (0%)", "50%", "TTS (100%)", "+5 kb", "+10 kb")


d_plot_B$bin_global %>% table()

d_plot_B$cell_type %<>% factor(., levels = cell_type_order)

p_B <- ggplot(d_plot_B, aes(x = bin_index, y = signal, color = transition, group = transition)) +
    geom_line(linewidth = 1) +
    facet_wrap(~cell_type, ncol = 1) +
    scale_x_continuous(
        breaks = x_breaks,
        labels = x_labels,
        expand = c(0.01, 0.01)
    ) +
    labs(
        x = "Gene coordinate",
        y = "Fragment density (per kb per 1000 cells)",
        title = "Metagene profile across gene body and flanking regions (B lineage)"
    ) +
    scale_color_manual(values = c("Un->Active" = "#1972b5", "Bivalent->Active" = "#e18726")) +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    geom_vline(xintercept = c(20, 40), linetype = 2, color = "grey50")

print(p_B)

ggsave(
    "figures/251_metagene_profile_H3K27me2_CoCnT_B.pdf",
    p_B,
    width = 8, height = 10
)
## }}}

## H3K4me1: Un->Active and Bivalent->Active transitions {{{
res_B <- calc_metagene_profile(
    d_bed = d_bed,
    d_meta = d_meta,
    sample_sheet = sample_sheet,
    gene_body_ann = gene_body_ann,
    gene_category = gene_category_bi,
    cell_type_order = cell_type_order,
    antibody_target_x = "H3K4me1",
    sample_type_x = "CoCnT",
    lineage_x = "B",
    # marker_x = "me2",
    min_gene_length = 5000,
    n_body = 20,
    flank_bp = 10000,
    flank_bins = 20,
    mc.cores = 10
)

d_plot_B <- res_B$plot

d_plot_B <- d_plot_B[transition %in% c("Un->Active", "Bivalent->Active")]

d_plot_B$cell_type %<>% factor(., levels = cell_type_order)

# custom x-axis labels
x_breaks <- c(1, 10, 20, 30, 40, 50, 60)
x_labels <- c("-10 kb", "-5 kb", "TSS (0%)", "50%", "TTS (100%)", "+5 kb", "+10 kb")

d_plot_B$bin_global %>% table()

p_B <- ggplot(d_plot_B, aes(x = bin_index, y = signal, color = transition, group = transition)) +
    geom_line(linewidth = 1) +
    facet_wrap(~cell_type, ncol = 1) +
    scale_x_continuous(
        breaks = x_breaks,
        labels = x_labels,
        expand = c(0.01, 0.01)
    ) +
    labs(
        x = "Gene coordinate",
        y = "Fragment density (per kb per 1000 cells)",
        title = "Metagene profile across gene body and flanking regions (B lineage)"
    ) +
    scale_color_manual(values = c("Un->Active" = "#1972b5", "Bivalent->Active" = "#e18726")) +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    geom_vline(xintercept = c(20, 40), linetype = 2, color = "grey50")

print(p_B)

ggsave(
    "figures/251_metagene_profile_H3K27me1_CoCnT_B.pdf",
    p_B,
    width = 8, height = 10
)
## }}}


## H3K4me1: Un->Active and Bivalent->Active transitions {{{
res_B <- calc_metagene_profile(
    d_bed = d_bed,
    d_meta = d_meta,
    sample_sheet = sample_sheet,
    gene_body_ann = gene_body_ann,
    gene_category = gene_category_bi,
    cell_type_order = cell_type_order,
    antibody_target_x = "H3K27me3",
    sample_type_x = "CoCnT",
    lineage_x = "B",
    # marker_x = "me2",
    min_gene_length = 5000,
    n_body = 20,
    flank_bp = 10000,
    flank_bins = 20,
    mc.cores = 10
)

d_plot_B <- res_B$plot

d_plot_B <- d_plot_B[transition %in% c("Un->Active", "Bivalent->Active")]

d_plot_B$cell_type %<>% factor(., levels = cell_type_order)

# custom x-axis labels
x_breaks <- c(1, 10, 20, 30, 40, 50, 60)
x_labels <- c("-10 kb", "-5 kb", "TSS (0%)", "50%", "TTS (100%)", "+5 kb", "+10 kb")

d_plot_B$bin_global %>% table()

p_B <- ggplot(d_plot_B, aes(x = bin_index, y = signal, color = transition, group = transition)) +
    geom_line(linewidth = 1) +
    facet_wrap(~cell_type, ncol = 1) +
    scale_x_continuous(
        breaks = x_breaks,
        labels = x_labels,
        expand = c(0.01, 0.01)
    ) +
    labs(
        x = "Gene coordinate",
        y = "Fragment density (per kb per 1000 cells)",
        title = "Metagene profile across gene body and flanking regions (B lineage)"
    ) +
    scale_color_manual(values = c("Un->Active" = "#1972b5", "Bivalent->Active" = "#e18726")) +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    geom_vline(xintercept = c(20, 40), linetype = 2, color = "grey50")

print(p_B)

ggsave(
    "figures/251_metagene_profile_H3K27me3_CoCnT_B.pdf",
    p_B,
    width = 8, height = 10
)
## }}}


