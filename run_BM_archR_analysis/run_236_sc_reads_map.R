library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

## celltype order {{{
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

d_cell_id_sample = d_bed[, .(cell_id, sample_name)] %>% unique


# Helper function to prepare annotation and plot
plot_fragment <- function(target_x, library_type_x, pseudotime_x) {
    sample_v <- sample_sheet[antibody_target %in% target_x & library_type %in% library_type_x, sample_name]
    d_ann <- d_cell_id_sample[sample_name %in% sample_v]
    d_ann <- merge(
	d_ann,
	d_meta[, .(cell_id, assignment, pseudotime = get(pseudotime_x), ct2)],
	by = "cell_id"
    )
    d_ann <- d_ann[!is.na(pseudotime)]
    setorder(d_ann, ct2, pseudotime)
    d_ann[, y := -.I]
    d_plot <- merge(d_sub_r, d_ann, by = c("cell_id", "sample_name"))

    ## set the minimum width of the segment to 1% of the range
    d_plot[abs(end - start) < (gene_end - gene_begin) * 0.005, end := start + (gene_end - gene_begin) * 0.005]

    p <- ggplot(d_plot) +
	geom_segment(
	    aes(x = gene_begin, xend = gene_begin - 100, y = y, yend = y, color = ct2),
	    data = d_ann, linewidth = 1
	    ) +
	geom_segment(
	    aes(x = start, xend = end, y = y, yend = y),
	    linewidth = 0.5
	    ) +
	scale_color_manual(values = color_ct2) +
	# xlim(gene_begin, gene_end) +
	theme_classic() +
	theme(
	    axis.text.y  = element_blank(),
	    axis.ticks.y = element_blank(),
	    axis.title.y = element_blank(),
	    plot.margin = margin(5, 5, 5, 2),
	    legend.position = "none"
	    ) +
	labs(x = "Genomic position")
    p
}

## Genome region: chr9:36,824,080-37,089,361 {{{

gene_chr = "chr9"
gene_begin = 36824080
gene_end = 37089361

Cairo::CairoPDF("figures/236_fragment_plot_pax5.pdf", width = 6, height = 6)
d_sub_r = d_bed[chr == gene_chr & start >= gene_begin & end <= gene_end]
# Plot for H3K27me3 CoCnT
plot_fragment(target_x = "H3K27me3", library_type_x = "CoCnT", pseudotime_x = "Bcell_Trajectory")
# Plot for H3K4me2 CoCnT
plot_fragment(target_x = "H3K4me2", library_type_x = "CoCnT", pseudotime_x = "Bcell_Trajectory")
# Plot for H3K4me2 CoCnt1/CoCnt2
plot_fragment(target_x = "H3K4me2", library_type_x = c("CoCnt1", "CoCnt2"), pseudotime_x = "Bcell_Trajectory")
dev.off()

## }}}

## Genome region: chrX:48,781,333-48,799,150  {{{

gene_chr = "chrX"
gene_begin = 48781333
gene_end = 48799150

Cairo::CairoPDF("figures/236_fragment_plot_gata1.pdf", width = 6, height = 6)
d_sub_r = d_bed[chr == "chrX" & start >= gene_begin & end <= gene_end]
# Plot for H3K27me3 CoCnT
plot_fragment(target_x = "H3K27me3", library_type_x = "CoCnT", pseudotime_x = "Bcell_Trajectory")
# Plot for H3K4me2 CoCnT
plot_fragment(target_x = "H3K4me2", library_type_x = "CoCnT", pseudotime_x = "Bcell_Trajectory")
# Plot for H3K4me2 CoCnt1/CoCnt2
plot_fragment(target_x = "H3K4me2", library_type_x = c("CoCnt1", "CoCnt2"), pseudotime_x = "Bcell_Trajectory")
dev.off()

## }}}


## Genome region: chr4:16,039,470-16,092,386 {{{

gene_chr <- "chr4"
gene_begin <- 16039470
gene_end <- 16092386

d_sub_r <- d_bed[chr == gene_chr & start >= gene_begin & end <= gene_end]

Cairo::CairoPDF("figures/236_fragment_plot_prom1.pdf", width = 6, height = 6)
# Plot for H3K27me3 CoCnT
plot_fragment(target_x = "H3K27me3", library_type_x = "CoCnT", pseudotime_x = "Erythroid_Trajectory")

# Plot for H3K4me2 CoCnT
plot_fragment(target_x = "H3K4me2", library_type_x = "CoCnT", pseudotime_x = "Erythroid_Trajectory")

# Plot for H3K4me2 CoCnt1/CoCnt2
plot_fragment(target_x = "H3K4me2", library_type_x = c("CoCnt1", "CoCnt2"), pseudotime_x = "Erythroid_Trajectory")
dev.off()
## }}}


## Genome region: chr6:146824-188695 {{{

gene_chr <- "chr16"
gene_begin <- 146854
gene_end <- 188695

d_sub_r <- d_bed[chr == gene_chr & start >= gene_begin & end <= gene_end]

Cairo::CairoPDF("figures/236_fragment_plot_HBA1&2.pdf", width = 6, height = 6)
# Plot for H3K27me3 CoCnT
plot_fragment(target_x = "H3K27me3", library_type_x = "CoCnT", pseudotime_x = "Erythroid_Trajectory")

# Plot for H3K4me1 CoCnT
plot_fragment(target_x = "H3K4me1", library_type_x = "CoCnT", pseudotime_x = "Erythroid_Trajectory")

# Plot for H3K4me1 CoCnt1/CoCnt2
plot_fragment(target_x = "H3K4me1", library_type_x = c("CoCnt1", "CoCnt2"), pseudotime_x = "Erythroid_Trajectory")


# Plot for H3K4me2 CoCnT
plot_fragment(target_x = "H3K4me2", library_type_x = "CoCnT", pseudotime_x = "Erythroid_Trajectory")

# Plot for H3K4me2 CoCnt1/CoCnt2
plot_fragment(target_x = "H3K4me2", library_type_x = c("CoCnt1", "CoCnt2"), pseudotime_x = "Erythroid_Trajectory")
dev.off()
## }}}

## Genome region: chr20:62,437,561-62,528,319 {{{

gene_chr <- "chr20"
gene_begin <- 62437561
gene_end <- 62528319

d_sub_r <- d_bed[chr == gene_chr & start >= gene_begin & end <= gene_end]

Cairo::CairoPDF("figures/236_fragment_plot_GATA5.pdf", width = 6, height = 6)
# Plot for H3K27me3 CoCnT
plot_fragment(target_x = "H3K27me3", library_type_x = "CoCnT", pseudotime_x = "Erythroid_Trajectory")

# Plot for H3K4me1 CoCnT
plot_fragment(target_x = "H3K4me1", library_type_x = "CoCnT", pseudotime_x = "Erythroid_Trajectory")

# Plot for H3K4me1 CoCnt1/CoCnt2
plot_fragment(target_x = "H3K4me1", library_type_x = c("CoCnt1", "CoCnt2"), pseudotime_x = "Erythroid_Trajectory")


# Plot for H3K4me2 CoCnT
plot_fragment(target_x = "H3K4me2", library_type_x = "CoCnT", pseudotime_x = "Erythroid_Trajectory")

# Plot for H3K4me2 CoCnt1/CoCnt2
plot_fragment(target_x = "H3K4me2", library_type_x = c("CoCnt1", "CoCnt2"), pseudotime_x = "Erythroid_Trajectory")
dev.off()
## }}}


## Genome region: chr4:16,039,470-16,092,386 {{{

httpgd::hgd(port = 4322)

gene_chr = "chr4"
gene_begin = 16039470
gene_end = 16092386

d_sub_r = d_bed[chr == gene_chr & start >= gene_begin & end <= gene_end]

Cairo::CairoPDF("figures/236_fragment_plot_prom1.pdf", width = 6, height = 6)

sample_v = sample_sheet[antibody_target == "H3K27me3" & library_type %in% c("CoCnT"), sample_name]

d_ann = d_cell_id_sample[sample_name %in% sample_v]
d_ann = merge(d_ann, d_meta[, .(cell_id, assignment, pseudotime = Erythroid_Trajectory, ct2)], by = "cell_id")
d_ann = d_ann[!is.na(pseudotime)]
setorder(d_ann, ct2, pseudotime)
d_ann[, y := -.I]

d_plot = merge(d_sub_r, d_ann, by = c("cell_id", "sample_name"))

p <- ggplot(d_plot) +
    geom_segment(aes(x = gene_begin, xend = gene_begin - 100, y = y, yend = y, color = ct2), data = d_ann, linewidth = 1) +
    geom_segment(aes(x = start, xend = end, y = y, yend = y), linewidth = 0.5) +
    scale_color_manual(values = color_ct2) +
    theme_classic() +
    theme(
	axis.text.y  = element_blank(),
	axis.ticks.y = element_blank(),
	axis.title.y = element_blank(),
	plot.margin = margin(5, 5, 5, 2),
	legend.position = "none"
	) +
    labs(x = "Genomic position")

print(p)

sample_v = sample_sheet[antibody_target == "H3K4me2" & library_type %in% c("CoCnT"), sample_name]

d_ann = d_cell_id_sample[sample_name %in% sample_v]
d_ann = merge(d_ann, d_meta[, .(cell_id, assignment, pseudotime = Erythroid_Trajectory, ct2)], by = "cell_id")
d_ann = d_ann[!is.na(pseudotime)]
setorder(d_ann, ct2, pseudotime)
d_ann[, y := -.I]

d_plot = merge(d_sub_r, d_ann, by = c("cell_id", "sample_name"))

p <- ggplot(d_plot) +
    geom_segment(aes(x = gene_begin, xend = gene_begin - 100, y = y, yend = y, color = ct2), data = d_ann, linewidth = 1) +
    geom_segment(aes(x = start, xend = end, y = y, yend = y), linewidth = 0.5) +
    scale_color_manual(values = color_ct2) +
    theme_classic() +
    theme(
	axis.text.y  = element_blank(),
	axis.ticks.y = element_blank(),
	axis.title.y = element_blank(),
	plot.margin = margin(5, 5, 5, 2),
	legend.position = "none"
	) +
    labs(x = "Genomic position")
print(p)

sample_v = sample_sheet[antibody_target == "H3K4me2" & library_type %in% c("CoCnt1", "CoCnt2"), sample_name]

d_ann = d_cell_id_sample[sample_name %in% sample_v]
d_ann = merge(d_ann, d_meta[, .(cell_id, assignment, pseudotime = Erythroid_Trajectory, ct2)], by = "cell_id")
d_ann = d_ann[!is.na(pseudotime)]
setorder(d_ann, ct2, pseudotime)
d_ann[, y := -.I]

d_plot = merge(d_sub_r, d_ann, by = c("cell_id", "sample_name"))

p <- ggplot(d_plot) +
    geom_segment(aes(x = gene_begin, xend = gene_begin - 100, y = y, yend = y, color = ct2), data = d_ann, linewidth = 1) +
    geom_segment(aes(x = start, xend = end, y = y, yend = y), linewidth = 0.5) +
    scale_color_manual(values = color_ct2) +
    theme_classic() +
    theme(
	axis.text.y  = element_blank(),
	axis.ticks.y = element_blank(),
	axis.title.y = element_blank(),
	plot.margin = margin(5, 5, 5, 2),
	legend.position = "none"
	) +
    labs(x = "Genomic position")
print(p)

dev.off()
## }}}


## back up {{{
## Count Y reads per cell
d_sub_Y = d_bed[chr == "chrY"]
d_sub_Y
merge(d_meta[, .(cell_id, assignment)], d_sub_Y[, .N, cell_id])[
    , median(N), assignment]

## Count X reads per cell
d_sub_X = d_bed[chr == "chrX"]
d_sub_X
merge(d_meta[, .(cell_id, assignment)], d_sub_X[, .N, cell_id])[
    , median(N), assignment]


## H3K27me3 distribution on Y chrom
httpgd::hgd(port = 4323)

sample_v = sample_sheet[antibody_target == "H3K27me3" & library_type == "CoCnT", sample_name]
d_sub_Y_k27 = d_sub_Y[sample_name %in% sample_v]
d_sub_Y_k27 = merge(d_sub_Y_k27, d_meta[, .(cell_id, assignment)], by = "cell_id")
d_sub_Y_k27 = d_sub_Y_k27[assignment %in% c("0", "1", "2", "3")][order(assignment)]

## order by donnor

ggplot(d_sub_Y_k27) +
    geom_segment(aes(x = start, xend = end,
	    y = cell_id, yend = cell_id),
	linewidth = 0.5) +
    theme_classic() +
    theme(
	axis.text.y  = element_blank(),
	axis.ticks.y = element_blank(),
	axis.title.y = element_blank()
	) +
    labs(x = "Genomic position")

## H3K27me3 distribution on X chrom
sample_v = sample_sheet[antibody_target == "H3K4me2" & library_type %in% c("CoCnt1", "CoCnt2"), sample_name]
# sample_v = sample_sheet[antibody_target == "H3K4me2" & library_type %in% c("CoCnT"), sample_name]
d_sub_X = d_bed[chr == "chr4"]
d_sub_X_k27 = d_sub_X[sample_name %in% sample_v]
d_sub_X_k27 = merge(d_sub_X_k27, d_meta[, .(cell_id, assignment, Bcell_Trajectory)], by = "cell_id") %>% na.omit
# d_sub_X_k27 = d_sub_X_k27[assignment %in% c("0", "1", "2", "3")][order(assignment)] %>% na.omit

unique(d_sub_X_k27$cell_id) %>% length

## order by donnor

g = ggplot(d_sub_X_k27) +
    geom_segment(aes(x = start, xend = end,
	    y = -Bcell_Trajectory, yend = -Bcell_Trajectory),
	linewidth = 0.5) +
    # facet_wrap(~ assignment, ncol = 1) +
    theme_classic() +
    theme(
	axis.text.y  = element_blank(),
	axis.ticks.y = element_blank(),
	axis.title.y = element_blank()
	) +
    labs(x = "Genomic position")

ggsave("figures/test_x_H3K4me2_cooc.pdf", g, width = 4, height = 3)

## }}}
