library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(ggsci)

# Define cell type trajectories
TRAJECTORIES <- list(
    B = c("HSPC", "Pre-Pro-B cells", "Immature B cells", "Mature B cells1", 
                "Mature B cells2", "Memory B cells", "Plasma cells"),
    E = c("HSPC", "Erythroid progenitors1", "Erythroid progenitors2", 
                "Erythroid progenitors3"),
    M = c("HSPC", "Myelocyte/Classical Monocytes1", "Myelocyte/Classical Monocytes2",
                "Myelocyte/Classical Monocytes3")
)

# Load data
d_wide_l = list(
    me1 = fread("./tmp/table_gene_chromatin_status_by_celltype_me1.tsv"),
    me2 = fread("./tmp/table_gene_chromatin_status_by_celltype_me2.tsv"),
    me3 = fread("./tmp/table_gene_chromatin_status_by_celltype_me3.tsv")
)


# Function to analyze transitions
analyze_transitions <- function(d_wide, ct_set, mark, trajectory_name) {
    d_sub <- d_wide[cell_type %in% ct_set]
    d_sub$cell_type <- factor(d_sub$cell_type, levels = ct_set)
    setorder(d_sub, gene, cell_type)

    d_sub[, `:=`(next_cat = shift(cat, type = "lead"),
	next_cell = shift(as.character(cell_type), type = "lead")), by = gene]

    transitions <- d_sub[!is.na(next_cat)]
    transitions[, step := paste0(cell_type, "->", next_cell)]

    trans_probs <- transitions[
	, .N, by = .(step, from = cat, to = next_cat)
	][, total_from := sum(N), by = .(step, from)
	][, prob := N / total_from
	][order(step, from, to)]

    trans_probs$step <- factor(trans_probs$step, 
	levels = paste0(ct_set[-length(ct_set)], "->", ct_set[-1]))

    write_tsv(trans_probs, 
	paste0("./tmp/table_chromatin_status_transition_probabilities_by_step_",
	    trajectory_name, "_", mark, ".tsv"))

    return(trans_probs)
}

# Plot
plot_transitions <- function(trans_probs, title) {
    ## remove the unchanged
    trans_probs <- trans_probs[!(from == to) & from != "Other" & to != "Other"]
    ggplot(trans_probs, aes(x = step, y = prob, group = interaction(from, to), color = to)) +
	geom_line() + geom_point() +
	scale_color_manual(values = c("Active" = "#21854e", "Repressive" = "#bc3b2a",
		"Bivalent" = "#e18727", "Un" = "#1872b5")) +
	lims(y = c(0, 0.5)) +
	facet_wrap(~from) + theme_classic() +
	    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	    ggtitle(title)
}



# Analyze each trajectory
trans_probs_by_step_E <- analyze_transitions(d_wide_l$me2, TRAJECTORIES$E, "me2", "E")
trans_probs_by_step_M <- analyze_transitions(d_wide_l$me2, TRAJECTORIES$M, "me2", "M")
trans_probs_by_step_B <- analyze_transitions(d_wide_l$me2, TRAJECTORIES$B, "me2", "B")

d = rbindlist(list(
	trans_probs_by_step_E[, trajectory := "Erythroid"],
	trans_probs_by_step_M[, trajectory := "Myeloid"],
	trans_probs_by_step_B[, trajectory := "B cell"]
))

d = d[, .(prob = max(prob)), by = .(from, to, trajectory)]

d[, category := cut(prob, breaks = c(-Inf, 0.05, 0.1, 0.2, 0.5, Inf),
	labels = c("<0.05", "0.05-0.1", "0.1-0.2", "0.2-0.5", ">0.5"))]

d_summary = d[, .(from, to, category, trajectory)][category != "<0.05"][, .(trajectory = paste(sort(unique(trajectory)), collapse = ",")),
	by = .(from, to, category)][order(from, to, category)][order(from, to)][from != "Other" & to != "Other"] 

d_summary

write_tsv(d_summary, "./tmp/table_chromatin_status_transition_probability_summary_me2.tsv")

httpgd::hgd(port = 4322)

pdf("./figures/plot_chromatin_status_transition_probabilities_by_step.pdf", width = 12, height = 8)
plot_transitions(trans_probs_by_step_E, "Chromatin status transition probabilities along Erythroid trajectory")
plot_transitions(trans_probs_by_step_B, "Chromatin status transition probabilities along B cell trajectory")
plot_transitions(trans_probs_by_step_M, "Chromatin status transition probabilities along Myeloid trajectory")
dev.off()
