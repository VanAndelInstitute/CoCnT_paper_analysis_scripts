library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

httpgd::hgd(port = 4323)

d_meta <- fread("./tmp/table_metadata_h3k27me3_final3_with_trajectories.tsv")

m_list <- readRDS("./tmp/list_matrix_imputed_gene_score.rds")


## smooth the gene expression and find the peak

x = d_meta[, .(cell_id, Bcell_Trajectory, ct3, nFrags)] %>% na.omit


y = m_list$H3K4me3[c("RAG1", "RAG2", "DNTT", "TCF3", "EBF1", "PAX5", "IL7R", "ID1", "ID2", "ID3"), ]
cell_id = colnames(y)
y = as.data.table(as.matrix(t(y)))
y = data.table(cell_id = cell_id, y)

d = merge(x, y, by = "cell_id")
d = d[order(Bcell_Trajectory)]

## RAG1 {{{
d_plot = d[, .(x = Bcell_Trajectory, y = d[["RAG1"]])]

fit <- loess(y ~ x, data = d_plot, span = 0.05)

xg <- seq(min(d_plot$x), max(d_plot$x), length.out = 2000)
yg <- predict(fit, newdata = data.frame(x = xg))
x_peak <- xg[which.max(yg)]
y_peak <- max(yg, na.rm = TRUE)
x_peak

xg2 <- seq(min(d_plot$x), max(10), length.out = 2000)
yg2 <- predict(fit, newdata = data.frame(x = xg2))
x_peak2 <- xg2[which.max(yg2)]
y_peak2 <- max(yg2, na.rm = TRUE)

xg3 <- seq(x_peak2 + 1, x_peak - 1, length.out = 2000)
yg3 <- predict(fit, newdata = data.frame(x = xg3))
x_peak3 <- xg3[which.min(yg3)]
y_peak3 <- min(yg3, na.rm = TRUE)


## show the fitted curve
ggplot(d_plot, aes(x = x, y = y)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_line(data = data.frame(x = xg, y = yg), aes(x = x, y = y), color = "red", size = 1) +
    geom_vline(xintercept = x_peak, col = "blue", lty = 2) +
    geom_vline(xintercept = x_peak2, col = "green", lty = 2) +
    geom_vline(xintercept = x_peak3, col = "purple", lty = 2) +
    annotate("text", x = x_peak + 20, y = y_peak, label = paste0("Peak: ", round(x_peak, 2)), vjust = -1, color = "blue") +
    annotate("text", x = x_peak2 + 20, y = y_peak2, label = paste0("Peak2: ", round(x_peak2, 2)), vjust = -1, color = "green") +
    annotate("text", x = x_peak3 + 20, y = y_peak3, label = paste0("Trough: ", round(x_peak3, 2)), vjust = -1, color = "purple") +
    theme_classic() 

ggsave("./figures/242_rag1_trajectry_peak.pdf", width = 6, height = 4)

### }}}

## TdT {{{
d_plot = d[, .(x = Bcell_Trajectory, y = d[["DNTT"]])]

fit <- loess(y ~ x, data = d_plot, span = 0.05)

xg <- seq(min(d_plot$x), max(d_plot$x), length.out = 2000)
yg <- predict(fit, newdata = data.frame(x = xg))
x_peak <- xg[which.max(yg)]
y_peak <- max(yg, na.rm = TRUE)
x_peak


## show the fitted curve
ggplot(d_plot, aes(x = x, y = y)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_line(data = data.frame(x = xg, y = yg), aes(x = x, y = y), color = "red", size = 1) +
    geom_vline(xintercept = x_peak, col = "blue", lty = 2) +
    annotate("text", x = x_peak + 20, y = y_peak, label = paste0("Peak: ", round(x_peak, 2)), vjust = -1) +
    theme_classic()

ggsave("./figures/242_tdt_trajectry_peak.pdf", width = 6, height = 4)
## }}}

## the curve for every gene and put them together {{{

pdf("./figures/242_B_lienage_gene_trajectory.pdf", width = 6, height = 4)

names(m_list)
plot_trajectory <- function(mark, lineage, genes, plot_cols) {
  y <- m_list[[mark]][genes, ]
  y <- as.data.table(as.matrix(t(y)))
  y[, cell_id := colnames(m_list[[mark]])]
  setcolorder(y, c("cell_id", genes))
  d <- merge(x, y, by = "cell_id")
  d <- d[order(Bcell_Trajectory)]
  d_plot <- melt(d[, c("Bcell_Trajectory", plot_cols), with = FALSE],
                 id.vars = "Bcell_Trajectory", variable.name = "gene", value.name = "gene_score")
  ggplot(d_plot, aes(x = Bcell_Trajectory, y = gene_score, color = gene)) +
    geom_smooth(method = "loess", se = FALSE, span = 0.15, n = 2000) +
    theme_classic()
}

genes <- c("RAG1", "DNTT", "EBF1", "PAX5")

print(plot_trajectory("H3K27me3", genes, c("PAX5", "EBF1", "RAG1", "DNTT")))
print(plot_trajectory("H3K4me3", genes, c("PAX5", "EBF1", "RAG1", "DNTT")))
print(plot_trajectory("H3K4me2", genes, c("PAX5", "EBF1", "RAG1", "DNTT")))
print(plot_trajectory("H3K4me1", genes, c("PAX5", "EBF1", "RAG1", "DNTT")))
print(plot_trajectory("H3K4me2_cooc", genes, c("PAX5", "EBF1", "RAG1", "DNTT")))

dev.off()

pdf("./figures/242_B_lienage_gene_trajectory.pdf", width = 6, height = 4)

x = d_meta[, .(cell_id, Bcell_Trajectory, ct3, nFrags)] %>% na.omit
names(m_list)
plot_trajectory <- function(mark, genes, plot_cols) {
  y <- m_list[[mark]][genes, ]
  y <- as.data.table(as.matrix(t(y)))
  y[, cell_id := colnames(m_list[[mark]])]
  setcolorder(y, c("cell_id", genes))
  d <- merge(x, y, by = "cell_id")
  d <- d[order(Bcell_Trajectory)]
  d_plot <- melt(d[, c("Bcell_Trajectory", plot_cols), with = FALSE],
                 id.vars = "Bcell_Trajectory", variable.name = "gene", value.name = "gene_score")
  ggplot(d_plot, aes(x = Bcell_Trajectory, y = gene_score, color = gene)) +
    geom_smooth(method = "loess", se = FALSE, span = 0.15, n = 2000) +
    theme_classic()
}

genes <- c("RAG1", "DNTT", "EBF1", "PAX5")

print(plot_trajectory("H3K27me3", genes, c("PAX5", "EBF1", "RAG1", "DNTT")))
print(plot_trajectory("H3K4me3", genes, c("PAX5", "EBF1", "RAG1", "DNTT")))
print(plot_trajectory("H3K4me2", genes, c("PAX5", "EBF1", "RAG1", "DNTT")))
print(plot_trajectory("H3K4me1", genes, c("PAX5", "EBF1", "RAG1", "DNTT")))
print(plot_trajectory("H3K4me2_cooc", genes, c("PAX5", "EBF1", "RAG1", "DNTT")))
print(plot_trajectory("H3K4me2_cooc", genes, c("PAX5", "EBF1", "RAG1", "DNTT")))

dev.off()

pdf("./figures/242_E_lienage_gene_trajectory.pdf", width = 6, height = 4)

x = d_meta[, .(cell_id, Erythroid_Trajectory, ct3, nFrags)] %>% na.omit

names(m_list)
plot_trajectory <- function(mark, genes, plot_cols) {
  y <- m_list[[mark]][genes, ]
  y <- as.data.table(as.matrix(t(y)))
  y[, cell_id := colnames(m_list[[mark]])]
  setcolorder(y, c("cell_id", genes))
  d <- merge(x, y, by = "cell_id")
  d <- d[order(Erythroid_Trajectory)]
  d_plot <- melt(d[, c("Erythroid_Trajectory", plot_cols), with = FALSE],
                 id.vars = "Erythroid_Trajectory", variable.name = "gene", value.name = "gene_score")
  ggplot(d_plot, aes(x = Erythroid_Trajectory, y = gene_score, color = gene)) +
    geom_smooth(method = "loess", se = FALSE, span = 0.15, n = 2000) +
    theme_classic()
}

genes = c("GATA1", "HBA1", "HBA2", "KLF1")
print(plot_trajectory("H3K27me3", genes, genes))
print(plot_trajectory("H3K4me3", genes, genes))
print(plot_trajectory("H3K4me2", genes, genes))
print(plot_trajectory("H3K4me1", genes, genes))
print(plot_trajectory("H3K4me2_cooc", genes, genes))
print(plot_trajectory("H3K4me2_cooc", genes, genes))

dev.off()



# y = m_list$H3K4me2[c("RAG1", "RAG2", "DNTT", "TCF3", "EBF1", "PAX5", "IL7R"), ]
# cell_id = colnames(y)
# y = as.data.table(as.matrix(t(y)))
# y = data.table(cell_id = cell_id, y)
#
# d = merge(x, y, by = "cell_id")
# d = d[order(Bcell_Trajectory)]
#
#
# d_plot = melt(d[, .(Bcell_Trajectory, TCF3, PAX5, EBF1)], id.vars = "Bcell_Trajectory", variable.name = "gene", value.name = "expression")
# ggplot(d_plot, aes(x = Bcell_Trajectory, y = expression, color = gene)) +
# 	geom_smooth(method = "loess", se = FALSE, span = 0.05, n = 2000) +
# 	geom_vline(xintercept = 9.2, col = "blue", lty = 2) +
# 	geom_vline(xintercept = 17, col = "green", lty = 2) +
# 	geom_vline(xintercept = 12.26, col = "purple", lty = 2) +
# 	theme_classic()
#
# y = m_list$H3K4me1[c("RAG1", "RAG2", "DNTT", "TCF3", "EBF1", "PAX5", "IL7R"), ]
# cell_id = colnames(y)
# y = as.data.table(as.matrix(t(y)))
# y = data.table(cell_id = cell_id, y)
#
# d = merge(x, y, by = "cell_id")
# d = d[order(Bcell_Trajectory)]
#
#
# d_plot = melt(d[, .(Bcell_Trajectory, TCF3, PAX5, EBF1)], id.vars = "Bcell_Trajectory", variable.name = "gene", value.name = "expression")
# ggplot(d_plot, aes(x = Bcell_Trajectory, y = expression, color = gene)) +
# 	geom_smooth(method = "loess", se = FALSE, span = 0.05, n = 2000) +
# 	geom_vline(xintercept = 9.2, col = "blue", lty = 2) +
# 	geom_vline(xintercept = 17, col = "green", lty = 2) +
# 	geom_vline(xintercept = 12.26, col = "purple", lty = 2) +
# 	theme_classic()
#
# y = m_list$H3K27me3[c("RAG1", "RAG2", "DNTT", "TCF3", "EBF1", "PAX5", "IL7R"), ]
# cell_id = colnames(y)
# y = as.data.table(as.matrix(t(y)))
# y = data.table(cell_id = cell_id, y)
#
# d = merge(x, y, by = "cell_id")
# d = d[order(Bcell_Trajectory)]
#
# d_plot = melt(d[, .(Bcell_Trajectory, TCF3, PAX5, EBF1)], id.vars = "Bcell_Trajectory", variable.name = "gene", value.name = "expression")
# ggplot(d_plot, aes(x = Bcell_Trajectory, y = expression, color = gene)) +
# 	geom_smooth(method = "loess", se = FALSE, span = 0.05, n = 2000) +
# 	geom_vline(xintercept = 9.2, col = "blue", lty = 2) +
# 	geom_vline(xintercept = 17, col = "green", lty = 2) +
# 	geom_vline(xintercept = 12.26, col = "purple", lty = 2) +
# 	theme_classic()
#
# y = m_list$H3K27me3_cooc[c("RAG1", "RAG2", "DNTT", "TCF3", "EBF1", "PAX5", "IL7R"), ]
# cell_id = colnames(y)
# y = as.data.table(as.matrix(t(y)))
# y = data.table(cell_id = cell_id, y)
#
# d = merge(x, y, by = "cell_id")
# d = d[order(Bcell_Trajectory)]
#
# d_plot = melt(d[, .(Bcell_Trajectory, TCF3, PAX5, EBF1)], id.vars = "Bcell_Trajectory", variable.name = "gene", value.name = "expression")
# ggplot(d_plot, aes(x = Bcell_Trajectory, y = expression, color = gene)) +
# 	geom_smooth(method = "loess", se = FALSE, span = 0.05, n = 2000) +
# 	geom_vline(xintercept = 9.2, col = "blue", lty = 2) +
# 	geom_vline(xintercept = 17, col = "green", lty = 2) +
# 	geom_vline(xintercept = 12.26, col = "purple", lty = 2) +
# 	theme_classic()
#
#
#
# ggplot(d, aes(x = Bcell_Trajectory, y = nFrags, color = ct3)) +
#     geom_point(size = 0.5, alpha = 0.3) +
#     # geom_smooth(method = "loess", se = FALSE, span = 0.05, n = 2000) +
#     geom_vline(xintercept = 9.2, col = "blue", lty = 2) +
#     geom_vline(xintercept = 17, col = "green", lty = 2) +
#     geom_vline(xintercept = 12.26, col = "purple", lty = 2) +
#     theme_classic() +
#     scale_y_log10()
#
#
#
#
# # points(x_peak, y_peak, col = "blue", pch = 16, cex = 1.5)
#
# # with(d_plot, points(Bcell_Trajectory, TCF3, type = "p"))
#
# xg <- seq(min(df$x), max(df$x), length.out = 2000)
# yg <- predict(fit, newdata = data.frame(x = xg))
#
#
#
# with(d_plot, plot(Bcell_Trajectory, TCF3, type = "p"))
# with(d_plot, points(Bcell_Trajectory, ksmooth(TCF3, Bcell_Trajectory, kernel = "normal", bandwidth = 0.5)$y, type = "l"))
#
# ggplot(d_plot, aes(x = Bcell_Trajectory, y = TCF3)) +
#     geom_point() +
#     geom_smooth(method = "loess", se = FALSE, n = 100) +
#     theme_classic()
#
## }}}

# gat the peak for each donor
x = d_meta[, .(cell_id, Bcell_Trajectory, ct3, nFrags, assignment)] %>% na.omit

# y = m_list$H3K27me3[c("RAG1", "RAG2", "DNTT", "TCF3", "EBF1", "PAX5", "IL7R", "ID1", "ID2", "ID3"), ]
y = m_list$H3K27me3[c("RAG1", "RAG2", "DNTT", "TCF3", "EBF1", "PAX5", "IL7R", "ID1", "ID2", "ID3"), ]
cell_id = colnames(y)
y = as.data.table(as.matrix(t(y)))
y = data.table(cell_id = cell_id, y)

d = merge(x, y, by = "cell_id")
d = d[order(Bcell_Trajectory)]

d = d[assignment %in% c(0, 1, 2, 3)]


ggplot(d) + aes(x = Bcell_Trajectory, y = DNTT) +
	# geom_point(size = 0.5, alpha = 0.5) +
	geom_smooth(method = "loess", se = FALSE, span = 0.10, n = 2000) +
	theme_classic()
