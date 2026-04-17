ct_label_doublet = function(proj, doublet_percent_cuttoff = 10, redo_LSI = F) {
    if (redo_LSI) {
	proj <- addIterativeLSI(
	    ArchRProj = proj,
	    useMatrix = "TileMatrix",
	    name = "IterativeLSI",
	    force = T
	)
	proj = addHarmony(
	    ArchRProj = proj,
	    reducedDims = "IterativeLSI",
	    name = "Harmony",
	    groupBy = "library_id",
	    force = T
	)
	proj = addUMAP(
	    proj,
	    reducedDims = "Harmony",
	    name = "UMAP_Harmony",
	    nNeighbors = 30,
	    minDist = 0.5,
	    metric = "cosine",
	    force = T
	)
    }
    d_plot = data.table(data.frame(proj@cellColData))
    umap_df = proj@embeddings$UMAP_Harmony[[1]]
    colnames(umap_df) = c("UMAP1", "UMAP2")
    d_plot = cbind(d_plot, umap_df)

    # g_umap_after = ggplot(d_plot[order(status, decreasing = T)], aes(x = UMAP1, y = UMAP2, color = status)) +
	# facet_wrap(~library_id, ncol = 4) +
	# scale_color_nejm() +
	# geom_point(size = 0.5) +
	# theme_classic()

    lsi <- getReducedDims(proj, reducedDims = "Harmony")  # cells x components

    set.seed(2025)
    rs = kmeans(lsi, 100, iter.max = 30)
    d_plot$ks = rs$cluster
    d_plot = d_plot[order(status, decreasing = T),]
    ks_doublet_perc = d_plot[, 
	.(doublet_perc = sum(status == "doublet") / length(status) * 100, clone_size = .N)
	, by = ks][
	order(doublet_perc, decreasing = T)]

    g_barplot = ggplot(ks_doublet_perc, aes(x = rank(-doublet_perc, ties = "first"), y = doublet_perc)) +
	    geom_bar(stat = "identity") +
	    geom_hline(yintercept = doublet_percent_cuttoff, color = "red", linetype = "dashed") +
	    theme_classic() + labs(y = "Doublet percentage (%)", x = "Rank K-means clusters (K=100)") +
	    scale_y_continuous(sec.axis = sec_axis(~ . / 20, name = "Clone size% versus total cell")) +
	    geom_point(aes(y = clone_size / sum(clone_size) * 100 * 20), color = "blue") +
	    theme(axis.title.y.right = element_text(color = "blue"))

    include_ks = ks_doublet_perc[doublet_perc < doublet_percent_cuttoff]$ks

    g_umap_before = ggplot(d_plot, aes(x = UMAP1, y = UMAP2, color = status)) +
	scale_color_nejm() +
	geom_point(size = 0.5, alpha = 0.5) +
	theme_classic()

    g_umap_after = ggplot(d_plot[ks %in% include_ks]) +
	aes(x = UMAP1, y = UMAP2, color = status) +
	scale_color_nejm() +
    	geom_point(size = 0.5, alpha = 0.5) +
	theme_classic()

    # g_merge = ggarrange(g_umap_before, g_barplot, g_umap_after, ncol = 3, widths = c(1, 1, 1))

    # print(g_merge)
    print(g_umap_before)
    print(g_barplot)
    print(g_umap_after)

    d_plot[, doublet_label := "keep"]
    d_plot[ks %in% ks_doublet_perc[doublet_perc > doublet_percent_cuttoff, ks], doublet_label := "removed"]
}


