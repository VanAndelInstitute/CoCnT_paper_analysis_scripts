run_cluster = function(proj, minDist = 0.5, nNeighbors = 30, seed = 123, harmony_corcutoff = 0.75, cluster_corcutoff = 0.75, cluster_resolusion = 0.2, groupBy = "library_id", ...) {
    set.seed(123)
    proj <- addIterativeLSI(
	ArchRProj = proj,
	useMatrix = "TileMatrix",
	name = "IterativeLSI",
	force = T
    )

    proj = addUMAP(
	proj,
	reducedDims = "IterativeLSI",
	name = "UMAP",
	# nNeighbors = 30,
	# minDist = 0.5,
	metric = "cosine",
	force = T,
	seed = seed
    )

    proj = addHarmony(
	ArchRProj = proj,
	reducedDims = "IterativeLSI",
	name = "Harmony",
	groupBy = groupBy,
	force = T,
	seed = seed,
	corCutOff = harmony_corcutoff,
	...
    )

    proj = addClusters(
	proj,
	reducedDims = "Harmony",
	method = "Seurat",
	name = "Clusters",
	resolusion = cluster_resolusion,
	force = T,
	seed = seed
    )

    proj = addUMAP(
	proj,
	reducedDims = "Harmony",
	name = "UMAP_Harmony",
	# nNeighbors = 30,
	# minDist = minDist,
	metric = "cosine",
	force = T,
	seed = seed
    )

    proj
}
