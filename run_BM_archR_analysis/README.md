# run_BM_archR_analysis

ArchR-based analysis scripts for the bone marrow coCnT histone-mark dataset. 
Each `run_*.R` file performs one analysis step and writes intermediate objects to `./tmp/` and figures to `./figures/`.

## Repository layout

- `run_*.R`: analysis scripts in numeric order
- `tmp/`: generated ArchR projects, matrices, metadata tables, and other intermediate outputs
- `figures/`: exported plots and figure panels

## Dependencies

- `../data/`: sample sheets, Seurat objects
    You need to download the data folder from Zenodo, please refer to the main README for details.
- `../run_BM_coCnT_create_arrow/`: Arrow files used when building ArchR projects
- `../run_BM_souporcell/`: donor-assignment and barcode-refinement inputs
- `../run_BM_bed_profile/`: fragment / bed-profile summaries used by downstream plotting scripts
- `../lib/R/`: helper R functions sourced by some scripts

## R packages

Common packages used across the current scripts include:

- `ArchR`
- `data.table`
- `ggplot2`
- `ggpubr`
- `magrittr`
- `readr`
- `stringr`
- `ggsci`
- `Seurat`
- `GenomicRanges`
- `GenomicFeatures`
- `mclust`
- `clusterProfiler`
- `org.Hs.eg.db`

## Current script ranges

The numbering is historical rather than continuous. In the current checkout, the available scripts are:

- `run_101` to `run_199`: project creation, QC, clustering, annotation, imputation, and pseudo-bulk analysis
- `run_201` to `run_245`: HSPC-focused re-analysis, chromatin-state analyses, and pseudotime / RNA comparison

The scripts are not fully standalone. In practice, later scripts expect files produced by earlier scripts to already exist in `./tmp/`.

## Main workflow

### 1. Project creation

- `run_101_create_proj.R`: build the primary single-mark ArchR projects
- `run_101_gene_cpg_island_promoter.R`: derive the promoter/CpG-island gene list used in chromatin-state annotation
- `run_121_create_proj_cooc.R`: build co-occurrence ArchR projects

### 2. Whole-BM QC, clustering, and annotation

- `run_102_cluster_cell_round0.R`
- `run_103_label_doublet.R`
- `run_104_remove_doublet.R`
- `run_105_cluster_cell_round1.R`
- `run_106_identify_bad_cluster_round1.R`
- `run_107_remove_bad_cluster1.R`
- `run_108_cluster_cell_round2.R`
- `run_109_join_umap1.R`
- `run_110_get_signature_genes_round1.R`
- `run_111_rna_map1.R`
- `run_112_remove_bad_cluster2.R`
- `run_113_cluster_cell_round3.R`
- `run_114_join_umap2.R`
- `run_115_get_signature_genes_round2.R`
- `run_116_rna_map2.R`
- `run_117_join_umap3.R`
- `run_118_annotation_qc.R`
- `run_135_qc_final.R`

### 3. Gene-score imputation and pseudo-bulk summaries

- `run_131_imputation_gene_score.R`
- `run_132_imputation_gene_score_cooc.R`
- `run_133_draw_imputed_gene_score.R`
- `run_141_update_metadata.R`
- `run_142_gene_score_pseudo_bulk.R`
- `run_199_rna_pseudo_bulk_match_cocnt_cluster.R`

### 4. HSPC annotation refinement and trajectory setup

- `run_201_re-cluster_HSPC.R`
- `run_204_imputation_gene_score_HSPC.R`
- `run_205_draw_imputed_gene_score_HSPC.R`
- `run_206_update_whole_umap_with_HSPC.R`
- `run_207_output_new_celltype_annotation_HSPC.R`
- `run_211_update_metadata_add_trajectory.R`
- `run_212_trajectory_umap.R`
- `run_213_three_lineage_umap.R`

### 5. Chromatin-state assignment and QC

- `run_221_gaussian_mixture_model_cutoff.R`
- `run_222_chromatin_status_assignment.R`
- `run_223_chromatin_status_dynamics_stack_barplot.R`
- `run_224_chromatin_status_between_marker_venn_plot.R`
- `run_227_x2active.R`
- `run_228_chromatin_status_qc.R`

### 6. Pseudotime and RNA comparison

- `run_233_pseudotime_single_dim.R`
- `run_236_sc_reads_map.R`
- `run_240_gene_score_scRNA_gene.R`
- `run_241_gene_score_pseudo_bulk_RNA_genesets.R`
- `run_244_chromatin_state_enrichment.R`
- `run_245_pseudotime_heatmap.R`

These scripts compare chromatin-state behavior, gene scores, and scRNA-seq expression along the lineage trajectories defined in the HSPC re-analysis.
