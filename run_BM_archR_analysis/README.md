# run_BM_archR_analysis

ArchR-based analysis scripts for the bone marrow coCnT histone-mark dataset. This checkout is script-centric: each `run_*.R` file performs one analysis step and writes intermediate objects to `./tmp/` and figures to `./figures/`.

## Repository layout

- `run_*.R`: analysis scripts in historical numeric order
- `tmp/`: generated ArchR projects, matrices, metadata tables, and other intermediate outputs
- `figures/`: exported plots and figure panels
- `bk/`: local backup or scratch files kept alongside the main analysis
- `ArchRLogs/`: ArchR log output

## External dependencies

This directory depends on data and helper code stored in sibling locations.

- `../data/`: sample sheets, Seurat objects, gene lists, and reference tables
- `../run_BM_coCnT_create_arrow/`: Arrow files used when building ArchR projects
- `../run_BM_souporcell/`: donor-assignment and barcode-refinement inputs
- `../run_BM_bed_profile/`: fragment / bed-profile summaries used by downstream plotting scripts
- `../app_trajectory_shiny/`: metadata tables reused when exporting trajectory annotations
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

- `run_101` to `run_142`: project creation, QC, clustering, annotation, imputation, and pseudo-bulk summaries
- `run_199` to `run_245`: RNA pseudo-bulk mapping, HSPC-focused re-analysis, chromatin-state analyses, and pseudotime / RNA comparison

No `run_299*` or `run_999*` scripts are present in this checkout.

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

This block covers doublet labelling/removal, iterative clustering, UMAP joins, marker-gene identification, RNA mapping, and final QC summaries across the cleaned single-mark and co-occurrence projects.

### 3. Gene-score imputation and pseudo-bulk summaries

- `run_131_imputation_gene_score.R`
- `run_132_imputation_gene_score_cooc.R`
- `run_133_draw_imputed_gene_score.R`
- `run_141_update_metadata.R`
- `run_142_gene_score_pseudo_bulk.R`
- `run_199_rna_pseudo_bulk_match_cocnt_cluster.R`

Representative outputs from this stage include imputed gene-score matrices, refreshed UMAP / metadata tables, `./tmp/table_gene_marker_average_score_by_celltype.tsv.gz`, and RNA pseudo-bulk tables derived from the Triana 2021 dataset.

### 4. HSPC re-analysis and trajectory setup

- `run_201_re-cluster_HSPC.R`
- `run_204_imputation_gene_score_HSPC.R`
- `run_205_draw_imputed_gene_score_HSPC.R`
- `run_206_update_whole_umap_with_HSPC.R`
- `run_207_output_new_celltype_annotation_HSPC.R`
- `run_211_update_metadata_add_trajectory.R`
- `run_212_trajectory_umap.R`
- `run_213_three_lineage_umap.R`

These scripts subset HSPC-like populations, recluster them, propagate HSPC annotations back into the whole-BM view, and add lineage trajectory information to the final H3K27me3 project.

### 5. Chromatin-state assignment and QC

- `run_221_gaussian_mixture_model_cutoff.R`
- `run_222_chromatin_status_assignment.R`
- `run_223_chromatin_status_dynamics_stack_barplot.R`
- `run_224_chromatin_status_between_marker_venn_plot.R`
- `run_227_x2active.R`
- `run_228_chromatin_status_qc.R`

Key outputs include the `table_gene_chromatin_status_*` matrices in `./tmp/` plus QC and overlap figures summarizing state composition across cell types.

### 6. Pseudotime and RNA comparison

- `run_233_pseudotime_single_dim.R`
- `run_236_sc_reads_map.R`
- `run_240_gene_score_scRNA_gene.R`
- `run_241_gene_score_pseudo_bulk_RNA_genesets.R`
- `run_244_chromatin_state_enrichment.R`
- `run_245_pseudotime_heatmap.R`

These scripts compare chromatin-state behavior, gene scores, and scRNA-seq expression along the lineage trajectories defined in the HSPC re-analysis.

## Running scripts

Scripts are generally run directly with `Rscript` from this directory:

```bash
Rscript run_101_create_proj.R
Rscript run_135_qc_final.R
Rscript run_222_chromatin_status_assignment.R
Rscript run_245_pseudotime_heatmap.R
```

Recommended practice:

- run scripts from the repository root so relative paths resolve correctly
- inspect the top of each script for required input files and expected output locations
- do not assume every script is idempotent; many overwrite `./tmp/` outputs
- treat gaps in the numbering as historical artifacts, not missing mandatory pipeline stages

## Notes

- This README documents the scripts that are present in the current checkout.
- Earlier revisions of the project referenced additional script ranges that are not in this directory anymore.
- Some paths and object names still use iterative labels such as `final2` and `final3`; those reflect analysis snapshots rather than separate repositories.
