# run_BM_archR_analysis

ArchR-based analysis scripts for the bone marrow coCnT histone-mark dataset. This repository is script-centric: each `run_*.R` file performs one analysis step and writes intermediate objects to `./tmp/` and figures to `./figures/`.

## Repository layout

- `run_*.R`: analysis scripts, roughly ordered by stage and script number
- `tmp/`: generated ArchR projects, matrices, metadata tables, and other intermediate outputs
- `figures/`: exported plots and figure panels

## External dependencies

This repository depends on data and helper code outside this directory.

- Input sample sheets are read from `../data/`
- Arrow files are read from `../run_BM_coCnT_create_arrow/`
- Some scripts read outputs from sibling projects such as `../run_BM_souporcell/` and `../app_trajectory_shiny/`
- Utility functions are sourced from `../lib/R/`

## R packages

Commonly used packages include:

- `ArchR`
- `data.table`
- `ggplot2`
- `ggpubr`
- `magrittr`
- `readr`
- `stringr`
- `ggsci`
- `httpgd`
- `Cairo`

## Numbering convention

The script prefixes map to broad analysis stages:

- `run_101` to `run_199`: project creation, QC, clustering, annotation, imputation, and early summary analyses
- `run_201` to `run_244`: HSPC-focused re-analysis, trajectory inference, chromatin-state assignment, and downstream comparisons
- `run_299*`: extended chromatin-state and transition analyses
- `run_999*`: exploratory, ad hoc, or supplementary analyses

The scripts are not fully standalone. In practice, later scripts expect files produced by earlier scripts to already exist in `./tmp/`.

## Main workflow

### 1. Build ArchR projects

Use these scripts to assemble ArchR projects from Arrow files and metadata.

- `run_101_create_proj.R`: create single-mark ArchR projects and subset by histone mark
- `run_101_create_proj_cooc.R`: create co-occurrence ArchR projects and subset by antibody combination
- `run_101_gene_cpg_island_promoter.R`: generate promoter CpG-island gene annotations used downstream

### 2. QC and iterative clustering

These scripts perform the main iterative cleanup and annotation workflow for the whole dataset.

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
- `run_119_qc_final.R`

Typical outputs from this stage are cleaned ArchR projects such as `./tmp/H3K27me3_coCnT_final3/` and matching metadata and QC figures.

### 3. Gene-score imputation and matched-mark analyses

These scripts align cells across marks and generate imputed gene-score matrices.

- `run_131_imputation_gene_score.R`
- `run_132_imputation_gene_score_cooc.R`
- `run_133_draw_imputed_gene_score.R`
- `run_134_draw_imputed_gene_score_cooc.R`
- `run_141_update_metadata.R`
- `run_142_gene_score_pseudo_bulk.R`
- `run_151_nfrag_cooc.R`
- `run_152_nfrag.R`
- `run_153_qc_final_cooc.R`
- `run_199_gene_score_correlation_bulk.R`
- `run_199_gene_score_rna_correlation.R`
- `run_199_rna_pseudo_bulk_match_cocnt_cluster.R`

### 4. HSPC re-analysis and trajectory setup

These scripts subset the main projects to HSPC-related lineages, recluster them, and add lineage trajectories.

- `run_201_re-cluster_HSPC.R`
- `run_202_rna_map_HSPC.R`
- `run_203_qc_final_HSPC.R`
- `run_204_imputation_gene_score_HSPC.R`
- `run_205_draw_imputed_gene_score_HSPC.R`
- `run_206_update_whole_umap_with_HSPC.R`
- `run_207_output_new_celltype_annotation_HSPC.R`
- `run_211_update_metadata.R`
- `run_211_update_metadata_add_trajectory.R`
- `run_212_trajectory_umap.R`
- `run_213_three_lineage_umap.R`

Key outputs include trajectory annotations in the ArchR project metadata and exported tables such as `table_metadata_h3k27me3_final3_with_trajectories.tsv`.

### 5. Chromatin-state analysis

These scripts assign chromatin states from imputed gene scores and compare state dynamics across markers and lineages.

- `run_221_gaussian_mixture_model_cutoff.R`
- `run_222_chromatin_status_assignment.R`
- `run_223_chromatin_status_dynamics_stack_barplot.R`
- `run_223_chromatin_status_transition_probability.R`
- `run_224_chromatin_status_between_marker.R`
- `run_224_chromatin_status_between_marker_venn_plot.R`
- `run_225_chromatin_status_consensus_trajectory.R`
- `run_225_chromatin_status_trajectory_heatmap.R`
- `run_226_chromatin_status_state.R`
- `run_227_x2active.R`

Representative outputs:

- `./tmp/table_gene_chromatin_status_by_celltype_me1.tsv`
- `./tmp/table_gene_chromatin_status_by_celltype_me2.tsv`
- `./tmp/table_gene_chromatin_status_by_celltype_me3.tsv`
- `./tmp/table_gene_chromatin_status_wide_me1.tsv`
- `./tmp/table_gene_chromatin_status_wide_me2.tsv`
- `./tmp/table_gene_chromatin_status_wide_me3.tsv`

### 6. Pseudotime and downstream comparisons

These scripts summarize gene-score behavior along trajectories and compare chromatin patterns with RNA.

- `run_231_update_metadata.R`
- `run_232_gene_score_correlation.R`
- `run_233_pseudotime_single_dim.R`
- `run_234_pseudotime_single_gene_v3_plus.R`
- `run_235_pseudotime_single_gene_v3.R`
- `run_236_sc_reads_map.R`
- `run_237_active2repressive_rna.R`
- `run_238_active2repressive_cpg_percentage.R`
- `run_239_promoter_reads_percent.R`
- `run_240_gene_score_scRNA_gene.R`
- `run_241_gene_score_pseudo_bulk_RNA_genesets.R`
- `run_242_timing_of_trajectory.R`
- `run_243_chromatin_vs_RNA.R`
- `run_244_chromatin_state_enrichment.R`

### 7. Extended and exploratory scripts

The `run_299*` and `run_999*` scripts extend the main analysis with supplementary views, overlap analyses, gene-module summaries, and one-off figure generation. They are useful references for downstream analyses but should be treated as dependent on prior outputs rather than as a clean linear pipeline.

## Running scripts

Scripts are generally run directly with `Rscript` from this directory:

```bash
Rscript run_101_create_proj.R
Rscript run_102_cluster_cell_round0.R
Rscript run_222_chromatin_status_assignment.R
```

Recommended practice:

- run scripts from the repository root so relative paths resolve correctly
- inspect the top of each script for required input files and expected output directories
- do not assume every script is idempotent; many overwrite `./tmp/` outputs

## Notes

- The README intentionally documents the current script inventory rather than a manuscript-style figure map, because the codebase has evolved beyond the original figure-specific list.
- Some script headers still refer to older output names such as `final2` and `final3`; those names reflect the iterative project versions used during analysis.
