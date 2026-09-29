# Package index

## Simulation

Generate synthetic tissues with a planted cell-cell interaction, and
convert them to Seurat objects.

- [`generate_sim()`](https://juninamo.github.io/cohalu/reference/generate_sim.md)
  : Simulate Spatial Coordinates and Cell Types
- [`generate_sim_groups()`](https://juninamo.github.io/cohalu/reference/generate_sim_groups.md)
  : Generate multi-sample simulated data with disease-group structure
- [`sim_to_seurat()`](https://juninamo.github.io/cohalu/reference/sim_to_seurat.md)
  : Convert simulated cells to a Seurat object with one FOV per sample

## Single-sample analysis

Spatial neighborhood analysis (SNA), spatial co-localization score
(sCLS), radius-based co-occurrence and interaction spots for one tissue.

- [`nhood_enrichment()`](https://juninamo.github.io/cohalu/reference/nhood_enrichment.md)
  : Neighborhood Enrichment (Generic method)
- [`plot_nhood_heatmap()`](https://juninamo.github.io/cohalu/reference/plot_nhood_heatmap.md)
  : Heatmap of a single-sample neighbourhood enrichment with
  significance
- [`nhood_enrichment.Seurat()`](https://juninamo.github.io/cohalu/reference/nhood_enrichment.Seurat.md)
  : Neighborhood Enrichment (Seurat Method)
- [`cooccur_local()`](https://juninamo.github.io/cohalu/reference/cooccur_local.md)
  : Local Co-occurrence Score (Generic method)
- [`cooccur_local.Seurat()`](https://juninamo.github.io/cohalu/reference/cooccur_local.Seurat.md)
  : Local Co-occurrence Score (Seurat Method)
- [`cooccur_local_oe()`](https://juninamo.github.io/cohalu/reference/cooccur_local_oe.md)
  : Local co-localization score with a permutation-calibrated
  expectation
- [`calc_co_occurrence_for_radius()`](https://juninamo.github.io/cohalu/reference/calc_co_occurrence_for_radius.md)
  : Calculate Co-occurrence Matrix for a Given Radius
- [`compute_co_occurrence_ratio()`](https://juninamo.github.io/cohalu/reference/compute_co_occurrence_ratio.md)
  : Compute Enrichment Ratios from Count Matrix
- [`search_interaction_spot()`](https://juninamo.github.io/cohalu/reference/search_interaction_spot.md)
  : Search for Spatial Interaction Spots

## Multi-sample scores

Compute a score per image (or patient) and return a tidy table ready for
group comparison.

- [`build_sample_design()`](https://juninamo.github.io/cohalu/reference/build_sample_design.md)
  : Build a sample design table for disease-group comparisons
- [`nhood_enrichment_per_sample()`](https://juninamo.github.io/cohalu/reference/nhood_enrichment_per_sample.md)
  : Per-sample neighborhood enrichment
- [`cooccur_ratio_per_sample()`](https://juninamo.github.io/cohalu/reference/cooccur_ratio_per_sample.md)
  : Per-sample radius-based co-occurrence ratio
- [`cooccur_local_per_sample()`](https://juninamo.github.io/cohalu/reference/cooccur_local_per_sample.md)
  : Per-sample local co-occurrence summary
- [`interaction_spot_per_sample()`](https://juninamo.github.io/cohalu/reference/interaction_spot_per_sample.md)
  : Per-sample interaction-spot summary
- [`summarize_by_patient()`](https://juninamo.github.io/cohalu/reference/summarize_by_patient.md)
  : Aggregate per-image scores to one row per patient

## Segmentation-free analysis (experimental)

Work directly on transcript coordinates: bin transcripts, estimate the
cross pair correlation of gene sets model-free, or fit a random-feature
log-Gaussian Cox process factor model (after Gundersen et al. 2021).

- [`read_xenium_transcripts()`](https://juninamo.github.io/cohalu/reference/read_xenium_transcripts.md)
  : Read a Xenium transcript table
- [`bin_transcripts()`](https://juninamo.github.io/cohalu/reference/bin_transcripts.md)
  : Bin transcript coordinates on a square grid
- [`pcf_cross()`](https://juninamo.github.io/cohalu/reference/pcf_cross.md)
  : Empirical cross pair correlation of two gene sets
  (segmentation-free)
- [`pcf_matrix()`](https://juninamo.github.io/cohalu/reference/pcf_matrix.md)
  : Cross pair correlation for every pair of gene sets
- [`colocalization_gene_matrix()`](https://juninamo.github.io/cohalu/reference/colocalization_gene_matrix.md)
  : Gene-by-gene co-localization of transcripts
- [`colocalization_modules()`](https://juninamo.github.io/cohalu/reference/colocalization_modules.md)
  : Cluster genes into co-localization modules
- [`module_enrichment()`](https://juninamo.github.io/cohalu/reference/module_enrichment.md)
  : Interpret co-localization modules with gene sets
- [`module_enrichr()`](https://juninamo.github.io/cohalu/reference/module_enrichr.md)
  : Query enrichR for each co-localization module
- [`rff_programs()`](https://juninamo.github.io/cohalu/reference/rff_programs.md)
  : Find spatial gene programs (recommended entry point)
- [`fit_spatial_rff()`](https://juninamo.github.io/cohalu/reference/fit_spatial_rff.md)
  : Fit a random-feature spatial factor model to binned transcripts
- [`rff_fields()`](https://juninamo.github.io/cohalu/reference/rff_fields.md)
  : Latent fields of a fitted random-feature LGCP at cells, transcripts
  or any point
- [`rff_offset()`](https://juninamo.github.io/cohalu/reference/rff_offset.md)
  : Offset for a residual random-feature model: expression explained by
  known structure
- [`rff_factor_test()`](https://juninamo.github.io/cohalu/reference/rff_factor_test.md)
  : Significance of residual spatial factors by parametric bootstrap
- [`rff_program_test()`](https://juninamo.github.io/cohalu/reference/rff_program_test.md)
  : Gene-program test on smoothed residuals, with the fitted RFLVM as
  null model
- [`fit_spatial_rff_joint()`](https://juninamo.github.io/cohalu/reference/fit_spatial_rff_joint.md)
  : Random-feature factor model across several windows with shared
  loadings
- [`rff_program_test_joint()`](https://juninamo.github.io/cohalu/reference/rff_program_test_joint.md)
  : Gene-program test across several tissues or windows with shared
  loadings
- [`rff_expected_offset()`](https://juninamo.github.io/cohalu/reference/rff_expected_offset.md)
  : Offset from the cell types that own each bin's transcripts
- [`rff_report()`](https://juninamo.github.io/cohalu/reference/rff_report.md)
  : Interactive HTML report of a residual RFLVM analysis
- [`rff_pair_correlation()`](https://juninamo.github.io/cohalu/reference/rff_pair_correlation.md)
  : Model-based cross pair correlation from a fitted random-feature LGCP
- [`colocalization_per_sample()`](https://juninamo.github.io/cohalu/reference/colocalization_per_sample.md)
  : Per-sample segmentation-free co-localization of gene sets
- [`simulate_transcripts()`](https://juninamo.github.io/cohalu/reference/simulate_transcripts.md)
  : Simulate transcript coordinates from a multivariate log-Gaussian Cox
  process
- [`simulate_transcripts_groups()`](https://juninamo.github.io/cohalu/reference/simulate_transcripts_groups.md)
  : Simulate a multi-sample case-control transcript dataset
- [`lgcp_true_pair_correlation()`](https://juninamo.github.io/cohalu/reference/lgcp_true_pair_correlation.md)
  : True cross pair correlation of a simulated LGCP

## Group comparison

Test every cell-type pair between two groups.

- [`compare_groups()`](https://juninamo.github.io/cohalu/reference/compare_groups.md)
  : Compare disease groups across samples
- [`associate_continuous()`](https://juninamo.github.io/cohalu/reference/associate_continuous.md)
  : Associate co-localization with a continuous clinical variable

## Visualization

- [`plot_group_delta_heatmap()`](https://juninamo.github.io/cohalu/reference/plot_group_delta_heatmap.md)
  : Heatmap of the disease-group effect for every cluster pair
- [`plot_pair_boxplot()`](https://juninamo.github.io/cohalu/reference/plot_pair_boxplot.md)
  : Per-sample boxplot of one (or several) cluster pair(s) across groups
- [`plot_volcano_groups()`](https://juninamo.github.io/cohalu/reference/plot_volcano_groups.md)
  : Volcano plot of a compare_groups result
- [`manual_colors`](https://juninamo.github.io/cohalu/reference/manual_colors.md)
  : Default Manual Colors for Clusters

## Low-level helpers

- [`compute_count()`](https://juninamo.github.io/cohalu/reference/compute_count.md)
  : Compute Co-occurrence Count Matrix
- [`permute_clusters()`](https://juninamo.github.io/cohalu/reference/permute_clusters.md)
  : Permute Cluster Assignments and Recompute Counts
