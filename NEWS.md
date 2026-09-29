# cohalu 0.99.3

## Package renamed

* The package is renamed from `spatialCooccur` to `cohalu` (COHALU:
  CO-localization, Hotspots And sample-Level Units; read "koharu"). Use
  `library(cohalu)`; function names and arguments are unchanged. The internal
  class `spatialCooccurSample` is now `cohaluSample`.

## Changes that affect results

* `nhood_enrichment()`, `nhood_enrichment.Seurat()` and
  `nhood_enrichment_per_sample()`: the default is now `transformation = FALSE`
  (every kNN link counts once, as in squidpy). The previous default weighted
  each link from cell u by 1 / (1 + d_u), d_u = number of cells that chose u.
  With a kNN graph every cell sends k links at any density, so this was not a
  density correction; in simulations it kept the calibration (false positives
  4-5% per pair, 2-7% family-wise) but lowered power (planted pair, 1,500
  cells: 14% with the weighting vs 28% without). Use `transformation = TRUE`
  to reproduce earlier results. The docs now describe what the weighting does.
* `nhood_enrichment()`: `log2_oe` is now centred on the label shuffles
  (the mean of the same log ratio over the shuffles is subtracted), so it is
  0 on average without interaction for any number of cells. The log of a
  ratio of small counts was biased below 0 for rare cell types (about -0.03
  for 25 cell types of 60 cells, -0.05 for pairs with < 50 cells). The
  previous value is returned as `log2_oe_raw`.
* `generate_sim(test_type = "distribute")` no longer places relocated cells
  outside the tissue. Such cells sat in empty space where their k nearest
  neighbours reached their partner cell even at 100 um, so co-localization
  leaked to long planted distances.
* `compare_groups()`: new `unit = c("patient", "image")`, default
  `"patient"`. With `patient_key` and `method = "wilcox"` or `"t"`, images
  are averaged within patient before testing, so the patient is the unit of
  analysis. Previously these tests used image-level rows (pseudoreplication,
  13.8% false positives at 5% in simulations) unless the data were first
  passed through `summarize_by_patient()`; use `unit = "image"` for the old
  behaviour.
* `compare_groups()` now defaults to `value = "log2_oe"` and uses a
  `patient` column automatically when `patient_key` is not given and
  patients have several images.
* `compare_groups(symmetric = TRUE)` averages the (i, j) and (j, i) values
  of each sample instead of keeping only the `cluster_i <= cluster_j` row
  (degree-normalised neighbourhood scores are slightly directional).
* `summarize_by_patient()` keeps distances (`r`) of
  `colocalization_per_sample()` output separate instead of averaging them;
  new `pair_keys` argument.

## New features

* Several windows at once (experimental, in development):
  `rff_program_test_joint()` tests gene programs across several tissues or
  windows with shared loadings, `fit_spatial_rff_joint()` fits the factor
  model across windows with shared loadings, and `rff_expected_offset()`
  builds an offset from the cell types that own each bin's transcripts.
  `rff_offset(base = )` adds covariates on top of such an offset, and
  `fit_spatial_rff(loadings = )` maps given programs in a new tissue (fixed
  loadings; per-tissue `amplitude`). In simulations (8 windows, 30 genes,
  4-gene program, gene-own residual fields, `highpass = 20`), the joint test
  found a program of amplitude 0.35 / 0.5 in 13/20 / 20/20 window sets while
  per-window tests found it in 0% / 10% of windows; without a program it
  called 0/20 (with or without gene-own fields); held-out confirmation
  (`loadings =`) was significant in 18/20 / 20/20. On 195 real TLS windows,
  gene-shifted sanity data were called in 30/585 replicates (5.1%).
  `rff_expected_offset(neighbour_bandwidth = )` also returns the counts
  expected from the cell types of neighbouring bins (attribute
  `neighbour`), and `rff_offset(gene_covariates = )` adds gene-specific
  covariates (one coefficient per gene) - together they absorb segmentation
  spill-over / mixed cells, which otherwise appear as cell-type-marker
  "programs" in residual analyses.
  `rff_report()` accepts a `fit_spatial_rff_joint()` fit (new argument
  `window`) and the discovery output of `rff_program_test_joint()`.
  `rff_program_test_joint()` now reduces each surrogate draw inside the
  worker (memory) and recomputes draws lost in a parallel worker serially.
* Residual random-feature model (experimental): `fit_spatial_rff()` accepts
  a per-bin, per-gene log `offset` matrix describing structure that is
  already known (cell-type composition, domains, or an embedding such as
  PCA, Harmony or SCIGMA; built with the new `rff_offset()`), so that the
  factors capture only the remaining spatially coherent variation. New
  `ard` argument: a group penalty on each factor's loadings that shrinks
  unneeded factors, reported as `factor_strength`. New `rff_factor_test()`:
  parametric-bootstrap p-value per factor against the largest factor fitted
  to null data (family-wise over factors). The test statistic is the program
  strength (loading norm after removing the loading shared by all genes):
  factors that move all genes together are cellularity, not gene programs,
  and the raw loading norm called such factors significant in simulations
  without any program. `fit_spatial_rff()` also returns `program_strength`.
* `fit_spatial_rff()` is several times faster on binned data (about 5x for
  256 random features; same objective to numerical precision): on the bin
  grid, `cos/sin(w_x x + w_y y)` factorise into x- and y-parts, so the
  fields, their length-scale derivatives and the weight gradients are
  small matrix products instead of bins x features trigonometric
  evaluations. Irregular coordinates use the previous direct evaluation.
* `fit_spatial_rff()`: new `factor_init = "residual_pca"` (start the factors
  from a PCA of smoothed Pearson residuals under the offset), `init` (warm
  start of the intercepts, dispersions and cellularity field from an earlier
  fit) and `weights` (per-entry likelihood weights, e.g. for held-out
  cross-validation). The fit settings now include `factor_init`.
  `rff_factor_test()` does not warm-start its refits: in simulations,
  warm-started null refits had much smaller statistics than the cold-started
  observed fit (false positives in 5 of 10 null tissues).
* `rff_factor_test()`: new `n_cores` (parallel refits by forking; the null
  data sets are simulated first, so results do not depend on `n_cores`).
  A `max_iter` different from the fit's now triggers a warning: the
  optimiser usually stops at `max_iter` and the factor statistics grow with
  the number of iterations, so null refits with fewer iterations than the
  observed fit make the test anti-conservative. The run time is returned as
  attribute `elapsed`.
* `fit_spatial_rff()`: new `basis = "grid"` - each field is the RBF Gaussian
  process on the zero-padded bin grid, applied by FFT (circulant embedding),
  instead of `n_features` random frequencies; new `learn_lengthscales`
  (fixed factor length scales) and `l1` (L1 penalty on single loadings, one
  weight per gene allowed). In simulations the grid basis with a fixed 8 um
  length scale recovered 10 um programs better (median |r| 0.74 vs 0.60 at
  amplitude 1; 0.81 with `factor_init = "residual_pca"`, residual PCA 0.84);
  learned length scales shrank to the bin size with the grid basis. `l1`
  did not improve power. `rff_fields()` interpolates grid-basis fields.
* `rff_factor_test()`: new `sequential = TRUE` - the factor ranked `k` is
  compared with refits to data simulated from the null model plus the
  `k - 1` stronger factors (column `p_sequential`). With one planted
  structure, 5/6 simulated tissues called exactly one factor (2-3 with the
  max-null p-values); in a real Xenium synovium window where the max-null
  test called all 6 factors, 3-4 were called.
* New `rff_program_test()` (experimental): residual PCA of smoothed Pearson
  residuals with the fitted RFLVM (offset, intercepts, cellularity field and
  the all-gene part of every factor) as null model and parametric
  bootstrap, tested sequentially. In simulations it matched the power of
  residual PCA (20/20 at amplitude 0.5; `rff_factor_test()` 12-13/20) while
  keeping the calibration (3/60 without program; 1/40 with extra
  cellularity only, residual PCA 12/20).
* `rff_program_test()`: new `null = "shift"` - a null that keeps each gene's
  own spatial autocorrelation. Every gene's counts are moved together with
  its null mean by an independent random flip/transposition and a shift on
  the mirror-extended bin grid; only the alignment between genes is
  destroyed, so significance means a program shared by several genes rather
  than any structure the offset misses. Components are compared with the
  surrogates' spectra (parallel analysis, sequential). With gene-specific
  residual fields (10-40 um) the parametric null called 100% of simulated
  tissues, the gene-shift null 0-10% (0/40 without residual structure);
  power 80-100% for programs of amplitude 0.5-1. New `highpass` (band-pass
  filtering of the residuals) restores power when broad gene-specific fields
  mask small-scale programs. Toroidal shifts without mirror extension were
  miscalibrated (18-32% false calls) and are not offered.
* `rff_offset()`: new `method = c("covariates", "kmeans", "pca")` and `k`
  for data without labels: the offset is built from the smoothed
  composition of `k` k-means clusters of the bins' log-normalised expression
  (a stand-in for cell-type composition) or from its top `k` principal
  components. `covariates` is now optional for these methods (existing calls
  are unchanged). In simulations without labels, `"kmeans"` (`k = 6`) found a
  minor 4-gene program of amplitude 1 in 23/24 tissues (true labels 24/24, no
  offset 0/24); PCA offsets with `k >= 5` absorbed it.
* `rff_fields()` (experimental): evaluates the latent fields of a
  `fit_spatial_rff()` fit (cellularity and the `K` factors, on the
  unit-variance scale of `fit$field_grid`) at any coordinates - cell
  centroids, transcripts or other points - and averages them per cell with
  `by = "cell_id"`. `fit_spatial_rff()` now also stores the coordinate
  centre used for the random features (`center`); older fits still work.
* New `rff_programs()` (experimental; the recommended entry point for
  spatial gene programs): one call that builds the offset (a shortcut
  `"kmeans"` / `"pca"` / `"smoothed_total"` / `"area"`, covariates or a
  matrix; recorded in `offset_info`), fits `fit_spatial_rff()` (grid basis,
  fixed length scale, residual-PCA initialisation, ARD), runs
  `rff_program_test()` (gene-shift null, sequential; optionally
  `rff_factor_test()`) and returns one Programs table: significant detection
  axes (program-test PCs) matched to program maps (fit factors) by |r|,
  with status Confirmed / Candidate / Not supported / Cellularity-technical.
  `print()`, `summary()`, `as.data.frame()`; `rff_report()` accepts the
  result directly.
* New `rff_report()` (experimental): writes a single self-contained,
  interactive HTML report of a (residual) RFLVM analysis - input data and
  offset summary, auto-generated key findings with ok / caution / warning
  flags (significance, cellularity-like factors, calibration caveats of the
  null used, fit stopped at `max_iter`, ...), a sortable factor table, pan /
  zoom maps of every field and of selected genes (observed, null
  expectation, log ratio, factor contribution), loadings heatmap and gene
  search, per-cell summaries by label, overlap with known gene programs, the
  null distributions of `rff_factor_test()` / `rff_program_test()`, methods
  and code to reproduce. With `export = TRUE` (default) the factor table,
  loadings, fields per bin, cell values, test results and null statistics
  (CSV), the fit (RDS), the settings (JSON) and an R script are written to
  `<report>_files/`. The main result is a merged Programs table (as in
  `rff_programs()`); detection axes and program maps are in a collapsed
  "Details & diagnostics" section. No new dependencies.
* `associate_continuous()`: association of per-image / per-patient
  co-localization with a continuous clinical variable (CRP, disease
  activity, age): Spearman on patient means (default), linear model with
  covariates, mixed model on images, or permutation; BH over pairs.
* Unsupervised transcript-level co-localization (experimental):
  `colocalization_gene_matrix()` (gene x gene log2 O/E of transcript pairs
  within a radius, label-shuffling expectation in closed form, one FFT per
  gene), `colocalization_modules()` (clusters co-localizing genes),
  `module_enrichment()` (hypergeometric test of any gene sets, e.g. pathways
  or cell-type markers) and `module_enrichr()` (enrichR wrapper).
* `nhood_enrichment()` returns a within-sample test per unordered pair:
  `pvalue` (normal, from the shuffles), `padj` (Westfall-Young max-T,
  family-wise error rate, calibrated for any number of cell types) and
  `padj_bh`. New `plot_nhood_heatmap()` draws `log2_oe` with significance
  stars.
* `nhood_enrichment()` also returns directional statistics (row = centre
  cell type, column = neighbour type; not symmetric), because the pair-level
  `log2_oe` is symmetric by construction and cannot tell "A is surrounded by
  B" from "B is surrounded by A": `contact` (share of centre cells with at
  least one neighbour of the other type) and `dominance` (share of centre
  cells whose neighbours are at least half of the other type), each with
  `*_expected`, centred `*_log2_oe`, `*_pvalue` and max-T `*_padj` over all
  ordered pairs. `plot_nhood_heatmap(value = "dominance_log2_oe")` draws them
  without symmetrising; the pair-level heatmap is labelled as the average of
  both directions.
* `plot_nhood_heatmap()`: new default `triangle = "auto"` draws symmetric
  pair-level values once (lower triangle with the diagonal) and directional
  values in full; `triangle = "full"` restores the previous layout.
* Tutorials use `log2_oe` with the max-T adjusted `padj`; new sections on
  `cooccur_local_oe()`, `associate_continuous()` and gene-level modules. The
  algorithm reference covers all current methods.

# cohalu 0.99.2

## Bug fixes that change results

* `nhood_enrichment()` / `permute_clusters()`: the permutation null now applies
  a single label permutation to rows and columns. Previously rows and columns
  were shuffled independently, which under spatial randomness inflated
  same-type z-scores (about +11) and biased different-type z-scores (about -2).
  z-scores are not comparable with versions <= 0.99.1.
* `cooccur_local()` / `cooccur_local.Seurat()`: diffusion now iterates from the
  previous step, and the kurtosis early-stop rule compares consecutive steps.
  Previously every step restarted from the raw indicator, so any
  `maxnsteps >= 1` gave the one-step result. Results with `maxnsteps <= 1`
  are unchanged.

## Other bug fixes

* `nhood_enrichment()`, `nhood_enrichment.Seurat()`, `cooccur_local()` and
  `cooccur_local.Seurat()` build the kNN graph exactly with a kd-tree
  (`RANN::nn2`) instead of Seurat's approximate annoy search. On 200,000
  cells the graph takes 0.4 s instead of 33 s, and 99.995% of links are
  identical (on four RA sections, log2 O/E changed by at most 0.016 and the
  significant pairs were identical). `connectivity_key = "snn"` still uses
  `Seurat::FindNeighbors()`.
* `nhood_enrichment()` with `n_jobs > 1`: the workers did not load the Matrix
  methods, so every worker failed and the permutations silently fell back to
  sequential (no speed-up). The workers now load Matrix, and a fallback to
  sequential is reported as a warning.
* `nhood_enrichment()` is now reproducible for a given `seed`, both
  sequentially and with `n_jobs > 1` (worker RNG streams are seeded).
* `compare_groups()` and `nhood_enrichment()` no longer overwrite the caller's
  RNG state.
* Seurat input to the `*_per_sample()` helpers matches images to `meta.data`
  by cell name, so `sample_key` values need not equal the image names.
* `interaction_spot_per_sample()` returns `NA` (not 0 spots) when the spot
  search fails, and no longer requires a `cell` column in `meta.data`.
* Absent cell types give `NA` scores instead of silently dropped rows.

## New features

* `cooccur_local_oe()`: abundance-adjusted local co-localization. For every
  cell, the number of cluster_x-cluster_y pairs within `radius` is divided by
  its exact expectation under label permutation (closed form), smoothed with
  a Gaussian kernel of explicit width, with optional permutation hotspot
  p-values (O(n k) per permutation). `cooccur_local_per_sample()` gains the
  corresponding `log2_oe` summary, recommended for group comparison: the
  mean diffusion sCLS is unchanged by diffusion (mass-conserving) and grows
  with cell-type abundance.

* **Experimental segmentation-free analysis** of transcript coordinates
  (e.g. Xenium `transcripts.parquet`): `bin_transcripts()`, model-free cross
  pair correlation of gene sets `pcf_cross()` (the relative version is the
  label-permutation O/E, the continuous analogue of `log2_oe`), a
  random-feature log-Gaussian Cox process factor model `fit_spatial_rff()`
  after Gundersen, Zhang & Engelhardt (AISTATS 2021) with model-based
  `rff_pair_correlation()`, and `colocalization_per_sample()` whose output
  goes into `compare_groups()` with `pair_keys = c("cluster_i", "cluster_j", "r")`.
  Simulators `simulate_transcripts()` / `simulate_transcripts_groups()` and
  `lgcp_true_pair_correlation()` support validation.
  `read_xenium_transcripts()` reads Xenium transcript tables (binary gene
  names in older outputs, qv filter, gene filtering inside Arrow for 5K
  panels) and `pcf_matrix()` computes all gene-set pairs with cached FFTs.
  Any labelled point set can be analysed the same way, e.g. pixel-level
  factors from FICTURE / punkst (`bin_transcripts(gene_col = "K1")`).
* `compare_groups()` supports paired / repeated-measures designs
  (e.g. pre- vs post-treatment): new `method = "signrank"`, and
  `method = "perm"` now permutes labels within patients when patients appear
  in both groups (the previous between-patient permutation was invalid for
  such designs).

* `nhood_enrichment()` also returns `expected` (permutation mean) and
  `log2_oe` (log2 observed / expected), an effect size that does not grow with
  the number of cells. Recommended for between-group comparison.
* `*_per_sample()` outputs gain `n_cells`, `n_i`, `n_j` (abundance of the two
  cell types).
* `compare_groups()`: exact Wilcoxon p-values for small samples; LMM p-values
  with Satterthwaite df (via lmerTest) instead of Wald z; `covariates`;
  `symmetric`; `min_n_per_group`; exact enumeration for small permutation
  tests; a warning on pseudoreplication.
* `summarize_by_patient()` aggregates image-level results to patients.
* `generate_sim_groups()` supports several images per patient
  (`n_images_per_patient`, `within_patient_noise`) and group-specific
  `n_cells`.
* New tutorial: case-control comparison (`vignettes/case_control_tutorial.ipynb`).
* Faster counting with sparse matrix products.
