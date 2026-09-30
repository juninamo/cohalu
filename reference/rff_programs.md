# Find spatial gene programs (recommended entry point)

## Usage

``` r
rff_programs(
  binned,
  offset = NULL,
  fit = NULL,
  n_factors = 6,
  lengthscale = 8,
  ard = 100,
  fit_args = list(),
  null = c("shift", "parametric", "none"),
  highpass = NULL,
  bandwidth = NULL,
  n_boot = 99,
  alpha = 0.05,
  factor_test = FALSE,
  factor_n_boot = 19,
  min_r = 0.5,
  k = 6,
  offset_bandwidth = 10,
  seed = 1,
  n_cores = 1,
  report = NULL,
  crossfit = FALSE,
  control = NULL,
  min_effect = NULL,
  ls_grid = 5 * 2^(0:6),
  profile = list(),
  ...
)
```

## Arguments

- binned:

  A \`binned_transcripts\` object.

- offset:

  Known structure: \`NULL\` / \`"area"\` (none), \`"kmeans"\` or
  \`"pca"\` (\[rff_offset()\] with \`method\`, \`k\`),
  \`"smoothed_total"\` (the log of the Gaussian-smoothed total
  transcript density, bandwidth \`offset_bandwidth\`, as covariate), a
  data frame of covariates (one row per bin or in-tissue bin;
  \[rff_offset()\] with \`method = "covariates"\`), or a matrix from
  \[rff_offset()\] (bins x genes).

- fit:

  Optional existing \[fit_spatial_rff()\] fit to reuse (then \`offset\`,
  \`n_factors\`, \`lengthscale\`, \`ard\` and \`fit_args\` are not used
  for fitting; \`offset\` only documents how the fit's offset was made).

- n_factors, lengthscale, ard:

  Fit settings (see \[fit_spatial_rff()\]); \`lengthscale = "profile"\`
  estimates one length scale per program by held-out likelihood
  (\[rff_lengthscale_profile()\]).

- fit_args:

  Further arguments to \[fit_spatial_rff()\] (they override the defaults
  above, e.g. \`list(max_iter = 300)\`).

- null:

  Null of \[rff_program_test()\]: \`"shift"\` (recommended on real
  data), \`"parametric"\`, or \`"none"\` (no program test).

- highpass, bandwidth:

  Passed to \[rff_program_test()\].

- n_boot:

  Null data sets of the program test.

- alpha:

  Significance level.

- factor_test:

  Also run \[rff_factor_test()\] (sequential), with \`factor_n_boot\`
  refits.

- factor_n_boot:

  Refits for the factor test.

- min_r:

  Minimum \|r\| between a detection axis score and a program map field
  for the two to be matched.

- k:

  Number of clusters / components for \`offset = "kmeans"\` / \`"pca"\`.

- offset_bandwidth:

  Bandwidth for \`offset = "smoothed_total"\`.

- seed:

  Random seed.

- n_cores:

  Cores for the factor test refits.

- report:

  Optional path: write \[rff_report()\] of the result there.

- crossfit, control, min_effect:

  Calibrated program test (see \[rff_program_test()\]): cross-fitted
  held-out effect sizes, comparison with negative-control windows
  (\[rff_control_reference()\]) and an effect-size threshold. With any
  of them, a detection axis counts as significant only when its \`call\`
  is \`TRUE\`. Recommended on real tissue, where the plain gene-shift
  test is rejected by almost any shared residual structure.

- ls_grid, profile:

  With \`lengthscale = "profile"\`: length scales to profile and
  settings for \[rff_lengthscale_profile()\] (see
  \[fit_spatial_rff()\]); the programs table then reports each program
  map's profiled length scale with its 95

  ...Passed to \[rff_report()\] when \`report\` is given.

An object of class \`rff_programs\`: a list with \`programs\` (data
frame, one row per program; see Details), \`fit\`, \`program_test\`,
\`factor_test\`, \`offset_info\` (how the offset was made),
\`settings\`, \`call\` and \`elapsed\` (seconds). The programs table has
columns \`program\` (P1, P2, ... for Confirmed / Candidate / Exploratory
programs), \`status\`, \`detected_by\`, \`p\`, \`p_factor_test\`,
\`detection_axis\` (PC of the program test), \`program_map\` (factor of
the fit), \`r\` (\|r\| between the two), \`best_map\` / \`best_r\` (best
map of an axis even if not matched), \`lengthscale\`,
\`program_strength\`, \`uniform_share\`, \`variance_share\`, \`n_genes\`
(genes with \|loading\| \>= 25 of the largest), \`top_up\`,
\`top_down\`, \`key\` (unique row name) and \`map_label\` /
\`axis_label\` (M\* / D\* names of the report), \`lengthscale_lower\` /
\`lengthscale_upper\` / \`heldout_dev_explained\` (profiled length
scales), \`excess\` / \`excess_lower\` / \`excess_upper\` /
\`p_control\` (calibrated test; \`NA\` otherwise). Attributes
\`loadings\` (genes x rows, the loadings used for interpretation) and
\`fields\` (all bins x rows: the program map field, or the axis score
for Candidates; \`NA\` outside the tissue). Per-cell values: use
\`rff_fields(res\$fit, cells)\` and the \`program_map\` column
(Confirmed programs), see Examples. \*\*Experimental.\*\* Runs the
recommended residual-RFLVM pipeline in one call and returns one table of
\*\*programs\*\*: 1. known structure as offset (\[rff_offset()\], from a
shortcut, covariates or a ready matrix); 2. the random-feature factor
model (\[fit_spatial_rff()\]; grid basis, fixed length scale,
residual-PCA initialisation, ARD); 3. the program test
(\[rff_program_test()\], gene-shift null by default, sequential),
optionally also \[rff_factor_test()\]; 4. matching of the two kinds of
components into one Programs table. Which function should I use? Most
users call \`rff_programs()\` and then \[rff_report()\] (or \`report =
"file.html"\`). \[fit_spatial_rff()\], \[rff_program_test()\],
\[rff_factor_test()\] and \[rff_offset()\] are the building blocks, for
custom pipelines.

When to use / limitations Use it to \*\*discover\*\* multi-gene spatial
programs beyond known structure without naming genes or structures in
advance, then validate them by \[rff_transfer_test()\] on held-out
patients or cohorts. On real tissue, always pass \`control\`
(negative-control regions, \[rff_control_reference()\]) with \`crossfit
= TRUE\`: the gene-shift test alone called almost every real window. The
programs found are only as new as the offset is rich: with lineage
labels only, discovered programs were cell-type / sub-lineage residuals
or technical axes (use fine labels, \`gene_covariates\` for spill-over
and the nuclear share in \[rff_offset()\]) - inspect each with
\[rff_report()\]. For a known question (distance to a structure, a known
gene set) prefer supervised tools (per-gene GLM on distance, C-SIDE,
\[pcf_cross()\]). Report length scales as ranks only (see
\[rff_lengthscale_profile()\]). See the article "When to use the
residual RFLVM (and when not)".

Detection axes and program maps The program test finds \*\*detection
axes\*\* (principal components of smoothed residuals, \`PC1\`, ...) and
decides whether and how many programs exist; the fit gives \*\*program
maps\*\* (factors, \`factor1\`, ...): a smooth field and gene loadings
used for interpretation. Each significant detection axis is matched to
the program map whose field correlates best with its score (greedy by
\|r\|, every map used once): \* \`Confirmed\`: significant axis matched
with \|r\| \>= \`min_r\`; p from the axis, genes, length scale and map
from the program map. \* \`Candidate\`: significant axis without such a
map - the RFLVM does not represent it well; genes and map come from the
axis itself. \* \`Not supported\`: a strong program map (large program
strength, or significant in the factor test) without a significant
matching axis. \* \`Cellularity/technical\`: program maps whose loadings
are mostly shared by all genes (uniform share \> 0.5). With only the
factor test (\`null = "none"\`, \`factor_test = TRUE\`), programs are
the significant program maps (\`detected_by = "factor test"\`); without
any test they are \`Exploratory\`.

\[rff_report()\], \[fit_spatial_rff()\], \[rff_program_test()\],
\[rff_factor_test()\], \[rff_offset()\]
