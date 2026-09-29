# Find spatial gene programs (recommended entry point)

\*\*Experimental.\*\* Runs the recommended residual-RFLVM pipeline in
one call and returns one table of \*\*programs\*\*: 1. known structure
as offset (\[rff_offset()\], from a shortcut, covariates or a ready
matrix); 2. the random-feature factor model (\[fit_spatial_rff()\]; grid
basis, fixed length scale, residual-PCA initialisation, ARD); 3. the
program test (\[rff_program_test()\], gene-shift null by default,
sequential), optionally also \[rff_factor_test()\]; 4. matching of the
two kinds of components into one Programs table.

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

  Fit settings (see \[fit_spatial_rff()\]).

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

- ...:

  Passed to \[rff_report()\] when \`report\` is given.

## Value

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
\`axis_label\` (M\* / D\* names of the report). Attributes \`loadings\`
(genes x rows, the loadings used for interpretation) and \`fields\` (all
bins x rows: the program map field, or the axis score for Candidates;
\`NA\` outside the tissue). Per-cell values: use \`rff_fields(res\$fit,
cells)\` and the \`program_map\` column (Confirmed programs), see
Examples.

## Which function should I use?

Most users call \`rff_programs()\` and then \[rff_report()\] (or
\`report = "file.html"\`). \[fit_spatial_rff()\],
\[rff_program_test()\], \[rff_factor_test()\] and \[rff_offset()\] are
the building blocks, for custom pipelines.

## Detection axes and program maps

The program test finds \*\*detection axes\*\* (principal components of
smoothed residuals, \`PC1\`, ...) and decides whether and how many
programs exist; the fit gives \*\*program maps\*\* (factors,
\`factor1\`, ...): a smooth field and gene loadings used for
interpretation. Each significant detection axis is matched to the
program map whose field correlates best with its score (greedy by \|r\|,
every map used once): \* \`Confirmed\`: significant axis matched with
\|r\| \>= \`min_r\`; p from the axis, genes, length scale and map from
the program map. \* \`Candidate\`: significant axis without such a map -
the RFLVM does not represent it well; genes and map come from the axis
itself. \* \`Not supported\`: a strong program map (large program
strength, or significant in the factor test) without a significant
matching axis. \* \`Cellularity/technical\`: program maps whose loadings
are mostly shared by all genes (uniform share \> 0.5). With only the
factor test (\`null = "none"\`, \`factor_test = TRUE\`), programs are
the significant program maps (\`detected_by = "factor test"\`); without
any test they are \`Exploratory\`.

## See also

\[rff_report()\], \[fit_spatial_rff()\], \[rff_program_test()\],
\[rff_factor_test()\], \[rff_offset()\]

## Examples

``` r
# \donttest{
tx <- simulate_transcripts(size = 100, rate = 0.02, n_genes_per_set = 3)
b <- bin_transcripts(tx, bin_size = 5)
res <- rff_programs(b, n_factors = 3, lengthscale = 10, n_boot = 19,
                    fit_args = list(max_iter = 60))
res
#> <rff_programs> 2 confirmed, 0 candidate program(s); program test, shift null, 19 draws; alpha 0.05; 0 s
#>   P1  Confirmed   p 0.05   [axis PC1; map factor1, r 0.92]  up: B_2, B_3, B_1  down: A_1, A_2, A_3, C_2, C_3, C_1
#>   P2  Confirmed   p 0.05   [axis PC2; map factor2, r 0.85]  up: C_2, C_1, C_3  down: A_3, A_1, A_2, B_1, B_3, B_2
as.data.frame(res)
#>   program    status  detected_by    p p_factor_test detection_axis program_map
#> 1      P1 Confirmed program test 0.05            NA            PC1     factor1
#> 2      P2 Confirmed program test 0.05            NA            PC2     factor2
#>           r best_map    best_r lengthscale program_strength uniform_share
#> 1 0.9234563  factor1 0.9234563          10        1.2105994   0.003451173
#> 2 0.8470046  factor2 0.8470046          10        0.6229164   0.008901251
#>   variance_share n_genes        top_up                     top_down key
#> 1      0.4028847       8 B_2, B_3, B_1 A_1, A_2, A_3, C_2, C_3, C_1  P1
#> 2      0.2473152       7 C_2, C_1, C_3 A_3, A_1, A_2, B_1, B_3, B_2  P2
#>   map_label axis_label
#> 1        M1         D1
#> 2        M2         D2
# per-cell values of the confirmed programs (via their program maps)
cells <- data.frame(x = runif(50, 0, 100), y = runif(50, 0, 100))
pm <- res$programs$program_map[res$programs$status == "Confirmed"]
if (length(pm)) head(rff_fields(res$fit, cells)[, pm, drop = FALSE])
#>      factor1     factor2
#> 1  0.3860252 -0.23832339
#> 2  0.2558472 -1.47706134
#> 3  2.1339783 -1.26043143
#> 4  0.2703143 -0.09860804
#> 5 -0.6231752 -2.58577076
#> 6 -1.4092947  0.48003989
# }
```
