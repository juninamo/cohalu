# Gene-program test across several tissues or windows with shared loadings

\*\*Experimental.\*\* Multi-window version of \`rff_program_test(null =
"shift")\`. Each window (for example one window per tertiary lymphoid
structure) has its own \[fit_spatial_rff()\] fit, which supplies its
null mean (offset, gene intercepts, cellularity field and the all-gene
part of the factors). The processed residuals of every window (Poisson
Pearson residuals, Gaussian smoothing, removal of the shared
multiplicative direction, per-bin normalisation over genes, gene
standardisation within the window) are pooled into one gene x gene
covariance, \`sum_w c_w R_w' R_w\`, whose eigenvectors are gene programs
with \*\*shared loadings\*\*; each window keeps its own score map.
Recurrent weak programs thus gain power from all windows, and no
clustering of per-window components is needed.

## Usage

``` r
rff_program_test_joint(
  fits,
  binned,
  loadings = NULL,
  bandwidth = NULL,
  highpass = NULL,
  n_components = 6,
  n_boot = 99,
  sequential = TRUE,
  alpha = 0.05,
  seed = 1,
  weights = c("equal", "bins"),
  n_cores = 1,
  crossfit = FALSE,
  crossfit_block = NULL,
  control = NULL,
  n_group_null = 499
)
```

## Arguments

- fits:

  List of \[fit_spatial_rff()\] fits (matrix or \`"area"\` offset), all
  with the same genes in the same order.

- binned:

  List of the matching \`binned_transcripts\` objects.

- loadings:

  \`NULL\` (discovery) or a numeric matrix, genes x programs (rows
  matched by gene name when named; missing genes get 0), on the scale of
  the processed residuals, e.g. the \`loadings\` attribute of an earlier
  call (confirmatory test).

- bandwidth, highpass:

  As in \[rff_program_test()\].

- n_components:

  Number of components reported and tested (discovery).

- n_boot:

  Number of surrogate data sets (each one shifts every window).

- sequential, alpha:

  Sequential testing as in \[rff_program_test()\] (discovery).

- seed:

  Random seed.

- weights:

  \`"equal"\` (default): every window contributes the same total weight
  to the pooled covariance (\`c_w = 1 / N_w\`, \`N_w\` its in-tissue
  bins), so that a few large windows do not dominate; \`"bins"\`: every
  bin counts equally. A numeric vector (one value per window) gives
  relative window weights on top of \`"equal"\`, e.g. \`1 / (number of
  windows of the patient)\` so that every patient counts the same.

- n_cores:

  Cores for the surrogates (forked with \[parallel::mclapply()\]; serial
  on Windows). Each surrogate data set has its own seed, so results do
  not depend on \`n_cores\`.

- crossfit:

  Discovery only. \`TRUE\`: cross-fitted version (as
  \`rff_program_test(crossfit = TRUE)\`). The bins of every window are
  split into two halves (checkerboard of \`crossfit_block\` blocks);
  components are found on the pooled covariance of one half and their
  pooled variance share is measured on the other half, both ways; the
  same pipeline on gene-shift surrogates gives \`p\`; \`excess\` is the
  weighted mean over windows of each window's held-out excess (with a
  window-bootstrap interval, \`excess_lower\` / \`excess_upper\`);
  \`call\` is the decision. Robust residual processing is used.

- crossfit_block:

  Side of the checkerboard blocks (default 12 bins or that of
  \`control\`).

- control:

  Output of \[rff_control_reference()\] for negative-control windows
  (implies \`crossfit = TRUE\`): \`excess_vs_control\` is the median
  over the windows of each window's held-out excess along the component
  direction, standardised by the mean and standard deviation of the
  control windows' excess along the same direction (robust to single
  windows with idiosyncratic structure), and \`p_control\` compares it
  with \`n_group_null\` pseudo-target groups of control windows (same
  number of windows, tested against the remaining controls in the same
  way); \`call\` also requires \`p_control \<= alpha\`.

- n_group_null:

  Number of pseudo-target groups for \`p_control\`.

## Value

Discovery: a data frame with one row per component (\`component\`,
\`share\`, \`p\`, \`top_genes\`) and attributes \`loadings\` (genes x
components, unit norm), \`scores\` (list with one in-tissue bins x
components matrix per window), \`window_share\` (windows x components:
share of each window's residual variance along the component),
\`null_shares\` and \`null = "shift_joint"\`. Confirmatory: a list with
\`windows\` (one row per window x program: \`window\`, \`program\`,
\`share\`, \`null_mean\`, \`null_sd\`, \`z\`, \`p\`), \`combined\` (one
row per program: \`share\` (weighted mean over windows), \`null_mean\`,
\`z\`, \`p\`, \`n_windows\`) and \`scores\` (per-window score maps of
the programs).

## Details

The null distribution comes from gene-shift surrogates generated
independently in every window (each gene's counts and null mean moved by
an independent rigid transformation of that window's mirror-extended bin
grid), which keep every gene's own spatial structure but destroy the
alignment between genes. Component \`k\` is compared with the \`k\`-th
variance share of the pooled surrogate covariance (parallel analysis;
sequential p-values are cumulative maxima).

With \`loadings\` given, the test is confirmatory instead: for every
window and program, the statistic is the share of the window's processed
residual variance along the (unit-norm) program direction, compared with
the same share in that window's surrogates (\`p\`, and \`z\` =
(observed - surrogate mean) / surrogate sd, a calibrated measure of how
strongly the program is spatially patterned there). The window
statistics are also summed (with the window weights) and compared with
the sum over windows of the surrogates (\`combined\`). Use it to
validate programs found in other data (held-out patients) and to score
program activity per window.

## See also

\[rff_program_test()\], \[fit_spatial_rff_joint()\],
\[fit_spatial_rff()\]

## Examples

``` r
# \donttest{
b <- lapply(1:3, function(s) bin_transcripts(
  simulate_transcripts(size = 80, rate = 0.02, n_genes_per_set = 3, seed = s), bin_size = 6))
fits <- lapply(b, fit_spatial_rff, n_factors = 2, basis = "grid", lengthscales = 10,
               learn_lengthscales = FALSE, max_iter = 30)
jt <- rff_program_test_joint(fits, b, n_boot = 19, n_components = 3)
jt
#>   component      share    p
#> 1       PC1 0.39447651 0.05
#> 2       PC2 0.35148244 0.05
#> 3       PC3 0.06122777 1.00
#>                                                         top_genes
#> 1 A_2 (+0.46), A_1 (+0.45), A_3 (+0.44), C_3 (-0.33), C_2 (-0.33)
#> 2 B_2 (+0.47), B_3 (+0.45), B_1 (+0.44), C_3 (-0.35), C_2 (-0.35)
#> 3 C_1 (+0.51), B_1 (+0.47), A_1 (-0.44), C_2 (-0.35), B_3 (-0.30)
# confirmatory scoring of the first component in each window
rff_program_test_joint(fits, b, loadings = attr(jt, "loadings")[, 1, drop = FALSE],
                       n_boot = 19)$windows
#>    window program     share null_mean    null_sd        z    p
#> 1 window1     PC1 0.3817189 0.1269353 0.03209430 7.938593 0.05
#> 2 window2     PC1 0.4115511 0.1121875 0.03288952 9.102097 0.05
#> 3 window3     PC1 0.3901595 0.1144506 0.03332763 8.272682 0.05
# }
```
