# Length scale of each spatial program by held-out likelihood

\*\*Experimental.\*\* \[fit_spatial_rff()\] does not estimate length
scales: learned length scales stay near their initial values, and its
MAP objective favours the smallest one (in spike-in tests the best
objective was always at the shortest length scale).
\`rff_lengthscale_profile()\` profiles each program instead. For program
\`k\` of \`fit\`, everything else in the fit (offset, gene intercepts,
cellularity field, the other factors) is held fixed as an offset, the
program's loadings are held fixed, and only its field (plus gene
intercepts and dispersions) is refitted at every length scale of
\`ls_grid\`, with the in-tissue bins split into \`folds\` groups: each
group is left out once (likelihood weight 0) and the fitted field
predicts its counts from the neighbouring bins. The profile is the
held-out log-likelihood gain over the same model without the program;
its maximum (refined by a parabola in log length scale) is the estimate.
A spatial block bootstrap of the per-bin held-out gains gives a 95
reflects the uncertainty of the held-out criterion for this field
realisation only and was much narrower than the spread between replicate
spike-ins.

## Usage

``` r
rff_lengthscale_profile(
  fit,
  binned,
  factors = NULL,
  ls_grid = 5 * 2^(0:6),
  folds = 4,
  holdout_size = NULL,
  max_iter = 150,
  n_boot = 200,
  boot_block = NULL,
  gene_share = 0.95,
  seed = 1,
  n_cores = 1,
  reach_kernel = "exponential"
)
```

## Arguments

- fit:

  Output of \[fit_spatial_rff()\] (preferably \`basis = "grid"\`).

- binned:

  The \`binned_transcripts\` object used for the fit.

- factors:

  Factors to profile (names or indices of \`fit\$L\`). Default: factors
  with program strength at least 10 loading pattern not dominated by a
  part shared by all genes (uniform share \<= 0.5).

- ls_grid:

  Length scales (coordinate units) to profile.

- folds:

  Number of held-out groups.

- holdout_size:

  Side (coordinate units) of the square units that are held out
  together; default 3 bins. For programs whose length scale is close to
  the bin size, use single bins (\`holdout_size\` = bin size): a field
  cannot predict across a block much larger than its length scale.

- max_iter:

  Iterations for each refit.

- n_boot:

  Block-bootstrap replicates for the intervals.

- boot_block:

  Side (coordinate units) of the bootstrap blocks (default 8 bins).

- gene_share:

  Only genes that carry this share of the program's squared loading norm
  (at least 5 genes) enter the refits; the others hardly inform the
  field.

- seed:

  Random seed (fold assignment, bootstrap).

- n_cores:

  Cores for the refits (forked; serial on Windows).

- reach_kernel:

  Decay kernel for the implied reach (\[rff_reach()\]):
  \`"exponential"\` (default), \`"gaussian"\`, or \`NULL\` (no reach).

## Value

An object of class \`rff_ls_profile\`: \`summary\` (one row per profiled
factor: \`factor\`, \`lengthscale\` (estimate), \`lower\`, \`upper\` (95
gain at the best grid length scale), \`gain_se\`,
\`heldout_dev_explained\` (share of the held-out deviance of the program
genes explained by the program), \`at_boundary\` (maximum at the edge of
\`ls_grid\`), \`n_genes\`, \`top_genes\`, \`reach\` / \`reach_lower\` /
\`reach_upper\` (decay length implied by the length scale,
\[rff_reach()\]) and \`identifiable\` (the length scale is at least 2
bins, at most a fifth of the smaller side of the tissue, and not at the
grid edge; otherwise the window or the bins limit it and it is a bound,
not an estimate)), \`curves\` (one row per factor and length scale:
\`gain\`, \`se\`, \`delta_se\` (standard error of the difference to the
best length scale), \`within_1se\`, \`heldout_dev_explained\`), \`boot\`
(bootstrap estimates) and \`settings\`.

## Details

Held-out units of 3 x 3 bins (default) make the field interpolate over a
short distance. Single held-out bins compressed the estimates towards
20-40 um in spike-ins (every smooth field interpolates one bin well);
larger units favour longer length scales and became unstable (estimates
at the grid edge). The other factors and the cellularity field come from
the full-data fit and therefore also saw the held-out bins; this is the
same for every length scale and does not move the maximum, but it makes
the gains conservative.

\*\*Bin size and window decide what is measured.\*\* Responses are
carried by single cells, so at fine bins (8 um) the dominant correlation
of a response map is the cell (half-correlation distance about 5 um in
spike-ins), and the profile returns 25-50 um whatever the reach. The
envelope (the reach) is seen at bins that average several cells, and the
window must be much larger than it. Compression diagnosis on spike-ins
(6-gene programs added to non-T cells of a real RA section around hidden
producers, decaying as \\e^{-d/\lambda}\\, same number of added
transcripts for every \\\lambda\\; 2 seeds): with 32-um bins in a 3.2-mm
window the profiled length scale converted by \[rff_reach()\] was within
a factor 1.5 of \\\lambda\\ for \\\lambda\\ = 40-160 um (6/6; Spearman
0.98 for \\\lambda\\ \>= 40), about half of it at \\\lambda\\ = 320
(window limit) and meaningless for \\\lambda\\ \<= 20 (below the bin
resolution). With 16-um bins in a 1.6-mm window the estimates ranked
\\\lambda\\ = 10-320 (Spearman 0.95) and ordered co-localised programs
with 4x different reach in 6/6 runs, but compressed (reach 25-95 um).
Recommended: profile at two or three bin sizes; treat a reach as
measured only when it is at least 1.5 bins and at most about a twentieth
of the window at that bin size, and otherwise report the ranking.

## When to use / limitations

Use the profile to \*\*rank\*\* the spatial scales of programs and to
say which ranges are identifiable, not to measure a signalling distance.
Report \`reach\` with \`identifiable\`, from profiles at 16-32 um bins
(at 8-um bins single cells dominate and the answer is 25-50 um for any
reach). Reaches of 20 um or less are not recoverable without known
producers; when producers are known, a producer-conditioned regression
on distance is the more direct measurement. The bootstrap interval is
too narrow to describe the reach itself (it did not cover the spread
between replicate spike-ins). \[rff_lengthscale_vi()\] gave no better
estimates on real spike-ins and is slower.

## See also

\[fit_spatial_rff()\] (\`lengthscales = "profile"\`), \[rff_programs()\]

## Examples

``` r
# \donttest{
set.seed(1); n <- 30; res <- 5
f1 <- cohalu:::.grf_fft(n, n, res, 15)
eta <- log(0.3) + outer(as.vector(f1), c(rep(0.8, 4), rep(0, 6)))
Y <- matrix(rpois(length(eta), exp(eta) * res^2), n * n, 10, dimnames = list(NULL, paste0("g", 1:10)))
gx <- rep(seq_len(n) - 1, times = n); gy <- rep(seq_len(n) - 1, each = n)
b <- structure(list(counts = Matrix::Matrix(Y, sparse = TRUE),
  coords = data.frame(x = (gx + 0.5) * res, y = (gy + 0.5) * res, in_tissue = TRUE),
  grid = list(nx = n, ny = n, bin_size = res, xmin = 0, ymin = 0), genes = colnames(Y)),
  class = "binned_transcripts")
fit <- fit_spatial_rff(b, n_factors = 1, lengthscales = 10, basis = "grid",
                       learn_lengthscales = FALSE, factor_init = "residual_pca", max_iter = 50)
pr <- rff_lengthscale_profile(fit, b, ls_grid = c(5, 10, 20, 40), folds = 3, max_iter = 50, n_boot = 50)
pr
#> <rff_ls_profile> 1 program(s); grid 5/10/20/40; 3 folds, held-out units 15
#>   factor1  length scale 10.1 (95% CI 9.15-12.2), reach 5.76 | held-out gain 2440.0 (se 321.5), dev. explained 74.69% | g4 (+0.69), g1 (+0.68), g3 (+0.68), g2 (+0.68), g10 (-0.02)
# }
```
