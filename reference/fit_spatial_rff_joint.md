# Random-feature factor model across several windows with shared loadings

\*\*Experimental.\*\* Fits one set of gene loadings shared by several
tissues or windows, with a window-specific Gaussian-process field for
every program: \$\$\log\mu\_{wbj} = \log a + o\_{wbj} + \alpha\_{wj} +
\sigma\_{0w} f\_{0w}(u_b) + \sum_k L\_{jk} f\_{wk}(u_b).\$\$ Estimation
alternates between (1) fitting every window with the loadings held fixed
(\[fit_spatial_rff()\] with \`loadings\`, which estimates the fields,
intercepts and cellularity field) and (2) updating the loadings of every
gene by a ridge-penalised Poisson regression of its counts, pooled over
all windows, on the fitted program fields (with each window's null part
as offset). Loadings are rescaled to unit norm after every update; the
size of a program in a window is its \`amplitude\` (standard deviation
of its field there), a per-window measure of program activity.

## Usage

``` r
fit_spatial_rff_joint(
  binned,
  offsets = NULL,
  loadings,
  init_scale = c("residual", "log"),
  n_iter = 3,
  lengthscales = 20,
  ridge = 1,
  n_cores = 1,
  ...
)
```

## Arguments

- binned:

  List of \`binned_transcripts\` objects with the same genes.

- offsets:

  \`NULL\` (bin area only) or a list with one offset per window (as
  accepted by \[fit_spatial_rff()\]; \`NULL\` entries mean \`"area"\`).

- loadings:

  Initial loadings, genes x programs (rows matched by gene name when
  named).

- init_scale:

  \`"residual"\` (default; the loadings come from the processed
  residuals of \[rff_program_test_joint()\] and are divided by the
  square root of each gene's mean count per bin to put them on the log
  scale) or \`"log"\`.

- n_iter:

  Number of alternating rounds.

- lengthscales:

  Fixed length scales of the program fields (recycled).

- ridge:

  Ridge penalty on the loadings in the pooled update (the single-window
  model uses 1, the standard normal prior).

- n_cores:

  Cores for the per-window fits (forked with \[parallel::mclapply()\];
  serial on Windows).

- ...:

  Further arguments for \[fit_spatial_rff()\] (e.g. \`ard\`,
  \`max_iter\`, \`family\`, \`density_lengthscale\`); \`basis = "grid"\`
  and fixed length scales are used.

## Value

An object of class \`spatial_rff_joint\`: \`loadings\` (genes x
programs, unit-norm columns, log scale), \`fits\` (list of per-window
\`spatial_rff_fit\` objects with those loadings held fixed; \`fit\$L\` =
loadings x amplitude), \`amplitude\` (windows x programs), \`change\`
(largest change of the loadings in each round) and \`settings\`.

## Details

Start from directions found by \[rff_program_test_joint()\]; their
residual-scale loadings are converted to the log scale internally when
\`init_scale = "residual"\`.

## See also

\[rff_program_test_joint()\], \[fit_spatial_rff()\]

## Examples

``` r
# \donttest{
b <- lapply(1:3, function(s) bin_transcripts(
  simulate_transcripts(size = 80, rate = 0.02, n_genes_per_set = 3, seed = s), bin_size = 6))
fits <- lapply(b, fit_spatial_rff, n_factors = 2, basis = "grid", lengthscales = 10,
               learn_lengthscales = FALSE, max_iter = 30)
jt <- rff_program_test_joint(fits, b, n_boot = 9, n_components = 2)
jf <- fit_spatial_rff_joint(b, loadings = attr(jt, "loadings")[, 1, drop = FALSE],
                            n_iter = 2, lengthscales = 10, max_iter = 30)
jf$amplitude
#>              PC1
#> window1 1.262702
#> window2 1.277615
#> window3 1.375074
# }
```
