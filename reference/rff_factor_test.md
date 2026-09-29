# Significance of residual spatial factors by parametric bootstrap

\*\*Experimental.\*\* Tests whether the factors of a
\[fit_spatial_rff()\] fit are gene programs stronger than those fitted
to data without any factor structure. Counts are simulated from the
fitted null model (bin area, offset, gene intercepts, cellularity field
and dispersion, but no factors), the same model is refitted to each
simulated data set, and the largest factor statistic is recorded. Each
observed factor is compared with this max-null distribution, which
controls the family-wise error over factors.

## Usage

``` r
rff_factor_test(
  fit,
  binned,
  n_boot = 19,
  max_iter = NULL,
  statistic = c("program", "strength"),
  seed = 1,
  n_cores = 1,
  sequential = FALSE,
  alpha = 0.05
)
```

## Arguments

- fit:

  Output of \[fit_spatial_rff()\].

- binned:

  The \`binned_transcripts\` object used for the fit.

- n_boot:

  Number of simulated null data sets.

- max_iter:

  Iterations for the refits (default: those of the fit). The fits are
  usually stopped by \`max_iter\`, and the factor statistics grow with
  the number of iterations, so a different value gives a miscalibrated
  test (a warning is issued).

- statistic:

  \`"program"\` (default; loading norm after removing the mean over
  genes) or \`"strength"\` (raw loading norm).

- seed:

  Random seed.

- n_cores:

  Number of cores for the refits (forked with \[parallel::mclapply()\];
  serial on Windows). The null data sets are simulated before the
  refits, so results do not depend on \`n_cores\`.

- sequential:

  If \`TRUE\`, factors are also tested sequentially: the factor with the
  \`k\`-th largest statistic is compared with the \`k\`-th largest
  statistic of refits to data simulated from the null model plus the
  \`k - 1\` stronger factors of the fit (so that structure already
  explained is part of the null). Testing stops at the first factor with
  \`p \> alpha\`; the p-values (\`p_sequential\`, non-decreasing) are
  \`NA\` after that. Costs up to one bootstrap per tested factor.

- alpha:

  Level at which sequential testing stops.

## Value

A data frame with one row per factor: \`factor\`, \`strength\` (raw
loading norm), \`program_strength\`, \`uniform_share\` (share of the
squared loading norm that is common to all genes; near 1 =
cellularity-like), \`lengthscale\`, \`p\` (share of null data sets whose
largest statistic is at least as large, with the +1 correction), and the
top genes by absolute loading, and \`p_sequential\` when \`sequential =
TRUE\`. The null maxima are attached as attribute \`null_max\` (and the
sequential null statistics as \`null_sequential\`).

## Details

The default statistic is the program strength: the norm of a factor's
loadings after removing their mean over genes. A factor that moves all
genes together is extra cellularity (e.g. fine-scale density that the
smooth cellularity field cannot follow), not a gene program; such
factors are reported (\`uniform_share\`) but not declared significant.
In simulations the raw loading norm (\`statistic = "strength"\`) called
such factors significant in data without any program.

## See also

\[fit_spatial_rff()\], \[rff_offset()\], \[rff_fields()\]

## Examples

``` r
# \donttest{
tx <- simulate_transcripts(size = 100, rate = 0.02, n_genes_per_set = 3)
b <- bin_transcripts(tx, bin_size = 6)
fit <- fit_spatial_rff(b, n_factors = 2, ard = 2, n_features = 24, max_iter = 40)
rff_factor_test(fit, b, n_boot = 4)
#>    factor strength program_strength uniform_share lengthscale   p
#> 2 factor2 1.479139        1.2516185     0.2839791    21.97312 0.2
#> 1 factor1 1.773590        0.3997862     0.9491900    10.17890 0.6
#>                                                         top_genes
#> 2 B_2 (+0.94), B_3 (+0.82), B_1 (+0.69), C_1 (+0.20), A_1 (-0.19)
#> 1 A_3 (-0.71), B_3 (-0.70), A_1 (-0.66), A_2 (-0.65), B_1 (-0.63)
# }
```
