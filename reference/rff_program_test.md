# Gene-program test on smoothed residuals, with the fitted RFLVM as null model

\*\*Experimental.\*\* A hybrid of residual PCA and the random-feature
model. The fit of \[fit_spatial_rff()\] supplies the null model: bin
area, \`offset\`, gene intercepts, the cellularity field and the part of
every factor that is shared by all genes (the mean of its loadings),
i.e. everything except gene programs. Pearson residuals of the counts
under this null mean (negative binomial variance with the fitted
dispersions) are Gaussian-smoothed; the per-bin direction of a
multiplicative effect shared by all genes (proportional to
\\\sqrt{\mu}\\) is projected out, genes are standardised, and the
variance shares of the principal components are compared with those of
data simulated from the null model (parametric bootstrap through the
same pipeline).

## Usage

``` r
rff_program_test(
  fit,
  binned,
  bandwidth = NULL,
  n_components = 6,
  n_boot = 99,
  sequential = TRUE,
  alpha = 0.05,
  seed = 1,
  null = c("parametric", "shift"),
  highpass = NULL
)
```

## Arguments

- fit:

  Output of \[fit_spatial_rff()\] with a matrix or \`"area"\` offset.

- binned:

  The \`binned_transcripts\` object used for the fit.

- bandwidth:

  Standard deviation of the Gaussian smoothing (coordinate units);
  default the bin size.

- n_components:

  Number of components reported (and at most tested).

- n_boot:

  Number of simulated null data sets per tested component.

- sequential:

  Sequential testing (see Details); \`FALSE\` tests every component
  against the largest variance share of the null.

- alpha:

  Level at which sequential testing stops.

- seed:

  Random seed.

- null:

  How the null distribution is generated. \`"parametric"\` (default):
  counts simulated from the fitted null model (independent noise around
  a smooth mean), as described above. It assumes that residual variation
  not shared by genes is spatially independent; when genes have their
  own spatially autocorrelated residual structure (usual in real tissue)
  it calls components in almost every tissue (100 simulated tissues with
  gene-specific residual fields) - significance then means only
  "structure the offset misses". \`"shift"\`: gene-shift surrogates that
  keep each gene's own spatial autocorrelation. Every gene's counts are
  moved together with its null mean by an independent random rigid
  transformation of the bin grid (flips, transposition for square grids
  and a shift on the mirror-extended grid, so that no seams are
  created), which destroys only the alignment between genes. Residuals
  use the Poisson variance, the shared multiplicative direction is
  removed and each bin is normalised over genes; component \`k\` is
  compared with the \`k\`-th variance share of the surrogates (parallel
  analysis) and sequential p-values are cumulative maxima. Significance
  then means a program shared by several genes. In simulations (300 x
  300 um, 30 genes) it called 0/40 tissues without any residual
  structure, 1/40, 1/40 and 4/40 with independent gene-specific residual
  fields of length scale 10, 20 and 40 um, 0/20 with extra cellularity,
  and 80 tissues with a 2-8-gene program. Gene-specific fields at scales
  larger than the program mask it (amplitude 1 program with 20 um fields
  of sd 0.25: 4/30 called); use \`highpass\` then. A covariate missing
  from the offset that moves several genes together is a shared
  structure and is called. On real Xenium windows, data with every gene
  shifted independently were called in 0/30 replicates for three windows
  but in 3/10 and 5/10 for two (an irregular tissue mask; a window
  dominated by a large lymphoid aggregate), so borderline p-values on
  real data should be read with care.

- highpass:

  Only for \`null = "shift"\`: if given (coordinate units), structure at
  scales above \`highpass\` is removed from the smoothed residuals
  (difference of Gaussians with standard deviations \`bandwidth\` and
  \`highpass\`), so that broad gene-specific residual fields do not mask
  small-scale programs. With \`highpass = 4 \* bandwidth\` (20 um),
  power for an amplitude-1 program with 20 um gene-specific fields rose
  from 4/30 to 26/30 tissues, with 0-3/30 calls without a shared
  program.

## Value

A data frame with one row per component: \`component\`, \`share\`
(variance share), \`p\`, \`factor\` and \`r_factor\` (the fitted factor
whose field correlates best with the component score, and that
correlation) and the top genes by loading. Attributes \`scores\`
(in-tissue bins x components), \`loadings\` (genes x components),
\`null_shares\` and, for \`null = "shift"\`, \`null\`.

## Details

With \`sequential = TRUE\` (default), component \`k\` is compared with
the \`k\`-th variance share of data simulated from the null model plus
the gene programs (loadings minus their mean) of the \`k - 1\` strongest
factors of the fit, so that structure already explained by earlier
components is part of the null; testing stops at the first component
with \`p \> alpha\` (p-values are made non-decreasing).

In simulations (300 x 300 um, 30 genes, 2 known domains in the offset;
fit with \`basis = "grid"\`, fixed length scale 8 um, \`ard = 100\`)
this test was as powerful as residual PCA with its own bootstrap (20/20
tissues at program amplitude 0.5, vs 12-13/20 for
\[rff_factor_test()\]), called 1/40 tissues whose only extra structure
was cellularity shared by all genes (residual PCA 12/20) and 3/60
tissues without any program.

## See also

\[fit_spatial_rff()\], \[rff_factor_test()\], \[rff_offset()\]

## Examples

``` r
# \donttest{
tx <- simulate_transcripts(size = 100, rate = 0.02, n_genes_per_set = 3)
b <- bin_transcripts(tx, bin_size = 6)
fit <- fit_spatial_rff(b, n_factors = 2, ard = 2, n_features = 24, max_iter = 40)
rff_program_test(fit, b, n_boot = 19)
#>   component      share    p  factor  r_factor
#> 1       PC1 0.54307310 0.05 factor2 0.8601089
#> 2       PC2 0.27262779 0.05 factor1 0.4313078
#> 3       PC3 0.05733013 1.00 factor1 0.2522480
#> 4       PC4 0.04647964   NA factor1 0.2007802
#> 5       PC5 0.02562339   NA factor1 0.2459856
#> 6       PC6 0.02226978   NA factor2 0.0422166
#>                                                         top_genes
#> 1 B_3 (+0.42), B_2 (+0.42), B_1 (+0.41), A_1 (-0.41), A_2 (-0.40)
#> 2 C_2 (-0.54), C_3 (-0.53), C_1 (-0.53), A_3 (+0.21), A_2 (+0.18)
#> 3 C_1 (-0.76), C_3 (+0.42), C_2 (+0.30), A_1 (+0.22), B_1 (+0.20)
#> 4 C_2 (+0.71), C_3 (-0.67), B_1 (+0.18), B_3 (-0.10), B_2 (-0.06)
#> 5 A_2 (-0.73), A_1 (+0.60), C_3 (-0.18), B_1 (-0.18), C_1 (+0.13)
#> 6 A_3 (-0.74), A_1 (+0.38), B_1 (-0.36), A_2 (+0.34), B_3 (+0.19)
# }
```
