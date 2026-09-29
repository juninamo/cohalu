# Fit a random-feature spatial factor model to binned transcripts

\*\*Experimental.\*\* Segmentation-free model of transcript counts. For
in-tissue bin \`b\` (centre \`u_b\`, area \`a\`) and gene \`j\`,
\$\$y\_{bj} \sim \mathrm{Poisson\\ or\\ NB}(\mu\_{bj}), \quad
\log\mu\_{bj} = \log a + \alpha_j + \sigma_0 f_0(u_b) + \sum\_{k=1}^K
L\_{jk} f_k(u_b),\$\$ where \\f_0\\ is a cellularity field shared by all
genes and \\f_1,\dots,f_K\\ are spatial factors. Each field is a
Gaussian process with RBF covariance \\\exp(-h^2/2\ell_k^2)\\
approximated by random Fourier features, \\f_k(u) =
\phi\_{\ell_k}(u)^\top\gamma_k\\ with \\\phi\_\ell(u) =
M^{-1/2}\[\cos(\Omega u/\ell), \sin(\Omega u/\ell)\]\\, \\\Omega \sim
N(0, I)\\ and \\\gamma_k \sim N(0, I)\\ (Gundersen et al. 2021, with the
coordinates as observed inputs). Loadings, weights, length scales and
dispersions are estimated jointly by MAP with closed-form gradients
(L-BFGS-B).

## Usage

``` r
fit_spatial_rff(
  binned,
  n_factors = 3,
  lengthscales = c(10, 25, 50, 100),
  density_lengthscale = 100,
  offset = c("area", "smoothed_total"),
  offset_bandwidth = 10,
  offset_genes = NULL,
  ard = 0,
  n_features = 64,
  family = c("nb", "poisson"),
  max_iter = 500,
  seed = 1,
  verbose = FALSE,
  init = NULL,
  weights = NULL,
  factor_init = c("random", "residual_pca"),
  basis = c("random", "grid"),
  learn_lengthscales = TRUE,
  l1 = 0,
  loadings = NULL
)
```

## Arguments

- binned:

  A \`binned_transcripts\` object.

- n_factors:

  Number of spatial factors \`K\`.

- lengthscales:

  Initial length scales (recycled to \`K\`); spreading them over scales
  breaks the symmetry between factors.

- density_lengthscale:

  Initial length scale of the cellularity field (\`NULL\` to omit the
  field). Ignored when \`offset = "smoothed_total"\`.

- offset:

  How overall cellularity enters the model. \`"area"\` (bin area only;
  cellularity is modelled by the latent field \\f_0\\) or
  \`"smoothed_total"\`: the log of the observed total transcript
  density, Gaussian-smoothed with bandwidth \`offset_bandwidth\`, is
  used as offset and no cellularity field is fitted. The latter
  conditions on where transcripts are, like the label permutation of
  \[nhood_enrichment()\], so that the factors describe composition only;
  it is the recommended choice for \`rff_pair_correlation(type =
  "composition")\` in group comparisons.

  A numeric matrix (bins x genes, natural log scale) can be given
  instead: a per-bin, per-gene log offset describing structure that is
  already known, e.g. from \[rff_offset()\] fitted on cell-type
  composition, domains or an embedding (PCA, Harmony, SCIGMA). The
  factors then describe only the spatially coherent variation that this
  known structure does not explain ("residual RFLVM"). Rows are all bins
  of \`binned\` or only its in-tissue bins; columns are genes (matched
  by name when named). The bin area is added internally.

- offset_bandwidth:

  Standard deviation (same unit as the coordinates) of the Gaussian
  kernel used for \`offset = "smoothed_total"\`.

- offset_genes:

  Genes whose transcripts define the smoothed total (default: all
  genes). With few cell types, conditioning on all transcripts also
  removes the co-localization of the tested sets; use reference genes
  that are not part of the tested pair (e.g. broadly expressed genes)
  instead.

- ard:

  Strength of a group penalty \\\lambda\sum_k \lVert L\_{\cdot
  k}\rVert_2\\ on the loading vector of each factor (automatic relevance
  determination): factors that are not needed shrink to (near) zero as a
  whole, so \`n_factors\` can be set generously. \`0\` (default) keeps
  the Gaussian prior only. After the fit, \`factor_strength\` gives the
  norm of each factor's loadings.

- n_features:

  Number of random frequencies \`M\` (2M features).

- family:

  \`"nb"\` (negative binomial, gene-specific dispersion) or
  \`"poisson"\`.

- max_iter:

  Maximum L-BFGS-B iterations.

- seed:

  Random seed for the random features and initialization.

- verbose:

  Print optimizer progress.

- init:

  Optional earlier fit to the same bins and genes, made with the same
  \`seed\` and \`n_features\` (hence the same random features), whose
  null-model parts (gene intercepts, dispersions, cellularity field) are
  used as starting values; the factors keep their usual initialisation.
  Not used by \[rff_factor_test()\]: in simulations, bootstrap refits
  warm-started this way gave much smaller null statistics than the
  cold-started observed fit (family-wise false positives 5 of 10 null
  tissues instead of 1 of 20), so the refits start from scratch.

- weights:

  Optional non-negative matrix (bins or in-tissue bins x genes) of
  likelihood weights, e.g. 0 for entries held out for cross-validation.
  In simulations, held-out-entry likelihood was not a reliable guide for
  choosing \`ard\`: with a weak program it preferred penalties that
  removed the program.

- factor_init:

  How the factors are initialised: \`"random"\` (default; small random
  loadings and weights) or \`"residual_pca"\`: from a PCA of the
  Gaussian-smoothed Pearson residuals under the starting null model (bin
  area, \`offset\`, gene intercepts), after removing the direction of a
  cellularity effect shared by all genes. Each score map is projected
  onto the random features. \[rff_factor_test()\] uses the same
  initialisation for its refits.

- basis:

  \`"random"\` (default): random Fourier features (\`n_features\`
  frequencies). \`"grid"\`: the deterministic Fourier basis of the
  (zero-padded) bin grid - each field is the RBF Gaussian process on a
  torus, applied by FFT (circulant embedding), with one weight per
  padded grid cell. The capacity then does not depend on \`n_features\`
  (ignored), which matters for small-scale structure over large tissue;
  fields at other points are interpolated bilinearly by
  \[rff_fields()\]. In simulations with 10 um programs, the grid basis
  with fixed length scales recovered the program map better (median
  \|r\| 0.74 vs 0.60).

- learn_lengthscales:

  If \`FALSE\`, the factor length scales stay at \`lengthscales\` (the
  cellularity length scale is still estimated). With \`basis = "grid"\`,
  learned length scales tend to shrink to the bin size (fitting noise),
  so fixed scales are recommended there.

- l1:

  L1 penalty on single loadings, \\\sum\_{jk} \lambda_j \|L\_{jk}\|\\
  (smoothed at 0), in addition to the group penalty \`ard\`: a number or
  one weight per gene (e.g. proportional to the square root of the
  gene's counts). Encourages programs that involve few genes; in
  simulations it did not improve the power of \[rff_factor_test()\].

- loadings:

  Optional numeric matrix (genes x factors, log scale; rows matched by
  gene name when named, missing genes get 0) of loadings that are held
  fixed: only the fields (and the null-model parts) are estimated, and
  \`n_factors\` is set to \`ncol(loadings)\`. Use it to map gene
  programs found elsewhere (e.g. shared loadings from
  \[rff_program_test_joint()\] or \[fit_spatial_rff_joint()\]) in a new
  tissue. The fitted field of each program is rescaled to unit variance
  as usual, so \`fit\$L\` is \`loadings\` times the program's
  \`amplitude\` in this tissue (the standard deviation of the fitted
  field, returned as \`amplitude\`).

## Value

An object of class \`spatial_rff_fit\` with the estimates (\`alpha\`,
\`L\`, \`sigma0\`, \`lengthscales\`, \`density_lengthscale\`, \`gamma\`,
\`dispersion\`), \`factor_strength\` (norm of each factor's loadings),
\`program_strength\` (norm after removing the loading shared by all
genes, which cannot be told apart from cellularity), the random
frequencies, the genes, the settings (\`settings\`, used by
\[rff_factor_test()\]) and convergence information; with fixed
\`loadings\`, also \`amplitude\` (standard deviation of each program's
fitted field; \`NULL\` otherwise).

## Details

The fields are evaluated on the bin grid through the factorisation
\\\cos(\omega_x x + \omega_y y) = \cos\omega_x x \cos\omega_y y -
\sin\omega_x x \sin\omega_y y\\, which turns the \\N \times M\\
trigonometric evaluations into small matrix products (same result,
several times faster; irregular coordinates fall back to the direct
evaluation). The optimiser usually stops at \`max_iter\` rather than at
convergence, so the loadings - and hence \`factor_strength\` - depend on
\`max_iter\`; compare fits (and bootstrap refits) only at the same
\`max_iter\`.

## References

Gundersen GW, Zhang MM, Engelhardt BE (2021). Latent variable modeling
with random features. AISTATS, PMLR 130.

## Examples

``` r
tx <- simulate_transcripts(size = 120, rate = 0.02, n_genes_per_set = 3)
b <- bin_transcripts(tx, bin_size = 6)
fit <- fit_spatial_rff(b, n_factors = 2, n_features = 32, max_iter = 50)
fit
#> <spatial_rff_fit> nb model, 9 genes, 400 bins, 2 factors
#>   factor length scales: 11.2, 20 
#>   cellularity field: sd 0.79, length scale 34.2
#>   factor strength (loading norm): 2.46, 1.8 
#>   convergence: 1 (NEW_X)
```
