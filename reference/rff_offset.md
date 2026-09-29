# Offset for a residual random-feature model: expression explained by known structure

\*\*Experimental.\*\* Fits, for every gene, a Poisson regression of the
bin counts on covariates that describe structure already known -
cell-type composition of the bins, tissue domains (one-hot), or an
embedding such as PCA, Harmony or SCIGMA - with the bin area as offset,
and returns the fitted log expectation per bin and gene (without the bin
area). Passed as \`offset\` to \[fit_spatial_rff()\], the spatial
factors then capture only the spatially coherent variation that the
covariates do not explain.

## Usage

``` r
rff_offset(
  binned,
  covariates = NULL,
  genes = NULL,
  ridge = 1e-04,
  method = c("covariates", "kmeans", "pca"),
  k = 6,
  seed = 1,
  base = NULL,
  gene_covariates = NULL
)
```

## Arguments

- binned:

  A \`binned_transcripts\` object.

- covariates:

  Numeric matrix or data frame with one row per bin of \`binned\` (all
  bins, or only the in-tissue bins). Factors (e.g. domain or section
  labels) are expanded to indicator columns. An intercept is added.

- genes:

  Genes to model (default: all genes of \`binned\`).

- ridge:

  Small ridge penalty stabilising the per-gene fits (collinear or sparse
  covariates).

- method:

  Where the covariates come from. \`"covariates"\` (default):
  \`covariates\` as given. When no labels are available: \`"kmeans"\` -
  the locally smoothed composition of \`k\` k-means clusters of the
  bins' log-normalised expression (top 10 principal components), a
  stand-in for cell-type composition; \`"pca"\` - the top \`k\`
  principal component scores of the bins' log-normalised expression.
  Given \`covariates\` are added to the data-driven ones. In simulations
  without labels (2 domains or 4 cell-type territories, minor 4-gene
  program of amplitude 1, tested with \`rff_program_test(null =
  "shift")\`), \`"kmeans"\` with \`k = 6\` found the minor program in
  11/12 and 12/12 tissues (true labels: 12/12; no offset: 0/12), but
  only 8/12 and 1/12 at amplitude 0.5 (true labels: 10/12 and 12/12).
  PCA offsets with \`k \>= 5\` absorbed the program (0-2/12).

- k:

  Number of clusters (\`"kmeans"\`) or components (\`"pca"\`).

- seed:

  Random seed for k-means.

- base:

  Optional numeric matrix (bins or in-tissue bins x genes, natural log
  scale, per unit area like the returned offset) of expected expression
  that is already known - e.g. from \[rff_expected_offset()\], the
  expected counts of each bin given the cell types that own its
  transcripts. It enters every gene's regression as a fixed offset, so
  the covariates only adjust it; with \`covariates = NULL\` (and
  \`method = "covariates"\`) only a gene intercept is fitted. The
  returned offset is \`base\` plus the fitted adjustment.

- gene_covariates:

  Optional matrix, or list of matrices, of gene-specific covariates
  (bins or in-tissue bins x genes, columns matched by gene name when
  named): for gene \`j\`, column \`j\` of every matrix is added to its
  regression, with its own coefficient. Use it for nuisance structure
  that differs between genes, e.g. the log ratio of the counts expected
  from neighbouring cells to those expected from the owning cells
  (segmentation spill-over; \`attr(, "neighbour")\` of
  \[rff_expected_offset()\] with \`neighbour_bandwidth\`). Covariates
  shared by all genes (e.g. the per-bin share of nuclear transcripts) go
  into \`covariates\`; each gene still gets its own coefficient.

## Value

A numeric matrix (in-tissue bins x genes, natural log scale) for
\`fit_spatial_rff(offset = )\`.

## Details

An offset estimated from the same data is not neutral: the richer it is,
the more of any program it absorbs, and structure it misses is called by
\[rff_program_test()\] like any other shared structure (in the
simulations, k-means offsets left 1-3 components of the known structure
to be called in tissues without a minor program). With data-driven
offsets, report the number of clusters or components, check that
conclusions hold for neighbouring values, and judge called components by
their genes and maps.

## See also

\[fit_spatial_rff()\], \[rff_factor_test()\]

## Examples

``` r
tx <- simulate_transcripts(size = 120, rate = 0.02, n_genes_per_set = 3)
b <- bin_transcripts(tx, bin_size = 6)
# known structure: here a smooth coordinate trend as a stand-in covariate
cv <- b$coords[b$coords$in_tissue, c("x", "y")]
off <- rff_offset(b, cv)
fit <- fit_spatial_rff(b, n_factors = 3, offset = off, ard = 2,
                       n_features = 32, max_iter = 50)
fit$factor_strength
#>  factor1  factor2  factor3 
#> 1.357359 1.733957 1.128547 
```
