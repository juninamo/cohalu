# Offset from the cell types that own each bin's transcripts

\*\*Experimental.\*\* Known structure for a residual model at cell
resolution. Every transcript of a \*reference\* gene (by default every
gene that is not modelled) is attributed to the cell type of the cell
that owns it (segmentation), or to its own class for unassigned or
unlabelled transcripts. For modelled gene \`j\` and bin \`b\`, the
expected count is \\E\_{bj} = \sum_t n\_{bt}\rho\_{jt}\\, with
\\n\_{bt}\\ the reference transcripts of type \`t\` in the bin and
\\\rho\_{jt}\\ the ratio of gene \`j\` to reference transcripts in cells
of type \`t\` (from \`profiles\`, e.g. a whole section, or from \`tx\`).
Unlike a smoothed cell-type composition, this follows the actual cells,
their size and their transcript capture, so that residuals describe
variation within cell types (cell states, programs); modelled genes
never enter their own expectation.

## Usage

``` r
rff_expected_offset(
  binned,
  tx,
  type_col = "type",
  x_col = "x",
  y_col = "y",
  gene_col = "gene",
  reference_genes = NULL,
  profiles = NULL,
  pseudo = 0.05
)
```

## Arguments

- binned:

  A \`binned_transcripts\` object (the modelled genes are
  \`binned\$genes\`).

- tx:

  Data frame of transcripts in the binned region, with coordinates, gene
  and the owner's cell type (\`NA\` = unlabelled).

- type_col, x_col, y_col, gene_col:

  Column names in \`tx\`.

- reference_genes:

  Genes whose transcripts measure the local amount of each cell type
  (default: all genes in \`tx\` or \`profiles\` that are not modelled).

- profiles:

  Optional numeric matrix, genes (rows) x cell types (columns), of
  transcript counts per type from a larger region (e.g. the whole
  section), used for \\\rho\\; default: computed from \`tx\`.

- pseudo:

  Added to the expected counts, as a fraction of each gene's mean
  expected count per in-tissue bin, to keep the log finite.

## Value

A numeric matrix (all bins x modelled genes, natural log scale per unit
area) for \`fit_spatial_rff(offset = )\` or \`rff_offset(base = )\`,
with attribute \`type_counts\` (bins x types reference counts).

## See also

\[rff_offset()\], \[fit_spatial_rff()\]

## Examples

``` r
tx <- simulate_transcripts(size = 80, rate = 0.03, n_genes_per_set = 3, seed = 1)
tx$type <- ifelse(tx$x < 40, "left", "right")      # stand-in for cell-type ownership
b <- bin_transcripts(tx, bin_size = 8)
b2 <- b; b2$genes <- b$genes[1:4]; b2$counts <- b$counts[, 1:4]
off <- rff_expected_offset(b2, tx)
dim(off)
#> [1] 100   4
```
