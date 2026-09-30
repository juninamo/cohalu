# Cluster genes into co-localization modules

\*\*Experimental.\*\* Hierarchical clustering of the gene x gene log2
O/E matrix from \[colocalization_gene_matrix()\].

## Usage

``` r
colocalization_modules(
  M,
  n_modules = NULL,
  min_oe = 0.5,
  min_size = 3,
  method = c("average", "profile"),
  k_range = 2:8
)
```

## Arguments

- M:

  Output of \[colocalization_gene_matrix()\].

- n_modules:

  Number of modules (\`cutree(k = )\`). If \`NULL\`: with \`method =
  "average"\` the tree is cut at height \`max(M) - min_oe\`, i.e. genes
  join a module when their average mutual log2 O/E exceeds \`min_oe\`;
  with \`method = "profile"\` the number of modules in \`k_range\` with
  the highest mean silhouette width is used.

- min_oe:

  Threshold used when \`n_modules\` is \`NULL\` and \`method =
  "average"\`.

- min_size:

  Modules with fewer genes are labelled \`NA\`.

- method:

  \`"average"\` or \`"profile"\`, see Details.

- k_range:

  Candidate numbers of modules for \`method = "profile"\` with
  \`n_modules = NULL\`.

## Value

A list with \`modules\` (data.frame: gene, module, connectivity = mean
log2 O/E with the other genes of its module), \`summary\` (data.frame:
module, size, mean_oe, top_genes), the \`tree\`, the \`method\`,
\`universe\` (all genes of \`M\`) and, for \`method = "profile"\`,
\`silhouette\` (data.frame over \`k_range\`: k, mean_silhouette,
within_oe, between_oe, sizes).

## Details

\* \`method = "average"\` (default): genes whose transcripts lie near
one another (high mutual O/E) end up in the same module. The distance
between genes is \`max(M) - M\` with average linkage. Large modules of
co-localizing populations can chain: lowering the cut then only peels
off single genes (see \[colocalization_submodules()\]). \* \`method =
"profile"\`: genes are grouped by the similarity of their
co-localization \*profiles\* (their rows of \`M\`, i.e. which genes they
co-localize with), with distance \`1 - cor\` and Ward (\`ward.D2\`)
linkage. This separates populations that share a niche but are not
intermixed. Every gene is assigned (also genes that co-localize with
nothing), so \`mean_oe\` in the summary tells co-localizing modules from
groups of background genes.

## Examples

``` r
tx <- simulate_transcripts(size = 200, rate = 0.02, seed = 1)
b <- bin_transcripts(tx, bin_size = 4, tissue_radius = Inf)
M <- colocalization_gene_matrix(b, radius = 12)
mods <- colocalization_modules(M, n_modules = 3)
mods$summary
#>   module size     mean_oe               top_genes
#> 1     M1    5 0.785425551 A_2, A_4, A_1, A_3, A_5
#> 2     M2    5 0.324309289 B_2, B_3, B_5, B_4, B_1
#> 3     M3    5 0.001650131 C_5, C_1, C_4, C_2, C_3
colocalization_modules(M, method = "profile", n_modules = 3)$summary
#>   module size     mean_oe               top_genes
#> 1     M1    5 0.785425551 A_2, A_4, A_1, A_3, A_5
#> 2     M2    5 0.324309289 B_2, B_3, B_5, B_4, B_1
#> 3     M3    5 0.001650131 C_5, C_1, C_4, C_2, C_3
```
