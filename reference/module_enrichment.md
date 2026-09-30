# Interpret co-localization modules with gene sets

\*\*Experimental.\*\* Over-representation of each gene set (pathways,
cell-type markers, ...) in each module by the hypergeometric test,
against the genes that were clustered (\`universe\`), with
Benjamini-Hochberg correction over all module x set tests. Gene sets can
come from anywhere, e.g. \`msigdbr::msigdbr()\` or enrichR libraries.

## Usage

``` r
module_enrichment(modules, gene_sets, universe = NULL, min_overlap = 1)
```

## Arguments

- modules:

  Output of \[colocalization_modules()\] or
  \[colocalization_submodules()\] (or its \`modules\` data.frame).

- gene_sets:

  Named list of character vectors.

- universe:

  Background genes (default: the \`universe\` element of \`modules\`,
  i.e. all genes of the O/E matrix, or else all genes in the \`modules\`
  data.frame).

- min_overlap:

  Minimum overlap to report.

## Value

A data.frame: module, gene_set, overlap, module_size, set_size (within
the universe), odds_ratio, p, padj, genes.

## Examples

``` r
tx <- simulate_transcripts(size = 200, rate = 0.02, seed = 1)
b <- bin_transcripts(tx, bin_size = 4, tissue_radius = Inf)
mods <- colocalization_modules(colocalization_gene_matrix(b, radius = 12), n_modules = 3)
truth <- attr(tx, "truth")
module_enrichment(mods, split(truth$genes, truth$set_of))
#>   module gene_set overlap module_size set_size odds_ratio            p
#> 1     M1        A       5           5        5        231 0.0003330003
#> 5     M2        B       5           5        5        231 0.0003330003
#> 9     M3        C       5           5        5        231 0.0003330003
#>          padj                   genes
#> 1 0.000999001 A_2, A_4, A_1, A_3, A_5
#> 5 0.000999001 B_2, B_3, B_5, B_4, B_1
#> 9 0.000999001 C_5, C_1, C_4, C_2, C_3
```
