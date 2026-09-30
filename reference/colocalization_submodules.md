# Split a co-localization module into sub-modules

\*\*Experimental.\*\* Large modules from \[colocalization_modules()\]
(average linkage) can hold several populations that share a niche: the
tree then chains, and cutting it lower only peels off single genes. This
function splits one module by the similarity of its genes'
co-localization profiles: each gene's row of \`M\` over \*all\* genes
(which genes it co-localizes with), distance \`1 - cor\`, Ward
(\`ward.D2\`) linkage. The number of sub-modules is chosen by the mean
silhouette width over \`k_range\` unless \`k\` is given. A mean
silhouette above about 0.25-0.5 indicates real sub-structure; compare
\`within_oe\` with \`between_oe\` in the \`silhouette\` table to see
whether sub-modules co-localize more within than between.

## Usage

``` r
colocalization_submodules(M, modules, module = "M1", k = NULL, k_range = 2:8)
```

## Arguments

- M:

  Output of \[colocalization_gene_matrix()\].

- modules:

  Output of \[colocalization_modules()\] (or its \`modules\`
  data.frame), or a character vector of the genes to split.

- module:

  Name of the module to split (ignored when \`modules\` is a gene
  vector).

- k:

  Number of sub-modules; \`NULL\` chooses it by mean silhouette.

- k_range:

  Candidate numbers of sub-modules when \`k\` is \`NULL\`.

## Value

A list with \* \`modules\`: data.frame gene, module (\`"\<module\>.1"\`,
\`"\<module\>.2"\`, ... by decreasing size), parent, connectivity (mean
log2 O/E with the other genes of its sub-module), silhouette (per gene);
\* \`summary\`: module, size, mean_oe (within the sub-module),
mean_silhouette, top_genes (by connectivity); \* \`silhouette\`:
data.frame over the candidate k: k, mean_silhouette, within_oe and
between_oe (mean log2 O/E of gene pairs within / between sub-modules),
sizes; \* \`k\`, \`tree\`, \`method = "profile"\`, \`parent\` and
\`universe\` (all genes of \`M\`, the default background of
\[module_enrichment()\]).

## Examples

``` r
tx <- simulate_transcripts(gene_sets = c("A", "B", "C"), size = 200,
                           coloc = c(A = 1.5, B = 1.5), set_sd = 0.6, seed = 1)
b <- bin_transcripts(tx, bin_size = 4, tissue_radius = Inf)
M <- colocalization_gene_matrix(b, radius = 12)
mods <- colocalization_modules(M, n_modules = 2)
sub <- colocalization_submodules(M, mods, module = "M1")
sub$silhouette
#>   k mean_silhouette within_oe between_oe           sizes
#> 1 2      0.99613075 0.4843502  0.2391916             5/5
#> 2 3      0.80291015 0.4676960  0.2821951           5/4/1
#> 3 4      0.53428912 0.4690667  0.3041816         4/4/1/1
#> 4 5      0.37288531 0.4478497  0.3232263       4/3/1/1/1
#> 5 6      0.25396227 0.4430232  0.3306745     4/2/1/1/1/1
#> 6 7      0.08848237 0.4463633  0.3411358   2/2/2/1/1/1/1
#> 7 8      0.06483124 0.4431173  0.3437339 2/2/1/1/1/1/1/1
sub$summary
#>   module size   mean_oe mean_silhouette               top_genes
#> 1   M1.1    5 0.5184361       0.9941057 A_1, A_4, A_3, A_2, A_5
#> 2   M1.2    5 0.4502643       0.9981558 B_4, B_1, B_5, B_3, B_2
```
