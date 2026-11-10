# Local co-localization score with a permutation-calibrated expectation

An abundance-adjusted replacement for the diffusion-smoothed sCLS of
\[cooccur_local()\]. For every cell \\v\\, the number of
\`cluster_x\`-\`cluster_y\` pairs among the cells within \`radius\`,
\\s_v = n_x(v)\\ n_y(v)\\, is compared with its exact expectation under
random relabelling of cells with positions fixed, \$\$E_v = N_v
(N_v - 1) \frac{K_x K_y}{M (M - 1)},\$\$ where \\N_v\\ is the number of
neighbours of \\v\\, \\M\\ the number of other cells and \\K_x, K_y\\
the numbers of cells of each type among them. Observed and expected
values are summed with Gaussian weights of width \`bandwidth\`
(micrometres) around each cell, giving the local \$\$\log_2
\mathrm{O/E}\_v = \log_2 \frac{\sum_u w\_{vu} s_u + c}{\sum_u w\_{vu}
E_u + c}.\$\$ Summed over the whole section this is the section-level
O/E, which, unlike the mean sCLS, does not grow with the abundance of
the two cell types or with cell density. Pair counts are used rather
than the 0/1 "both present" indicator of the original sCLS, because the
indicator saturates: strongly co-localized types occupy fewer
neighbourhoods than scattered ones, which would give them a negative
O/E.

## Usage

``` r
cooccur_local_oe(
  df,
  cluster_x,
  cluster_y,
  radius = 30,
  neighbors.k = NULL,
  bandwidth = radius,
  cluster_col = "cell_type",
  n_perms = 0,
  pseudocount = 0.5,
  seed = 1
)
```

## Arguments

- df:

  A data.frame with \`x\`, \`y\` and a cell-type column; row names are
  used as cell IDs.

- cluster_x, cluster_y:

  The two cell types (may be identical).

- radius:

  Neighbourhood radius for the indicator.

- neighbors.k:

  Optional cap on the neighbours used per cell, excluding the cell
  itself. \`NULL\` (default) uses every cell within \`radius\`. A number
  truncates dense neighbourhoods to the nearest cells, and a warning
  reports how many cells reached it.

- bandwidth:

  Width of the Gaussian smoothing window (the kernel SD is \`bandwidth /
  2\`, truncated at \`bandwidth\`). Defaults to \`radius\`.

- cluster_col:

  Column holding the cell types.

- n_perms:

  Number of label permutations for per-cell hotspot p-values (0 = none).
  Each permutation costs one pass over the neighbour index.

- pseudocount:

  Added to smoothed observed and expected pair counts.

- seed:

  Random seed for the permutations.

## Value

A data.frame (one row per cell) with \`n_x\`, \`n_y\`, \`pairs\`,
\`expected\`, \`n_neighbours\`, \`local_log2_oe\` and, if \`n_perms \>
0\`, one-sided \`p\` (enrichment) and BH-adjusted \`padj\`. The
attribute \`section_log2_oe\` holds log2(total pairs / total expected
pairs).

## Examples

``` r
df <- generate_sim(close_ratio = 1, n_types = 6, n_cells = 600,
                   max_loc = 400, test_type = "circle",
                   distance_param = 15, seed = 3)
rownames(df) <- paste0("cell", seq_len(nrow(df)))
res <- cooccur_local_oe(df, "cell_type_1", "cell_type_2", radius = 25,
                        n_perms = 99)
attr(res, "section_log2_oe")
#> [1] 1.375594
head(res)
#>       n_x n_y pairs  expected n_neighbours local_log2_oe    p       padj
#> cell1  11   0     0  6.633129           16    -6.4151188 1.00 1.00000000
#> cell2  10   5    50 16.582822           25     1.3740289 0.02 0.06030151
#> cell3  11   0     0  8.457239           18    -1.1325468 0.90 1.00000000
#> cell4  15   0     0 12.768773           22    -0.7459833 0.74 1.00000000
#> cell5  11   1    11  7.517546           17     1.0811611 0.02 0.06030151
#> cell6  14   1    14 13.984847           23     1.3343143 0.01 0.03508772
```
