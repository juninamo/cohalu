# Calculate Co-occurrence Matrix for a Given Radius

Calculate Co-occurrence Matrix for a Given Radius

## Usage

``` r
calc_co_occurrence_for_radius(
  seurat_obj,
  radius,
  sample_key,
  cluster_key,
  k = NULL
)
```

## Arguments

- seurat_obj:

  Seurat object with spatial coordinates.

- radius:

  Radius to define local neighborhood.

- sample_key:

  Metadata column specifying sample identity.

- cluster_key:

  Metadata column specifying cluster labels.

- k:

  Optional cap on the cells returned per cell, counting the cell itself.
  \`NULL\` (default) uses every cell within \`radius\`. A number keeps
  only the \`k\` nearest cells, so the result stops changing with
  \`radius\` once the cap is reached; a warning reports how many cells
  reached it.

## Value

A list with co-occurrence count and enrichment ratio matrices.

## Examples

``` r
df <- generate_sim(close_ratio = 0.8, n_types = 4, n_cells = 300,
                   max_loc = 300, test_type = "distribute",
                   distance_param = 10, seed = 1)
df$sample_id <- "fov1"
seu <- sim_to_seurat(df)
res <- calc_co_occurrence_for_radius(seu, radius = 20,
                                     sample_key = "sample_id",
                                     cluster_key = "cell_type")
round(res$ratio_mat, 2)
#>                    Clustercell_type_1 Clustercell_type_2 Clustercell_type_3
#> Clustercell_type_1               0.79               1.10               1.15
#> Clustercell_type_2               1.10               1.10               0.86
#> Clustercell_type_3               1.15               0.86               1.08
#> Clustercell_type_4               0.93               0.89               0.94
#>                    Clustercell_type_4
#> Clustercell_type_1               0.93
#> Clustercell_type_2               0.89
#> Clustercell_type_3               0.94
#> Clustercell_type_4               1.33
```
