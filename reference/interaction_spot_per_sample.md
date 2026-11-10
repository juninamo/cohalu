# Per-sample interaction-spot summary

Apply \[search_interaction_spot()\] to each image of a Seurat object and
summarize how many connected interaction spots were detected and their
mean size, ready for group comparison.

## Usage

``` r
interaction_spot_per_sample(
  seurat_object,
  sample_key,
  group_key,
  cluster_col,
  target_cluster,
  cell_id = NULL,
  patient_key = NULL,
  radius,
  n_min,
  neighbors.k = NULL
)
```

## Arguments

- seurat_object:

  A Seurat object (lists are not supported here).

- sample_key, group_key, patient_key:

  Column names; see \[build_sample_design()\].

- cluster_col:

  Cluster column in meta.data.

- target_cluster:

  Target cluster(s) of interest.

- cell_id:

  Vector of cell IDs to include. Defaults to all cells.

- radius:

  Radius defining neighborhood.

- n_min:

  Minimum number of cells per spot.

- neighbors.k:

  Optional cap on the cells returned per cell, counting the cell itself.
  \`NULL\` (default) uses every cell within \`radius\`.

## Value

A data.frame with one row per sample carrying \`n_spots\`,
\`mean_spot_size\` (number of cells), \`n_cells\` and
\`spots_per_1k_cells\`. Raw \`n_spots\` scales with image size; compare
\`spots_per_1k_cells\` between groups when images differ in size.
Samples for which the spot search failed get \`NA\` (with a warning),
not 0.

## Examples

``` r
df <- generate_sim_groups(n_samples_per_group = 3,
                          group_close_ratio = list(case = 0.8, control = 0.2),
                          n_types = 4, n_cells = 200, max_loc = 250,
                          test_type = "distribute", distance_param = 8,
                          seed = 1)
seu <- sim_to_seurat(df)
interaction_spot_per_sample(seu, sample_key = "sample_id",
                            group_key = "group", cluster_col = "cell_type",
                            target_cluster = c("cell_type_1", "cell_type_2"),
                            radius = 15, n_min = 3)
#>           sample_id   group patient          target_cluster n_spots
#> case_1       case_1    case    <NA> cell_type_1,cell_type_2      15
#> case_2       case_2    case    <NA> cell_type_1,cell_type_2      12
#> case_3       case_3    case    <NA> cell_type_1,cell_type_2      10
#> control_1 control_1 control    <NA> cell_type_1,cell_type_2      11
#> control_2 control_2 control    <NA> cell_type_1,cell_type_2      13
#> control_3 control_3 control    <NA> cell_type_1,cell_type_2      16
#>           mean_spot_size n_cells spots_per_1k_cells
#> case_1           9.40000     200                 75
#> case_2          10.66667     200                 60
#> case_3          13.30000     200                 50
#> control_1       10.00000     200                 55
#> control_2       10.30769     200                 65
#> control_3        8.62500     200                 80
```
