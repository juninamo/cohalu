# Search for Spatial Interaction Spots

Search for Spatial Interaction Spots

## Usage

``` r
search_interaction_spot(
  seurat_object,
  fov,
  radius,
  n_min,
  neighbors.k = NULL,
  cell_id = cell_id,
  cluster_col = cluster_col,
  target_cluster = target_cluster
)
```

## Arguments

- seurat_object:

  Seurat object.

- fov:

  Field of view identifier.

- radius:

  Radius threshold for neighborhood.

- n_min:

  Minimum number of cells to qualify as a spot.

- neighbors.k:

  Optional cap on the cells returned per cell, counting the cell itself.
  \`NULL\` (default) uses every cell within \`radius\`.

- cell_id:

  Vector of target cell IDs.

- cluster_col:

  Column name in metadata specifying cluster assignment.

- target_cluster:

  Target cluster(s) assumed to be interacting.

## Value

A data.frame of detected interaction clusters and metadata.

## Examples

``` r
df <- generate_sim(close_ratio = 0.8, n_types = 4, n_cells = 300,
                   max_loc = 300, test_type = "distribute",
                   distance_param = 10, seed = 1)
df$sample_id <- "fov1"
seu <- sim_to_seurat(df)
spots <- search_interaction_spot(seu, fov = "fov1", radius = 15, n_min = 3,
                                 cell_id = seu$cell,
                                 cluster_col = "cell_type",
                                 target_cluster = c("cell_type_1", "cell_type_2"))
length(unique(spots$cluster_id))
#> [1] 21
```
