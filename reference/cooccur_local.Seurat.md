# Local Co-occurrence Score (Seurat Method)

Local Co-occurrence Score (Seurat Method)

## Usage

``` r
cooccur_local.Seurat(
  seurat_obj,
  cluster_x,
  cluster_y,
  connectivity_key = "nn",
  cluster_key = "seurat_clusters",
  sample_key = "sample_id",
  neighbors.k = 20,
  radius = 30,
  maxnsteps = 15,
  radius_k = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object.

- cluster_x:

  First cluster of interest.

- cluster_y:

  Second cluster of interest.

- connectivity_key:

  Graph type to use.

- cluster_key:

  Metadata column with cluster info.

- sample_key:

  Metadata column with sample ID.

- neighbors.k:

  Number of neighbors of the diffusion graph. It no longer limits the
  radius search; see \`radius_k\`.

- radius:

  Radius for proximity-based interaction.

- maxnsteps:

  Maximum number of diffusion steps (each step starts from the previous
  one; early stop on the kurtosis criterion, see \[cooccur_local()\]).
  \`0\` returns the raw indicator.

- radius_k:

  Optional cap on the cells returned per cell by the radius search,
  counting the cell itself. \`NULL\` (default) uses every cell within
  \`radius\`.

## Value

A data.frame with local co-occurrence scores.

## Examples

``` r
df <- generate_sim(close_ratio = 0.8, n_types = 4, n_cells = 300,
                   max_loc = 300, test_type = "distribute",
                   distance_param = 10, seed = 1)
df$sample_id <- "fov1"
seu <- sim_to_seurat(df)
sc <- cooccur_local.Seurat(seu, cluster_x = "cell_type_1",
                           cluster_y = "cell_type_2",
                           cluster_key = "cell_type", sample_key = "sample_id",
                           neighbors.k = 10, radius = 20, maxnsteps = 1)
summary(sc[[1]])
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.0000  0.2912  0.6024  0.5367  0.7692  1.1405 
```
