# Local Co-occurrence Score (Generic method)

Local Co-occurrence Score (Generic method)

## Usage

``` r
cooccur_local(
  df,
  cluster_x,
  cluster_y,
  connectivity_key = "nn",
  neighbors.k = 20,
  radius = 30,
  maxnsteps = 1,
  radius_k = NULL
)
```

## Arguments

- df:

  Data.frame of coordinates and cluster.

- cluster_x:

  First cluster of interest.

- cluster_y:

  Second cluster of interest.

- connectivity_key:

  Graph type to use.

- neighbors.k:

  Number of neighbors of the diffusion graph. It no longer limits the
  radius search; see \`radius_k\`.

- radius:

  Radius for neighborhood.

- maxnsteps:

  Maximum number of diffusion steps. Each step computes \`s \<- (A + I)
  D^-1 s\` from the previous step; diffusion stops early once the
  kurtosis of the scores decreases by less than 3 between steps (checked
  after step 3). \`0\` returns the raw 0/1 indicator. Versions \<=
  0.99.1 always performed a single step regardless of \`maxnsteps\`.

- radius_k:

  Optional cap on the cells returned per cell by the radius search,
  counting the cell itself. \`NULL\` (default) uses every cell within
  \`radius\`.

## Value

Data.frame with scores.

## Examples

``` r
df <- generate_sim(close_ratio = 0.8, n_types = 4, n_cells = 300,
                   max_loc = 300, test_type = "distribute",
                   distance_param = 10, seed = 1)
sc <- cooccur_local(df, cluster_x = "cell_type_1", cluster_y = "cell_type_2",
                    neighbors.k = 10, radius = 20)
summary(sc[[1]])
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.0000  0.2912  0.6024  0.5367  0.7692  1.1405 
```
