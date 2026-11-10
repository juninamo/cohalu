# Per-sample local co-occurrence summary

Run \[cooccur_local()\] for one (cluster_x, cluster_y) pair on each
sample and summarize the per-cell scores at the sample level. Output is
a tidy data.frame suitable for group comparison.

## Usage

``` r
cooccur_local_per_sample(
  obj,
  sample_key,
  group_key,
  cluster_key,
  cluster_x,
  cluster_y,
  patient_key = NULL,
  unit = c("image", "patient"),
  neighbors.k = 20,
  radius = 30,
  maxnsteps = 1,
  summarize = c("mean", "q90", "pos_rate", "log2_oe"),
  radius_k = NULL
)
```

## Arguments

- obj:

  A Seurat object, a list of Seurat objects, or a data.frame.

- sample_key, group_key, cluster_key, patient_key:

  Sample / group / cluster / patient column names. See
  \[build_sample_design()\].

- cluster_x:

  First cluster of interest.

- cluster_y:

  Second cluster of interest.

- unit:

  "image" or "patient" (averaged across images of the same patient).

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

- summarize:

  Character vector of summary statistics to compute: "mean", "q90" (90th
  percentile), "pos_rate" (fraction of cells with score \> 0), and / or
  "log2_oe": log2 of the number of \`cluster_x\`- \`cluster_y\` pairs
  within \`radius\` of each cell, summed over cells, over its
  expectation under label permutation (see \[cooccur_local_oe()\]).
  \`log2_oe\` is adjusted for the abundance of the two cell types and
  for cell density and is the recommended summary for group comparison;
  "mean" is unchanged by the (mass-conserving) diffusion and grows with
  abundance.

- radius_k:

  Optional cap on the cells returned per cell by the radius search,
  counting the cell itself. \`NULL\` (default) uses every cell within
  \`radius\`.

## Value

A data.frame with one row per sample carrying the requested summary
statistics, plus \`n_cells\`, \`n_i\` and \`n_j\` (cells of
\`cluster_x\` / \`cluster_y\`).

## Caution

The local score has no permutation null, so its sample-level summaries
increase with the abundance of \`cluster_x\` and \`cluster_y\`. When the
two groups differ in cell-type composition, check \`n_i\` / \`n_j\` (or
adjust for them with \`covariates\` in \[compare_groups()\]) before
interpreting a group difference as a change in co-localization.

## Examples

``` r
df <- generate_sim_groups(n_samples_per_group = 3,
                          group_close_ratio = list(case = 0.8, control = 0.2),
                          n_types = 4, n_cells = 200, max_loc = 250,
                          test_type = "distribute", distance_param = 8,
                          seed = 1)
cooccur_local_per_sample(df, sample_key = "sample_id", group_key = "group",
                         cluster_key = "cell_type",
                         cluster_x = "cell_type_1", cluster_y = "cell_type_2",
                         patient_key = "patient", neighbors.k = 10,
                         radius = 20)
#>   sample_id   cluster_i   cluster_j  mean       q90 pos_rate   log2_oe n_cells
#> 1    case_1 cell_type_1 cell_type_2 0.515 0.8569268    0.920 0.8242564     200
#> 2    case_2 cell_type_1 cell_type_2 0.490 0.8304980    0.905 0.8304952     200
#> 3    case_3 cell_type_1 cell_type_2 0.595 0.8382076    0.995 0.5034124     200
#> 4 control_1 cell_type_1 cell_type_2 0.345 0.7153114    0.915 0.1183775     200
#> 5 control_2 cell_type_1 cell_type_2 0.410 0.7322789    0.875 0.3825044     200
#> 6 control_3 cell_type_1 cell_type_2 0.345 0.7429332    0.860 0.3201779     200
#>   n_i n_j   group   patient
#> 1  41  63    case    case_1
#> 2  44  48    case    case_2
#> 3  49  55    case    case_3
#> 4  47  46 control control_1
#> 5  46  55 control control_2
#> 6  46  40 control control_3
```
