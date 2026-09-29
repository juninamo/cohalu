# Algorithm reference: core COHALU functions

This vignette documents the mathematics behind each core function in
**COHALU**. The formulas mirror the implementation in
`R/cohalu_functions.R` line-by-line, not a stylized version.

``` r

library(cohalu)
#> Loading cohalu v0.99.3: COHALU - CO-localization, Hotspots And sample-Level Units (formerly spatialCooccur).
#> To cite this package in publications, please use:
#>   Inamo J, et al. (2026). Spatial transcriptomics reveals  immune-stromal crosstalk within the synovium of patients with  juvenile idiopathic arthritis. JCI Insight. 11(1):e198074.  doi:10.1172/jci.insight.198074
#> Developed by: Jun Inamo <juninamo@keio.jp>
```

## Notation

Throughout, a sample is a finite set of cells

``` math
\mathcal{V} = \{v_1, \dots, v_N\}, \quad
\mathbf{x}_v \in \mathbb{R}^{2}, \quad
c(v) \in \mathcal{C} = \{c_1, \dots, c_K\},
```

where $`\mathbf{x}_v`$ is the 2-D centroid and $`c(v)`$ the discrete
cluster label (e.g. cell type). For each function we additionally use:

- $`A \in \{0, 1\}^{N \times N}`$ — adjacency of a
  $`k`$-nearest-neighbor (kNN) graph on $`\{\mathbf{x}_v\}`$, with
  $`A_{uv} = 1`$ iff $`v`$ is one of the $`k`$ nearest neighbors of
  $`u`$ (the cell itself included). Built exactly with a kd-tree
  (`RANN::nn2(coords, k = k)`).
- $`d_v = \sum_u A_{uv}`$ — the (in-)degree of $`v`$ in the kNN graph.
- $`N_r(v) = \{u \in \mathcal{V}\setminus\{v\} : \lVert \mathbf{x}_u -
  \mathbf{x}_v \rVert \le r\}`$ — radius-$`r`$ neighborhood, computed
  with `RANN::nn2(..., searchtype = "radius", radius = r, k = k_\max)`.
- $`\mathcal{V}_c = \{v : c(v) = c\}`$ — the set of cells in cluster
  $`c`$.

## `generate_sim()` — synthetic spatial layouts

[`generate_sim()`](https://juninamo.github.io/cohalu/reference/generate_sim.md)
produces a single sample under one of three generative schemes. The user
supplies `close_ratio` $`\rho`$, `distance_param` $`\Delta`$, `n_cells`
$`N`$, `n_types` $`K`$, and `max_loc` $`L`$ (the side length of the
bounding box).

In every case the *non-special* cells are drawn uniformly from
$`\mathcal{U}\bigl([0, L]^2\bigr)`$ and the cluster labels are sampled
i.i.d. from a uniform distribution over $`\{c_1, \dots, c_K\}`$. What
differs between `test_type`s is how `cell_type_1` and `cell_type_2` are
positioned.

### `test_type = "distribute"`

Start from $`N`$ i.i.d. samples in $`[0, L]^2`$ with i.i.d. cluster
labels. Let $`I_k = \{v : c(v) = c_k\}`$ denote the index set of cluster
$`k`$. A fraction $`\rho`$ of `cell_type_2` cells, of size

``` math
n_\text{pair} = \mathrm{round}\bigl(\rho \cdot |I_2|\bigr),
```

is randomly chosen and *moved* to be near a uniformly random
`cell_type_1` cell:

``` math
\theta_m \stackrel{\text{iid}}{\sim} \mathcal{U}(0, 2\pi),
\quad
\boldsymbol{\eta}_m \stackrel{\text{iid}}{\sim} \mathcal{N}\!\left(\mathbf{0},\, (\Delta/5)^2 I_2\right),
```

``` math
\mathbf{x}_{v_m^{(2)}} \;\leftarrow\;
  \mathbf{x}_{v_m^{(1)}} \;+\; \Delta \binom{\cos \theta_m}{\sin \theta_m} \;+\; \boldsymbol{\eta}_m,
\quad m = 1, \dots, n_\text{pair},
```

where $`(v_m^{(1)}, v_m^{(2)})`$ are randomly paired indices from
$`I_1`$ and $`I_2`$. So $`\rho`$ controls *how many* type-2 cells are
seeded next to a type-1 cell, and $`\Delta`$ controls *how close* the
pair is.

### `test_type = "circle"`

[`generate_sim()`](https://juninamo.github.io/cohalu/reference/generate_sim.md)
sets the disk centre at $`\mathbf{c} = (L/2, L/2)`$ and the radius at
$`R = \lceil L/6 \rceil`$. Of $`n_1 = \lceil N/K \rceil`$`cell_type_1`
cells, a fraction $`\rho`$ is placed inside / on the disk:

``` math
\begin{aligned}
n_1^{\text{close}} &= \mathrm{round}(\rho \, n_1), \\
\text{outer ring:} \quad
\mathbf{x}_v &= \mathbf{c} + R\,(\cos\theta, \sin\theta),
\quad \theta \in \text{uniform grid on } [0, 2\pi], \\
\text{inner disk:} \quad
\mathbf{x}_v &= \mathbf{c} + \rho_v (\cos\theta_v, \sin\theta_v),
\quad \rho_v = \sqrt{U_v}\,R, \;
U_v \sim \mathcal{U}(0, 1), \;
\theta_v \sim \mathcal{U}(0, 2\pi),
\end{aligned}
```

with the remaining $`n_1 - n_1^\text{close}`$ cells dropped uniformly in
$`[0, L]^2`$. Of $`n_2 = \lceil N/K \rceil`$`cell_type_2` cells,
$`n_2^\text{close} = \mathrm{round}(\rho \, n_2)`$ are placed on the
outer ring $`R + \Delta`$:

``` math
\mathbf{x}_v = \mathbf{c} + (R + \Delta)\,(\cos\theta, \sin\theta) + \boldsymbol{\eta},
\quad \boldsymbol{\eta} \sim \mathcal{N}(\mathbf{0}, 5^2 I_2).
```

The remainder of `cell_type_2` is uniform, as are the other $`K - 2`$
clusters. The result is the familiar “type-1 disk + type-2 ring” layout:
the two clusters co-occur at the scale of the whole disk, but a
small-$`k`$ local neighborhood graph mostly sees same-cluster neighbors.

### `test_type = "line"`

Two thin horizontal layers separated vertically by $`\Delta`$. Layer 0
contains `cell_type_1` along an evenly spaced $`x`$-axis grid with
Gaussian jitter; layer 1 contains `cell_type_2` at $`y`$-position
shifted by $`\Delta`$. A fraction $`\rho`$ of layer-1 cells are
*snapped* to the $`(x, y + \Delta)`$ position of a randomly chosen
layer-0 cell, so that $`\rho`$ again controls how many `cell_type_2`
cells are tightly paired with a `cell_type_1` cell.

## `nhood_enrichment()` — neighbourhood enrichment

#### Step 1 — build the kNN adjacency

Construct the exact kNN graph with a kd-tree
(`RANN::nn2(coords, k = k)`; the cell counts as its own neighbour, as in
[`Seurat::FindNeighbors`](https://satijalab.org/seurat/reference/FindNeighbors.html)).
Let $`A \in \{0, 1\}^{N \times N}`$ be the resulting adjacency, or, with
`connectivity_key = "snn"`, Seurat’s shared-nearest-neighbour matrix.

#### Step 2 — degree weighting (optional, off by default)

By default (`transformation = FALSE`) every kNN link counts once, as in
squidpy. When `transformation = TRUE`, each link from cell $`u`$ is
divided by the number of cells that chose $`u`$ as a neighbour (plus
one), which mildly down-weights “hub” cells chosen by many others. With
a kNN graph every cell sends $`k`$ links at any density and $`d_u`$
averages $`k`$ in both dense and sparse tissue; it differs mainly at
density transitions (e.g. cells just outside a dense cluster are chosen
less often). Density itself is handled by the label shuffles, which keep
every position. In simulations the weighting kept the calibration but
lowered power (planted pair, 1,500 cells: 14% vs 28%), so it is off by
default:

``` math
\tilde A_{uv} \;=\; \frac{A_{uv}}{d_u + 1},
\qquad
d_u = \sum_{w} A_{wu}.
```

Otherwise $`\tilde A = A`$ and the function counts edges as binary.

#### Step 3 — cluster-pair count matrix

For every ordered cluster pair $`(c_i, c_j)`$, sum entries of
$`\tilde A`$ whose row cell belongs to $`c_i`$ and column cell to
$`c_j`$:

``` math
C_{ij} \;=\; \sum_{u \in \mathcal{V}_{c_i}} \sum_{v \in \mathcal{V}_{c_j}} \tilde A_{uv}.
```

This is the observed value, implemented in
[`compute_count()`](https://juninamo.github.io/cohalu/reference/compute_count.md)
as $`M^\top \tilde A M`$ with one-hot cluster indicators $`M`$.

#### Step 4 — permutation null

[`permute_clusters()`](https://juninamo.github.io/cohalu/reference/permute_clusters.md)
draws **one** uniform random permutation $`\pi^{(b)}`$ of the cluster
labels and applies it to rows and columns alike, so every cell keeps a
single (random) label, and recomputes the count:

``` math
C^{(b)}_{ij} \;=\; \sum_{u \in \mathcal{V}^{(\pi^{(b)})}_{c_i}}
                   \sum_{v \in \mathcal{V}^{(\pi^{(b)})}_{c_j}}
                   \tilde A_{uv},
\quad b = 1, \dots, B.
```

Positions, the graph and the number of cells of each type are fixed;
only which type sits where is random. (Versions $`\le`$ 0.99.1 shuffled
rows and columns independently, which inflated same-type z-scores.)

#### Step 5 — effect size: centred log2 O/E

With $`E_{ij} = \frac{1}{B}\sum_b C^{(b)}_{ij}`$ and a pseudocount $`c`$
of one mean edge weight,

``` math
\ell_{ij} = \log_2\frac{C_{ij} + c}{E_{ij} + c},
\qquad
\boxed{\;\mathrm{log2\_oe}_{ij} \;=\; \ell_{ij} - \frac{1}{B}\sum_{b}\log_2\frac{C^{(b)}_{ij} + c}{E_{ij} + c}\;}
```

The log of a ratio of small counts is biased below 0 (Jensen), so rare
cell types get $`\ell_{ij} < 0`$ without any interaction. Subtracting
the mean of the same statistic over the shuffles centres it at 0 under
the null for any number of cells. $`\ell_{ij}`$ is returned as
`log2_oe_raw`. `log2_oe` does not grow with the number of cells and is
the value to compare between samples.

#### Step 6 — within-sample test: `pvalue` and `padj`

For each unordered pair the two directions are summed,
$`S_{ij} = C_{ij} + C_{ji}`$ ($`S_{ii} = C_{ii}`$), and standardised
together with the shuffles, using the mean $`\mu_{ij}`$ and SD
$`s_{ij}`$ of the $`B + 1`$ values
$`\{S_{ij}, S^{(1)}_{ij}, \dots, S^{(B)}_{ij}\}`$:

``` math
z_{ij} = \frac{|S_{ij} - \mu_{ij}|}{s_{ij}},
\qquad
\mathrm{pvalue}_{ij} = 2\,\Phi(-z_{ij}).
```

`pvalue` is accurate for a single pre-specified pair. For all
$`K(K+1)/2`$ pairs at once, `padj` is the Westfall–Young single-step
max-T adjusted p-value

``` math
\boxed{\;\mathrm{padj}_{ij} \;=\; \frac{1 + \#\{b : \max_{kl} z^{(b)}_{kl} \ge z_{ij}\}}{B + 1}\;}
```

which controls the family-wise error rate and adapts to the skewed,
dependent null of rare pairs (Benjamini–Hochberg on `pvalue`, returned
as `padj_bh`, can exceed its level when there are many rare cell types).
The smallest possible `padj` is $`1/(B+1)`$. The unsymmetrised
$`z_{ij} = (C_{ij} - E_{ij}) / \mathrm{SD}(C^{(b)}_{ij})`$ is still
returned as `zscore`; it grows with the number of cells.

#### Step 7 — directional statistics: `contact` and `dominance`

Every $`i`$–$`j`$ contact is also a $`j`$–$`i`$ contact, so the
pair-level O/E is (nearly) symmetric:
$`\mathrm{O/E}_{i\to j} = \frac{\text{share of } j
\text{ among the neighbours of } i}{p_j} \approx \frac{\text{share of } i
\text{ among the neighbours of } j}{p_i}`$. It answers “is there more
$`i`$–$`j`$ contact than chance?”, not “from whose point of view”. The
heatmap of `plot_nhood_heatmap(value = "log2_oe")` shows the average of
the two directions, once per pair (lower triangle). Two directional
questions (row = centre cell type $`i`$, column = neighbour type $`j`$)
are answered on the unweighted $`k`$-NN graph without the cell itself,
with $`n_{u,j}`$ the number of type-$`j`$ neighbours of cell $`u`$ and
$`d_u`$ its number of neighbours:

``` math
\mathrm{contact}_{ij} = \frac{1}{n_i}\sum_{u \in i} 1[n_{u,j} \ge 1],
\qquad
\mathrm{dominance}_{ij} = \frac{1}{n_i}\sum_{u \in i} 1[n_{u,j} \ge d_u / 2].
```

`contact` asks how much of population $`i`$ touches $`j`$; `dominance`
asks whether the neighbourhoods of $`i`$ are dominated by $`j`$. Both
are compared with the same label shuffles (`*_expected`, centred
`*_log2_oe`), tested per ordered pair (`*_pvalue`) and adjusted by max-T
over all $`K^2`$ ordered pairs (`*_padj`). For a rare type $`A`$ that
always sits inside clusters of an abundant type $`B`$, `dominance[A, B]`
is enriched and `dominance[B, A]` is not, while `contact[B, A]` is
enriched (more $`B`$ cells than expected touch an $`A`$) and
`contact[A, B]` saturates near 1. Read the two together.

In code:

``` r

res <- nhood_enrichment(df, cluster_key = "cell_type", neighbors.k = 10,
                        n_perms = 200, seed = 1, n_jobs = 1)
res$log2_oe   # effect size (centred)
res$padj      # max-T adjusted, all pairs
plot_nhood_heatmap(res)                               # pair level (symmetric)
plot_nhood_heatmap(res, value = "dominance_log2_oe")  # directional
```

## `calc_co_occurrence_for_radius()` and `compute_co_occurrence_ratio()`

Where
[`nhood_enrichment()`](https://juninamo.github.io/cohalu/reference/nhood_enrichment.md)
works on a fixed-$`k`$ graph, these two functions work on a
fixed-*radius* neighborhood and report a conditional-probability ratio
instead of a z-score.

#### Step 1 — radius-based neighbor counts

For each cell $`v`$, find its radius neighborhood $`N_r(v)`$ via
[`RANN::nn2`](https://jefferislab.github.io/RANN/reference/nn2.html)
(with cap $`k`$). Tally pairs by cluster:

``` math
C_{ij} \;=\; \sum_{v \in \mathcal{V}_{c_i}}\; \bigl|N_r(v) \cap \mathcal{V}_{c_j}\bigr|.
```

Note that the same physical pair $`(v, u)`$ contributes once to
$`C_{c(v), c(u)}`$ and once to $`C_{c(u), c(v)}`$, so $`C`$ is symmetric
when accumulated over all cells. The matrix is stored under
`co_occur_count`.

#### Step 2 — enrichment ratio

[`compute_co_occurrence_ratio()`](https://juninamo.github.io/cohalu/reference/compute_co_occurrence_ratio.md)
then converts $`C`$ to a ratio of empirical probabilities (Lifshitz /
Giotto convention):

``` math
P(j \mid i) \;=\; \frac{C_{ij}}{\sum_k C_{ik}},
\qquad
P(j) \;=\; \frac{\sum_k C_{kj}}{\sum_{k, l} C_{kl}},
```

``` math
\boxed{\;
r_{ij} \;=\; \frac{P(j \mid i)}{P(j)}\;.
}
```

$`r_{ij} > 1`$ means a cell in cluster $`c_i`$ is more likely to have a
cluster-$`c_j`$ neighbor than chance predicts from the overall pair
composition. Unlike the z-score in
§[`nhood_enrichment()`](https://juninamo.github.io/cohalu/reference/nhood_enrichment.md),
this ratio is deterministic — no permutation is run.

## `cooccur_local()` — per-cell local co-occurrence

Where the previous two functions produce a cluster-pair matrix, this one
produces a *per-cell score* for a single pair $`(c_x, c_y)`$ chosen by
the caller. It is meant to highlight individual cells that sit at the
interface of two cluster populations.

#### Step 1 — local presence indicator

Take the radius-$`r`$ neighborhood $`N_r(v)`$ (with neighbor cap $`k`$).
The raw local co-occurrence score is

``` math
s_v^{(0)} \;=\; \mathbb{1}\!\left\{
  N_r(v) \cap \mathcal{V}_{c_x} \ne \emptyset
  \;\wedge\;
  N_r(v) \cap \mathcal{V}_{c_y} \ne \emptyset
\right\}.
```

So $`s_v^{(0)} = 1`$ exactly when at least one cell of each target
cluster appears within distance $`r`$ of $`v`$ (and $`v`$ itself is
excluded). The score is symmetric in $`(c_x, c_y)`$.

#### Step 2 — single graph-diffusion step

The raw indicator is sparse.
[`cooccur_local()`](https://juninamo.github.io/cohalu/reference/cooccur_local.md)
smooths it once over the kNN graph $`A`$ (same construction as in
[`nhood_enrichment()`](https://juninamo.github.io/cohalu/reference/nhood_enrichment.md)).
With $`D = \mathrm{diag}(d_v + 1)`$:

``` math
\boxed{\;
s_v^{(1)} \;=\; \frac{s_v^{(0)}}{d_v + 1} \;+\; \sum_{u} A_{vu} \, \frac{s_u^{(0)}}{d_u + 1},
\quad\text{i.e.}\quad
\mathbf{s}^{(1)} \;=\; (A + I)\, D^{-1}\, \mathbf{s}^{(0)}.
\;}
```

Concretely, every cell shares an equal portion of its own indicator with
each of its kNN neighbors plus itself. The output column is named
`cooccur_local_<cluster_x>_<cluster_y>`.

With `maxnsteps > 1` the step is iterated, each step starting from the
previous result:
$`\mathbf{s}^{(t)} = (A + I) D^{-1} \mathbf{s}^{(t-1)}`$. The total
$`\sum_v s_v`$ is conserved. Iteration stops early when the kurtosis of
the scores decreases by less than 3 between consecutive steps (checked
after step 3). Because $`\mathbf{s}^{(0)}`$ is a 0/1 indicator, its
kurtosis usually *increases* under smoothing, so in practice diffusion
stops after 4 steps whenever `maxnsteps >= 4`. `maxnsteps = 0` returns
the raw indicator. (Versions \<= 0.99.1 recomputed every step from
$`\mathbf{s}^{(0)}`$, so any `maxnsteps >= 1` gave the single step
above.)

## `cooccur_local_oe()` — calibrated local co-localization

For cell $`v`$ with $`N_v`$ other cells within radius $`r`$, of which
$`n_A(v)`$ are of type A and $`n_B(v)`$ of type B, the number of A–B
pairs is $`s_v = n_A(v)\, n_B(v)`$ ($`n_A(n_A - 1)`$ for A = B). Its
exact expectation when labels are shuffled over fixed positions is
multivariate hypergeometric:

``` math
E_v = N_v (N_v - 1)\,\frac{K_A K_B}{M (M - 1)},
```

with $`M`$ the other cells and $`K_A, K_B`$ the A and B cells among
them. Both are smoothed with Gaussian weights
$`w_{vu} = \exp(-d_{vu}^2 / 2\sigma^2)`$, $`\sigma = h/2`$ (`bandwidth`
$`= h`$, truncated at $`h`$, self weight 1):

``` math
\mathrm{local\_log2\_oe}_v = \log_2 \frac{\sum_u w_{vu} s_u + c}{\sum_u w_{vu} E_u + c}.
```

The section value is $`\log_2(\sum_v s_v / \sum_v E_v)`$. Hotspot
p-values use $`B`$ label shuffles:
$`p_v = (1 + \#\{b : \tilde s^{(b)}_v \ge \tilde s_v\}) / (B + 1)`$ with
$`\tilde s_v = \sum_u w_{vu} s_u`$, then BH over cells. One neighbour
search, $`O(nk)`$; each shuffle is one more pass.

## Group comparison — `compare_groups()`

Per-image scores $`y_{pm}`$ (image $`m`$ of patient $`p`$,
e.g. `log2_oe`) come from
[`nhood_enrichment_per_sample()`](https://juninamo.github.io/cohalu/reference/nhood_enrichment_per_sample.md)
and related helpers. For every cell-type pair, with the patient as the
unit:

- `"wilcox"`: patient means $`\bar y_p = \frac{1}{m_p}\sum_m y_{pm}`$,
  Wilcoxon rank-sum (exact for small samples).
- `"lmm"`:
  $`y_{pm} = \beta_0 + \beta_1\,\mathrm{case}_p + \gamma^\top c_p + u_p + \varepsilon_{pm}`$,
  $`u_p \sim N(0, \tau^2)`$, $`\varepsilon_{pm} \sim N(0, \sigma^2)`$;
  test of $`\beta_1`$ with Satterthwaite degrees of freedom (lmerTest).
  The group ($`\beta_0, \beta_1`$) and covariates ($`\gamma`$) are fixed
  effects; the patient is a random effect (random intercept $`u_p`$,
  i.e. `value ~ group + covariates + (1 | patient)`), so the images of a
  patient are treated as correlated observations. Patient cannot be a
  fixed effect: each patient belongs to one group, so patient dummies
  would be confounded with \$eta_1\$.
- `"perm"`:
  $`T = \overline{\bar y}_{\mathrm{case}} - \overline{\bar y}_{\mathrm{control}}`$
  compared with patient-level label shuffles (within-patient sign flips
  for paired data).
- `"signrank"`:
  $`d_p = \bar y_{p,\mathrm{post}} - \bar y_{p,\mathrm{pre}}`$, Wilcoxon
  signed-rank.

Benjamini–Hochberg over pairs.
[`summarize_by_patient()`](https://juninamo.github.io/cohalu/reference/summarize_by_patient.md)
averages images within patient, pair and distance `r`.

## Continuous variables — `associate_continuous()`

For a patient-level variable $`x_p`$ (CRP, disease activity, age):
Spearman correlation of $`\bar y_p`$ with $`x_p`$ (default; exact p for
small cohorts), linear regression
$`\bar y_p = \beta_0 + \beta_1 x_p + \gamma^\top c_p + \varepsilon_p`$,
the linear mixed model
$`y_{pm} = \beta_0 + \beta_1 x_p + \gamma^\top c_p + u_p + \varepsilon_{pm}`$
(patient as a random intercept $`u_p`$:
`value ~ x + covariates + (1 | patient)`), or a permutation test that
shuffles $`x_p`$ between patients. BH over pairs.

## Segmentation-free analysis

#### Pair correlation of transcripts — `pcf_cross()`, `pcf_matrix()`

Transcripts are counted in square bins ($`Y_A(u)`$ = transcripts of gene
set A in bin $`u`$, inside the tissue mask $`W`$). With the FFT
cross-correlation $`P_{AB}(h) = \sum_u Y_A(u)\, Y_B(u + h)`$ and the
mask overlap $`O(h) = \sum_u W(u) W(u + h)`$,

``` math
g_{AB}(r) = \frac{\sum_{|h| \approx r} P_{AB}(h)}{\lambda_A \lambda_B \sum_{|h| \approx r} O(h)},
\qquad
g^{\mathrm{rel}}_{AB}(r) = \frac{g_{AB}(r)}{g_{\mathrm{all}}(r)},
```

where $`\lambda`$ are mean intensities and $`g_{\mathrm{all}}`$ is the
pair correlation of all transcripts. $`g^{\mathrm{rel}}`$ equals
observed / expected under gene-label shuffling over fixed transcript
positions, the transcript counterpart of log2 O/E; cellularity cancels.

#### Random-feature log-Gaussian Cox process — `fit_spatial_rff()`

``` math
y_{bj} \sim \mathrm{NB}(\mu_{bj}),\quad
\log \mu_{bj} = \log a + \alpha_j + \sigma_0 f_0(u_b) + \sum_{k=1}^K L_{jk} f_k(u_b),
\quad f_k(u) = \phi_k(u)^\top \gamma_k,
```

with random Fourier features
$`\phi(u) = \sqrt{2/M}\,[\cos(\omega_m^\top u / \ell + b_m)]_m`$,
$`\omega_m \sim N(0, I)`$ (Gundersen, Zhang & Engelhardt, AISTATS 2021).
The coordinates $`u_b`$ are observed, so only the fields and loadings
are estimated (MAP, L-BFGS-B, analytic gradients). Model pair
correlations follow from
$`g_{AB}(r) = \exp \mathrm{Cov}[\log\lambda_A(u), \log\lambda_B(u + r)]`$;
removing $`f_0`$ gives the composition version. Because each $`f_k`$ is
a function of the coordinates,
[`rff_fields()`](https://juninamo.github.io/cohalu/reference/rff_fields.md)
evaluates the fitted (unit-variance) fields at any point - cell
centroids or transcripts - and averages them per cell when a cell id is
given.

**Residual factors.** To look for spatial programs that known structure
does not explain, a per-bin, per-gene log offset $`o_{bj}`$ is added,
$`\log \mu_{bj} = \log a + o_{bj} + \alpha_j + \sigma_0 f_0(u_b) + \sum_k L_{jk} f_k(u_b)`$.
[`rff_offset()`](https://juninamo.github.io/cohalu/reference/rff_offset.md)
builds $`o_{bj}`$ from covariates of the bins (cell-type composition,
domains, or an embedding such as PCA, Harmony or SCIGMA) by a
ridge-penalised Poisson regression per gene. The penalty
$`\lambda \sum_k \lVert L_{\cdot k} \rVert_2`$ (`ard`) shrinks whole
factors that are not needed.
[`rff_factor_test()`](https://juninamo.github.io/cohalu/reference/rff_factor_test.md)
simulates counts from the fitted model without factors, refits, and
compares each factor’s strength $`\lVert L_{\cdot k} \rVert_2`$ with the
largest null factor (parametric bootstrap, family-wise over factors).

#### Gene-level modules — `colocalization_gene_matrix()`, `colocalization_modules()`, `module_enrichment()`

With a disc kernel $`K_r`$ (one FFT convolution per gene),

``` math
P_{ab} = \sum_u Y_a(u)\, (K_r * Y_b)(u) \;(- n_a \text{ if } a = b),
\qquad
E_{ab} = \frac{n_a n_b}{N (N - 1)}\, P_{\mathrm{all}},
```

and $`\log_2((P_{ab} + c) / (E_{ab} + c))`$ is 0 without
co-localization. Genes are clustered (average linkage, distance
$`\max M - M`$) into modules;
[`module_enrichment()`](https://juninamo.github.io/cohalu/reference/module_enrichment.md)
tests any gene sets per module with the hypergeometric test and BH.

## `search_interaction_spot()` — connected-component spots

This function turns a radius-neighborhood graph into discrete
“interaction spots” — connected groups of cells that are likely to be
acting on each other in a tissue.

#### Step 1 — radius-neighbor graph

Restrict to a user-provided cell ID list $`\mathcal{V}^\star \subseteq
\mathcal{V}`$. For each $`v \in \mathcal{V}^\star`$ compute $`N_r(v)`$
via [`RANN::nn2`](https://jefferislab.github.io/RANN/reference/nn2.html)
(radius search, neighbor cap $`k`$). Build an undirected graph
$`G^\star = (\mathcal{V}^\star, E^\star)`$ with

``` math
E^\star \;=\; \bigl\{\{v, u\} : v \in \mathcal{V}^\star,\; u \in N_r(v) \cap \mathcal{V}^\star\bigr\}.
```

#### Step 2 — weak connected components

Run `igraph::components(G^\star, mode = "weak")`. Each component
$`\Gamma_q \subseteq \mathcal{V}^\star`$ is a candidate spot, and every
cell $`v`$ receives an integer label $`q(v)`$.

#### Step 3 — keep target-cluster, size-$`n_\min`$ spots

Define the bounding box of component $`q`$ from its cells’ coordinates,

``` math
\bigl[x_q^{\min}, x_q^{\max}\bigr] \times \bigl[y_q^{\min}, y_q^{\max}\bigr],
\quad
n_q = |\Gamma_q|,
```

and keep cell $`v`$ in the output iff

``` math
c(v) \in \mathcal{T} \quad\text{and}\quad n_{q(v)} > n_\min,
```

where $`\mathcal{T}`$ is the `target_cluster` argument. The function
returns a tidy data.frame of surviving cells, one row per cell, with
columns `cluster_id`, `x_min`, `x_max`, `y_min`, `y_max`, and
`n_all_cells`.

## Summary

| Function | Unit | Local geometry | Statistic | Test |
|----|----|----|----|----|
| [`nhood_enrichment()`](https://juninamo.github.io/cohalu/reference/nhood_enrichment.md) | cell-type pair (and ordered pair) | $`k`$-NN graph | centred log2 O/E; directional `contact`, `dominance` | shuffles: `pvalue`, max-T `padj` |
| [`cooccur_local_oe()`](https://juninamo.github.io/cohalu/reference/cooccur_local_oe.md) | cell, for a pair | radius $`r`$ + Gaussian smoothing | local log2 O/E | shuffles per cell, BH |
| [`calc_co_occurrence_for_radius()`](https://juninamo.github.io/cohalu/reference/calc_co_occurrence_for_radius.md) / [`compute_co_occurrence_ratio()`](https://juninamo.github.io/cohalu/reference/compute_co_occurrence_ratio.md) | cell-type pair | radius $`r`$ | $`P(j \mid i) / P(j)`$ | — |
| [`cooccur_local()`](https://juninamo.github.io/cohalu/reference/cooccur_local.md) | cell, for a pair | radius $`r`$ + graph diffusion | 0/1 indicator (sCLS) | — |
| [`search_interaction_spot()`](https://juninamo.github.io/cohalu/reference/search_interaction_spot.md) | cell | radius $`r`$ + components | component size | — |
| [`compare_groups()`](https://juninamo.github.io/cohalu/reference/compare_groups.md) | patient | per-image scores | group difference | Wilcoxon, LMM, permutation, signed-rank |
| [`associate_continuous()`](https://juninamo.github.io/cohalu/reference/associate_continuous.md) | patient | per-image scores | slope / Spearman $`\rho`$ | Spearman, lm, LMM, permutation |
| [`pcf_cross()`](https://juninamo.github.io/cohalu/reference/pcf_cross.md) / [`pcf_matrix()`](https://juninamo.github.io/cohalu/reference/pcf_matrix.md) | gene-set pair × $`r`$ | bins + FFT | (relative) $`g(r)`$ | via [`compare_groups()`](https://juninamo.github.io/cohalu/reference/compare_groups.md) |
| [`fit_spatial_rff()`](https://juninamo.github.io/cohalu/reference/fit_spatial_rff.md) | genes | random-feature LGCP | model $`g(r)`$ | — |
| [`colocalization_gene_matrix()`](https://juninamo.github.io/cohalu/reference/colocalization_gene_matrix.md) | gene pair | disc of radius $`r`$ | log2 O/E | modules + hypergeometric |

The multi-sample functions run the single-sample statistics per image
and test the per-image values with the patient as the unit; see the
[case-control
tutorial](https://juninamo.github.io/cohalu/articles/case_control_tutorial.md).

## Session information

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] cohalu_0.99.3
#> 
#> loaded via a namespace (and not attached):
#>   [1] deldir_2.0-4           pbapply_1.7-5          gridExtra_2.3.1       
#>   [4] rlang_1.3.0            magrittr_2.0.5         RcppAnnoy_0.0.23      
#>   [7] otel_0.2.0             spatstat.geom_3.8-3    matrixStats_1.5.0     
#>  [10] ggridges_0.5.7         compiler_4.6.1         png_0.1-9             
#>  [13] systemfonts_1.3.2      vctrs_0.7.3            reshape2_1.4.5        
#>  [16] stringr_1.6.0          pkgconfig_2.0.3        fastmap_1.2.0         
#>  [19] promises_1.5.0         rmarkdown_2.32         ragg_1.5.2            
#>  [22] purrr_1.2.2            xfun_0.61              cachem_1.1.0          
#>  [25] jsonlite_2.0.0         goftest_1.2-3          later_1.4.8           
#>  [28] spatstat.utils_3.2-5   irlba_2.3.7            parallel_4.6.1        
#>  [31] cluster_2.1.8.2        R6_2.6.1               ica_1.0-3             
#>  [34] spatstat.data_3.1-9    bslib_0.12.0           stringi_1.8.9         
#>  [37] RColorBrewer_1.1-3     reticulate_1.47.0      spatstat.univar_3.2-0 
#>  [40] parallelly_1.48.0      lmtest_0.9-40          jquerylib_0.1.4       
#>  [43] scattermore_1.2        Rcpp_1.1.2             knitr_1.52            
#>  [46] tensor_1.5.1           future.apply_1.20.2    zoo_1.9-1             
#>  [49] sctransform_0.4.3      httpuv_1.6.17          Matrix_1.7-5          
#>  [52] splines_4.6.1          igraph_2.3.3           tidyselect_1.2.1      
#>  [55] abind_1.4-8            yaml_2.3.12            spatstat.random_3.5-2 
#>  [58] spatstat.explore_3.8-3 codetools_0.2-20       miniUI_0.1.2          
#>  [61] listenv_1.0.0          lattice_0.22-9         tibble_3.3.1          
#>  [64] plyr_1.8.9             shiny_1.14.0           S7_0.2.2              
#>  [67] ROCR_1.0-12            evaluate_1.0.5         Rtsne_0.17            
#>  [70] future_1.76.0          fastDummies_1.7.6      desc_1.4.3            
#>  [73] survival_3.8-6         polyclip_1.10-7        fitdistrplus_1.2-6    
#>  [76] pillar_1.11.1          Seurat_5.5.1           KernSmooth_2.23-26    
#>  [79] plotly_4.12.1          generics_0.1.4         RcppHNSW_0.7.0        
#>  [82] sp_2.2-3               ggplot2_4.0.3          scales_1.4.0          
#>  [85] globals_0.19.1         xtable_1.8-8           glue_1.8.1            
#>  [88] tools_4.6.1            data.table_1.18.6.1    RSpectra_0.16-2       
#>  [91] RANN_2.6.3             fs_2.1.0               dotCall64_1.2         
#>  [94] cowplot_1.2.0          grid_4.6.1             tidyr_1.3.2           
#>  [97] nlme_3.1-169           patchwork_1.3.2        cli_3.6.6             
#> [100] spatstat.sparse_3.2-0  textshaping_1.0.5      spam_2.11-4           
#> [103] viridisLite_0.4.3      dplyr_1.2.1            uwot_0.2.5            
#> [106] gtable_0.3.6           sass_0.4.10            digest_0.6.39         
#> [109] progressr_1.0.0        ggrepel_0.9.8          htmlwidgets_1.6.4     
#> [112] SeuratObject_5.4.0     farver_2.1.2           htmltools_0.5.9       
#> [115] pkgdown_2.2.1          lifecycle_1.0.5        httr_1.4.9            
#> [118] mime_0.13              MASS_7.3-65
```
