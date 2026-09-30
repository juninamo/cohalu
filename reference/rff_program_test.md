# Gene-program test on smoothed residuals, with the fitted RFLVM as null model

## Usage

``` r
rff_program_test(
  fit,
  binned,
  bandwidth = NULL,
  n_components = 6,
  n_boot = 99,
  sequential = TRUE,
  alpha = 0.05,
  seed = 1,
  null = c("parametric", "shift"),
  highpass = NULL,
  crossfit = FALSE,
  crossfit_block = NULL,
  control = NULL,
  min_effect = NULL,
  n_boot_ci = 200,
  robust = FALSE,
  control_stat = c("direction", "rank")
)
```

## Arguments

- fit:

  Output of \[fit_spatial_rff()\] with a matrix or \`"area"\` offset.

- binned:

  The \`binned_transcripts\` object used for the fit.

- bandwidth:

  Standard deviation of the Gaussian smoothing (coordinate units);
  default the bin size.

- n_components:

  Number of components reported (and at most tested).

- n_boot:

  Number of simulated null data sets per tested component.

- sequential:

  Sequential testing (see Details); \`FALSE\` tests every component
  against the largest variance share of the null.

- alpha:

  Level at which sequential testing stops.

- seed:

  Random seed.

- null:

  How the null distribution is generated. \`"parametric"\` (default):
  counts simulated from the fitted null model (independent noise around
  a smooth mean), as described above. It assumes that residual variation
  not shared by genes is spatially independent; when genes have their
  own spatially autocorrelated residual structure (usual in real tissue)
  it calls components in almost every tissue (100 simulated tissues with
  gene-specific residual fields) - significance then means only
  "structure the offset misses". \`"shift"\`: gene-shift surrogates that
  keep each gene's own spatial autocorrelation. Every gene's counts are
  moved together with its null mean by an independent random rigid
  transformation of the bin grid (flips, transposition for square grids
  and a shift on the mirror-extended grid, so that no seams are
  created), which destroys only the alignment between genes. Residuals
  use the Poisson variance, the shared multiplicative direction is
  removed and each bin is normalised over genes; component \`k\` is
  compared with the \`k\`-th variance share of the surrogates (parallel
  analysis) and sequential p-values are cumulative maxima. Significance
  then means a program shared by several genes. In simulations (300 x
  300 um, 30 genes) it called 0/40 tissues without any residual
  structure, 1/40, 1/40 and 4/40 with independent gene-specific residual
  fields of length scale 10, 20 and 40 um, 0/20 with extra cellularity,
  and 80 tissues with a 2-8-gene program. Gene-specific fields at scales
  larger than the program mask it (amplitude 1 program with 20 um fields
  of sd 0.25: 4/30 called); use \`highpass\` then. A covariate missing
  from the offset that moves several genes together is a shared
  structure and is called. On real Xenium windows, data with every gene
  shifted independently were called in 0/30 replicates for three windows
  but in 3/10 and 5/10 for two (an irregular tissue mask; a window
  dominated by a large lymphoid aggregate), so borderline p-values on
  real data should be read with care.

- highpass:

  Only for \`null = "shift"\`: if given (coordinate units), structure at
  scales above \`highpass\` is removed from the smoothed residuals
  (difference of Gaussians with standard deviations \`bandwidth\` and
  \`highpass\`), so that broad gene-specific residual fields do not mask
  small-scale programs. With \`highpass = 4 \* bandwidth\` (20 um),
  power for an amplitude-1 program with 20 um gene-specific fields rose
  from 4/30 to 26/30 tissues, with 0-3/30 calls without a shared
  program.

- crossfit:

  Only for \`null = "shift"\`. \`FALSE\` (default): the test described
  above. \`TRUE\`: cross-fitted effect sizes. The in-tissue bins are
  split into two halves by a checkerboard of square blocks
  (\`crossfit_block\`); components are found on one half and their share
  of the processed residual variance is measured on the other half (bins
  next to the other half excluded), both ways and averaged. The same
  statistic is computed for gene-shift surrogates: \`p\` compares the
  held-out share with them, \`excess\` = held-out share / surrogate
  mean - 1 is the effect size (0 = no more shared structure than
  gene-shifted data), with a 95 \`excess_upper\`; \`n_boot_ci\`
  replicates), and \`call\` combines \`p\`, the effect-size threshold
  and, if given, the control comparison. Cross-fitting removes the
  selection bias of in-sample shares, but on real tissue shared nuisance
  structure (spill-over, mixing, cellularity) is real and replicates
  across halves: use \`control\` or \`min_effect\` to ask for more than
  that.

- crossfit_block:

  Side of the checkerboard blocks (coordinate units; default 12 bins).

- control:

  Output of \[rff_control_reference()\] for negative-control windows
  (same genes and settings): implies \`crossfit = TRUE\`. Adds
  \`p_control\` (share of control windows whose excess at the same rank
  is at least the target's), \`control_q95\`, and \`p_control_dir\` (the
  direction-specific version: excess of the residual variance along the
  target's component direction in each control window, standardised over
  the controls (\`z_control_dir\`) and compared with the same statistic
  of every control window against the others); \`call\` then also
  requires the control p-value chosen by \`control_stat\` to be \`\<=
  alpha\`.

- min_effect:

  Optional effect-size threshold: a component is called only if
  \`excess_lower\` exceeds it (implies \`crossfit = TRUE\`).

- n_boot_ci:

  Block-bootstrap replicates for the interval of \`excess\`.

- robust:

  Only for \`null = "shift"\` without cross-fitting (the cross-fitted
  test always uses it): floor the Poisson variance of the Pearson
  residuals at 5 give near-zero means for genes of absent cell types, so
  single spill-over transcripts would dominate) and take surrogate
  sources that fall outside the tissue from further independent rigid
  transformations of the same gene (instead of zeros, which created
  shared structure in the surrogates of windows with ragged tissue
  masks).

- control_stat:

  Which comparison with the controls enters \`call\`: \`"direction"\`
  (default; \`p_control_dir\`: is the target's component direction more
  active in the target than in the controls, calibrated by treating
  every control window the same way) or \`"rank"\` (\`p_control\`: is
  the target's k-th component stronger than the controls' k-th
  components?). On 24 non-TLS tiles of a real RA Xenium section both
  called 0-1/24 tiles (leave-one-out) versus 15/24 for the cross-fitted
  test alone; for 6-gene programs planted into the tiles (log-amplitude
  0.75 / 1) the direction-wise comparison found the program in 33 of
  tiles, the rank-wise one in 4

A data frame with one row per component: \`component\`, \`share\`
(variance share), \`p\`, \`factor\` and \`r_factor\` (the fitted factor
whose field correlates best with the component score, and that
correlation) and the top genes by loading; with \`crossfit = TRUE\` also
\`share_heldout\`, \`null_heldout\`, \`excess\`, \`excess_lower\`,
\`excess_upper\`, \`z_heldout\`, \`effect_threshold\`, \`call\` and,
with \`control\`, \`p_control\`, \`control_q95\`, \`z_control_dir\`,
\`p_control_dir\`. Attributes \`scores\` (in-tissue bins x components),
\`loadings\` (genes x components), \`null_shares\` and, for \`null =
"shift"\`, \`null\` (and \`crossfit\`). \*\*Experimental.\*\* A hybrid
of residual PCA and the random-feature model. The fit of
\[fit_spatial_rff()\] supplies the null model: bin area, \`offset\`,
gene intercepts, the cellularity field and the part of every factor that
is shared by all genes (the mean of its loadings), i.e. everything
except gene programs. Pearson residuals of the counts under this null
mean (negative binomial variance with the fitted dispersions) are
Gaussian-smoothed; the per-bin direction of a multiplicative effect
shared by all genes (proportional to \\\sqrt{\mu}\\) is projected out,
genes are standardised, and the variance shares of the principal
components are compared with those of data simulated from the null model
(parametric bootstrap through the same pipeline). With \`sequential =
TRUE\` (default), component \`k\` is compared with the \`k\`-th variance
share of data simulated from the null model plus the gene programs
(loadings minus their mean) of the \`k - 1\` strongest factors of the
fit, so that structure already explained by earlier components is part
of the null; testing stops at the first component with \`p \> alpha\`
(p-values are made non-decreasing).In simulations (300 x 300 um, 30
genes, 2 known domains in the offset; fit with \`basis = "grid"\`, fixed
length scale 8 um, \`ard = 100\`) this test was as powerful as residual
PCA with its own bootstrap (20/20 tissues at program amplitude 0.5, vs
12-13/20 for \[rff_factor_test()\]), called 1/40 tissues whose only
extra structure was cellularity shared by all genes (residual PCA 12/20)
and 3/60 tissues without any program. When to use / limitations The
simulation calibration above does not carry over to real tissue: there
the gene-shift null is valid (gene-shifted data were called in 0/46
windows), but real tissue always carries shared residual structure, and
the plain test called 95 windows and 42 TRUE\` and \`control\` from
\[rff_control_reference()\], built from at least 19 negative-control
regions of the same data type: 0-1 of 24 control tiles were then called
on four data types, at a power cost for weak programs and small windows.
A call means "more residual structure than in control regions", not a
new biological program: check cell-type / sub-lineage residuals and
technical axes with \[rff_report()\], and confirm on held-out data with
\[rff_transfer_test()\].

\[fit_spatial_rff()\], \[rff_factor_test()\], \[rff_offset()\]
