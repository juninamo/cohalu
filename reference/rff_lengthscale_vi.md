# Length scale of each program by the variational evidence (ELBO)

\*\*Experimental.\*\* Alternative to the held-out criterion of
\[rff_lengthscale_profile()\]: for each program (factor) of \`fit\`,
with its loading direction fixed and everything else in the fit as
offset, the program field gets a Gaussian variational posterior under
the grid-basis RBF prior and a Poisson likelihood (closed-form
expectations; Opper-Archambeau parametrisation, conjugate gradients,
Lanczos log-determinants and Hutchinson traces with FFT products), and
the ELBO - a lower bound on the marginal likelihood, which, unlike the
MAP objective, charges for the flexibility of short length scales - is
maximised over the length scale grid and the field amplitude (prior
variance, golden-section search in log scale). Only the program genes
(share \`gene_share\` of the squared loadings) enter.

## Usage

``` r
rff_lengthscale_vi(
  fit,
  binned,
  factors = NULL,
  ls_grid = 5 * 2^(0:6),
  gene_share = 0.95,
  amp_range = c(0.05, 5),
  n_probe = 8,
  seed = 1,
  n_cores = 1,
  reach_kernel = "exponential",
  nugget = TRUE
)
```

## Arguments

- fit:

  Output of \[fit_spatial_rff()\] (preferably \`basis = "grid"\`).

- binned:

  The \`binned_transcripts\` object used for the fit.

- factors:

  Factors to profile (names or indices of \`fit\$L\`). Default: factors
  with program strength at least 10 loading pattern not dominated by a
  part shared by all genes (uniform share \<= 0.5).

- ls_grid:

  Length scales (coordinate units) to profile.

- gene_share:

  Only genes that carry this share of the program's squared loading norm
  (at least 5 genes) enter the refits; the others hardly inform the
  field.

- amp_range:

  Range of the field amplitude (prior standard deviation times the
  loading norm) searched at each length scale.

- n_probe:

  Hutchinson / Lanczos probe vectors.

- seed:

  Random seed (fold assignment, bootstrap).

- n_cores:

  Cores for the refits (forked; serial on Windows).

- reach_kernel:

  Decay kernel for the implied reach (\[rff_reach()\]):
  \`"exponential"\` (default), \`"gaussian"\`, or \`NULL\` (no reach).

- nugget:

  If \`TRUE\` (default), the program field is a smooth RBF part plus an
  independent bin-level part (relative variance estimated with the
  amplitude by maximising the ELBO). Without it, on real tissue the ELBO
  chose the smallest length scale: responses are carried by single
  cells, so the dominant correlation scale of a response map is the
  cell, not the reach.

## Value

An object of class \`rff_ls_profile\` (as \[rff_lengthscale_profile()\],
\`criterion = "elbo"\`): \`summary\` with \`lengthscale\` (maximum of
the ELBO, parabola in log length scale), \`amplitude\`, \`elbo_gain\`
(ELBO at the best length scale minus the ELBO without the program),
reach and identifiability; \`curves\` with the ELBO per length scale
(\`gain\` = ELBO minus no-program ELBO).

## See also

\[rff_lengthscale_profile()\]
