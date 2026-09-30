# Posterior of a program's length scale by MCMC (reference)

\*\*Experimental; small windows.\*\* Reference sampler for checking
\[rff_lengthscale_profile()\] and \[rff_lengthscale_vi()\]: for one
program of \`fit\` with its loading direction fixed and everything else
in the fit as offset, the whitened field weights of the grid-basis RBF
prior are updated by elliptical slice sampling and the log length scale
and log amplitude by random-walk Metropolis (uniform priors on the log
scale), under a Poisson likelihood.

## Usage

``` r
rff_lengthscale_mcmc(
  fit,
  binned,
  factor = 1,
  ls_range = c(2, 640),
  amp_range = c(0.02, 10),
  n_iter = 3000,
  burn = 1000,
  gene_share = 0.95,
  seed = 1
)
```

## Arguments

- fit, binned:

  As in \[rff_lengthscale_profile()\].

- factor:

  The factor (name or index).

- ls_range, amp_range:

  Prior ranges (uniform on the log scale).

- n_iter, burn:

  Iterations and burn-in.

- gene_share:

  Share of the squared loadings carried by the genes used.

- seed:

  Random seed.

## Value

A list with \`samples\` (data frame of \`lengthscale\`, \`amplitude\`
after burn-in), \`summary\` (posterior median and 95 length scale and
amplitude, acceptance rates) and \`field_mean\` (posterior mean field at
the in-tissue bins).
