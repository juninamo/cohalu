# Reach (decay length) implied by a profiled length scale

\*\*Experimental.\*\* The length scale \`ell\` of
\[rff_lengthscale_profile()\] describes a Gaussian (RBF) correlation. If
a program is a response to point sources (producer cells) that decays as
\\e^{-d/\lambda}\\, the response map is shot noise whose correlation is
the self-convolution of that kernel, which is matched at half
correlation by an RBF length scale of about \\1.72\lambda\\ (more for
\`lambda\` near the bin size, where bin averaging adds to the width).
\`rff_reach()\` inverts this relation. The result is a reach under that
model only, and it is only identified when the length scale is well
inside the window (see \`rff_lengthscale_profile()\`'s \`identifiable\`
column).

## Usage

``` r
rff_reach(ell, bin_size, kernel = c("exponential", "gaussian"))
```

## Arguments

- ell:

  Length scale(s) (coordinate units).

- bin_size:

  Bin size used for the fit.

- kernel:

  \`"exponential"\` (decay \\e^{-d/\lambda}\\) or \`"gaussian"\` (decay
  \\e^{-d^2/2\lambda^2}\\: its shot noise has RBF length scale
  \\\sqrt{2}\lambda\\, bin averaging ignored).

## Value

Numeric vector of reaches (\`NA\` where \`ell\` is too small to be
produced by any reach at this bin size).

## Examples

``` r
rff_reach(c(20, 70, 140, 550), bin_size = 8)
#> [1]  11.44355  40.60708  81.29063 319.41827
```
