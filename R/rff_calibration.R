# Length-scale profiling, held-out effect sizes and calibration of the
# residual RFLVM (2026-09 method fixes):
#   * rff_lengthscale_profile(): one length scale per program by held-out
#     likelihood over a grid of fixed length scales;
#   * helpers shared with rff_program_test(crossfit / control) and
#     rff_transfer_test().

# ---- small helpers -----------------------------------------------------------

# Linear predictor of a fit at its in-tissue bins (N x J): bin area (count
# families), offset, gene intercepts, cellularity field and the factors in
# `factors` (column names or indices of fit$L; NULL = none).
.rff_linpred <- function(fit, factors = colnames(fit$L), genes = fit$genes) {
  keep <- fit$in_tissue; N <- sum(keep); J <- length(genes)
  binom <- identical(fit$family, "binomial")
  off <- if (is.null(fit$offset_matrix)) 0 else fit$offset_matrix[keep, genes, drop = FALSE]
  eta <- (if (binom || identical(fit$family, "multinomial")) 0 else 2 * log(fit$bin_size)) + off + matrix(fit$alpha[genes], N, J, byrow = TRUE) +
    (if (fit$has_density) fit$sigma0 * fit$field_grid[keep, "density"] else 0)
  if (length(factors)) {
    if (is.numeric(factors)) factors <- colnames(fit$L)[factors]
    eta <- eta + fit$field_grid[keep, factors, drop = FALSE] %*% t(fit$L[genes, factors, drop = FALSE])
  }
  if (identical(fit$family, "multinomial")) {
    # the bin effect of the full model (all genes, all factors), so that the result is a Poisson log mean
    full <- .rff_linpred_raw(fit)
    eta <- eta + log(pmax(fit$bin_total, 1e-300)) - log(rowSums(exp(full)))
  }
  eta
}
.rff_linpred_raw <- function(fit) {
  f <- fit; f$family <- "poisson"
  .rff_linpred(f, colnames(fit$L)) - 2 * log(fit$bin_size)
}

# Log-likelihood of every entry (N x J) under a family, given the linear
# predictor; `r` dispersions (nb), `n` trials (binomial).
.rff_loglik <- function(Y, eta, family, r = NULL, n = NULL) {
  if (family == "binomial") {
    return(lchoose(n, Y) + Y * eta - n * (pmax(eta, 0) + log1p(exp(-abs(eta)))))
  }
  mu <- exp(pmin(eta, 40))
  if (family == "poisson") return(Y * log(pmax(mu, 1e-300)) - mu - lgamma(Y + 1))
  rr <- matrix(r, nrow(Y), ncol(Y), byrow = TRUE)
  lgamma(Y + rr) - lgamma(rr) - lgamma(Y + 1) + rr * (log(rr) - log(rr + mu)) + Y * (log(pmax(mu, 1e-300)) - log(rr + mu))
}

# Saturated log-likelihood (mu = y) of every entry.
.rff_loglik_sat <- function(Y, family, r = NULL, n = NULL) {
  if (family == "binomial") {
    p <- ifelse(n > 0, Y / pmax(n, 1e-12), 0.5)
    return(lchoose(n, Y) + ifelse(Y > 0, Y * log(p), 0) + ifelse(n - Y > 0, (n - Y) * log(1 - p), 0))
  }
  if (family == "poisson") return(ifelse(Y > 0, Y * log(Y) - Y, 0) - lgamma(Y + 1))
  rr <- matrix(r, nrow(Y), ncol(Y), byrow = TRUE)
  lgamma(Y + rr) - lgamma(rr) - lgamma(Y + 1) + rr * (log(rr) - log(rr + Y)) + ifelse(Y > 0, Y * (log(Y) - log(rr + Y)), 0)
}

# Deviance of every entry: Poisson deviance for count families (not
# dependent on a dispersion estimate), binomial deviance for binomial.
.rff_deviance <- function(Y, eta, family, n = NULL) {
  if (family == "binomial") {
    p <- stats::plogis(eta); m <- n * p
    return(2 * (ifelse(Y > 0, Y * log(Y / pmax(m, 1e-300)), 0) +
                ifelse(n - Y > 0, (n - Y) * log((n - Y) / pmax(n - m, 1e-300)), 0)))
  }
  mu <- exp(pmin(eta, 40))
  2 * (ifelse(Y > 0, Y * log(Y / pmax(mu, 1e-300)), 0) - (Y - mu))
}

# Assign in-tissue bins to square spatial units of side `size` (coordinate
# units); returns an integer unit id per in-tissue bin.
.rff_units <- function(grid, keep, size) {
  s <- max(1L, as.integer(round(size / grid$bin_size)))
  i <- (which(keep) - 1L) %% grid$nx; j <- (which(keep) - 1L) %/% grid$nx
  u <- (i %/% s) + (j %/% s) * (grid$nx %/% s + 1L)
  as.integer(factor(u))
}

# Vertex of the parabola through three points (x on log scale), clamped.
.rff_parabola_max <- function(x, y, i) {
  if (i <= 1 || i >= length(x)) return(x[i])
  x1 <- x[i - 1]; x2 <- x[i]; x3 <- x[i + 1]; y1 <- y[i - 1]; y2 <- y[i]; y3 <- y[i + 1]
  den <- (x2 - x1) * (y2 - y3) - (x2 - x3) * (y2 - y1)
  if (!is.finite(den) || abs(den) < 1e-12) return(x2)
  v <- x2 - 0.5 * ((x2 - x1)^2 * (y2 - y3) - (x2 - x3)^2 * (y2 - y1)) / den
  min(max(v, x1), x3)
}

.rff_mclapply <- function(X, FUN, n_cores) {
  n_cores <- max(1L, as.integer(n_cores))
  if (n_cores > 1 && .Platform$OS.type == "windows") n_cores <- 1L
  if (n_cores == 1) return(lapply(X, FUN))
  res <- parallel::mclapply(X, FUN, mc.cores = n_cores, mc.preschedule = FALSE)
  bad <- which(vapply(res, function(r) is.null(r) || inherits(r, "try-error"), TRUE))
  for (b in bad) res[[b]] <- FUN(X[[b]])            # recompute lost jobs serially
  res
}

# ---- length-scale profile ----------------------------------------------------

#' Length scale of each spatial program by held-out likelihood
#'
#' **Experimental.** [fit_spatial_rff()] does not estimate length scales:
#' learned length scales stay near their initial values, and its MAP
#' objective favours the smallest one (in spike-in tests the best objective
#' was always at the shortest length scale). `rff_lengthscale_profile()`
#' profiles each program instead. For program `k` of `fit`, everything else
#' in the fit (offset, gene intercepts, cellularity field, the other
#' factors) is held fixed as an offset, the program's loadings are held
#' fixed, and only its field (plus gene intercepts and dispersions) is
#' refitted at every length scale of `ls_grid`, with the in-tissue bins split
#' into `folds` groups: each group is left out once (likelihood weight 0) and
#' the fitted field predicts its counts from the neighbouring bins. The
#' profile is the held-out log-likelihood gain over the same model without
#' the program; its maximum (refined by a parabola in log length scale) is
#' the estimate. A spatial block bootstrap of the per-bin held-out gains
#' gives a 95% interval and standard errors of the gains; the interval
#' reflects the uncertainty of the held-out criterion for this field
#' realisation only and was much narrower than the spread between
#' replicate spike-ins.
#'
#' Held-out units of 3 x 3 bins (default) make the field interpolate over a
#' short distance. Single held-out bins compressed the estimates towards
#' 20-40 um in spike-ins (every smooth field interpolates one bin well);
#' larger units favour longer length scales and became unstable (estimates
#' at the grid edge). The other factors and
#' the cellularity field come from the full-data fit and therefore also saw
#' the held-out bins; this is the same for every length scale and does not
#' move the maximum, but it makes the gains conservative.
#'
#' Validation (6-gene response programs spiked into a real 0.8 x 0.8 mm
#' RA Xenium region around hidden producer cells, decaying as
#' \eqn{e^{-d/\lambda}}, 8 um bins, 2 seeds per \eqn{\lambda}): Spearman
#' correlation of the estimate with \eqn{\lambda} = 10-320 um was 0.89, but
#' the estimates are strongly compressed (about 25 um at \eqn{\lambda} =
#' 10, 40-70 um at \eqn{\lambda} = 320): they rank reaches, they do not
#' measure them. The estimate is the RBF length scale that best
#' interpolates the field, not an exponential decay length, and scales
#' beyond about a tenth of the window are poorly identified.
#'
#' @param fit Output of [fit_spatial_rff()] (preferably `basis = "grid"`).
#' @param binned The `binned_transcripts` object used for the fit.
#' @param factors Factors to profile (names or indices of `fit$L`). Default:
#'   factors with program strength at least 10% of the largest and a
#'   loading pattern not dominated by a part shared by all genes (uniform
#'   share <= 0.5).
#' @param ls_grid Length scales (coordinate units) to profile.
#' @param folds Number of held-out groups.
#' @param holdout_size Side (coordinate units) of the square units that are
#'   held out together; default 3 bins. For programs whose length scale is
#'   close to the bin size, use single bins (`holdout_size` = bin size): a
#'   field cannot predict across a block much larger than its length scale.
#' @param max_iter Iterations for each refit.
#' @param n_boot Block-bootstrap replicates for the intervals.
#' @param boot_block Side (coordinate units) of the bootstrap blocks
#'   (default 8 bins).
#' @param gene_share Only genes that carry this share of the program's
#'   squared loading norm (at least 5 genes) enter the refits; the others
#'   hardly inform the field.
#' @param seed Random seed (fold assignment, bootstrap).
#' @param n_cores Cores for the refits (forked; serial on Windows).
#' @param reach_kernel Decay kernel for the implied reach ([rff_reach()]):
#'   `"exponential"` (default), `"gaussian"`, or `NULL` (no reach).
#'
#' @return An object of class `rff_ls_profile`: `summary` (one row per
#'   profiled factor: `factor`, `lengthscale` (estimate), `lower`, `upper`
#'   (95% bootstrap interval), `grid_best`, `gain` (held-out log-likelihood
#'   gain at the best grid length scale), `gain_se`, `heldout_dev_explained`
#'   (share of the held-out deviance of the program genes explained by the
#'   program), `at_boundary` (maximum at the edge of `ls_grid`), `n_genes`,
#'   `top_genes`, `reach` / `reach_lower` / `reach_upper` (decay length
#'   implied by the length scale, [rff_reach()]) and `identifiable` (the
#'   length scale is at least 2 bins, at most a fifth of the smaller side of
#'   the tissue, and not at the grid edge; otherwise the window or the bins
#'   limit it and it is a bound, not an estimate)), `curves` (one row per factor and length scale: `gain`,
#'   `se`, `delta_se` (standard error of the difference to the best length
#'   scale), `within_1se`, `heldout_dev_explained`), `boot` (bootstrap
#'   estimates) and `settings`.
#' @seealso [fit_spatial_rff()] (`lengthscales = "profile"`),
#'   [rff_programs()]
#' @export
#' @examples
#' \donttest{
#' set.seed(1); n <- 30; res <- 5
#' f1 <- cohalu:::.grf_fft(n, n, res, 15)
#' eta <- log(0.3) + outer(as.vector(f1), c(rep(0.8, 4), rep(0, 6)))
#' Y <- matrix(rpois(length(eta), exp(eta) * res^2), n * n, 10, dimnames = list(NULL, paste0("g", 1:10)))
#' gx <- rep(seq_len(n) - 1, times = n); gy <- rep(seq_len(n) - 1, each = n)
#' b <- structure(list(counts = Matrix::Matrix(Y, sparse = TRUE),
#'   coords = data.frame(x = (gx + 0.5) * res, y = (gy + 0.5) * res, in_tissue = TRUE),
#'   grid = list(nx = n, ny = n, bin_size = res, xmin = 0, ymin = 0), genes = colnames(Y)),
#'   class = "binned_transcripts")
#' fit <- fit_spatial_rff(b, n_factors = 1, lengthscales = 10, basis = "grid",
#'                        learn_lengthscales = FALSE, factor_init = "residual_pca", max_iter = 50)
#' pr <- rff_lengthscale_profile(fit, b, ls_grid = c(5, 10, 20, 40), folds = 3, max_iter = 50, n_boot = 50)
#' pr
#' }
rff_lengthscale_profile <- function(fit, binned, factors = NULL, ls_grid = 5 * 2^(0:6), folds = 4,
                                    holdout_size = NULL, max_iter = 150, n_boot = 200, boot_block = NULL,
                                    gene_share = 0.95, seed = 1, n_cores = 1, reach_kernel = "exponential") {
  if (!inherits(fit, "spatial_rff_fit")) stop("`fit` must come from fit_spatial_rff().")
  if (!inherits(binned, "binned_transcripts")) stop("`binned` must come from bin_transcripts().")
  if (identical(fit$offset, "smoothed_total")) stop("rff_lengthscale_profile() needs a fit with offset = \"area\" or a matrix offset.")
  ls_grid <- sort(unique(as.numeric(ls_grid)))
  if (length(ls_grid) < 2 || any(!is.finite(ls_grid)) || any(ls_grid <= 0)) stop("`ls_grid` needs >= 2 positive length scales.")
  folds <- as.integer(folds); if (folds < 2) stop("`folds` must be >= 2.")
  keep <- fit$in_tissue; g <- fit$grid; N <- sum(keep)
  if (length(keep) != nrow(binned$coords)) stop("`binned` does not match the fit.")
  genes <- fit$genes; L <- fit$L
  fam <- fit$family; if (fam == "multinomial") fam <- "poisson"   # bin effects of the full fit are in the offset
  strength <- sqrt(colSums(L^2)); pstr <- sqrt(colSums(sweep(L, 2, colMeans(L))^2))
  uni <- 1 - pstr^2 / pmax(strength^2, 1e-12)
  if (is.null(factors)) factors <- colnames(L)[pstr >= 0.1 * max(pstr) & uni <= 0.5]
  if (is.numeric(factors)) factors <- colnames(L)[factors]
  factors <- intersect(factors, colnames(L))
  if (!length(factors)) stop("No factor to profile.")
  .local_seed(seed)
  hs <- if (is.null(holdout_size)) 3 * g$bin_size else holdout_size
  unit <- .rff_units(g, keep, hs)
  fold <- sample(rep_len(seq_len(folds), max(unit)))[unit]
  cxy <- binned$coords[keep, c("x", "y")]
  extent <- min(diff(range(cxy$x)), diff(range(cxy$y))) + g$bin_size
  bb <- if (is.null(boot_block)) 8 * g$bin_size else boot_block
  bunit <- .rff_units(g, keep, bb)
  Yall <- as.matrix(binned$counts[keep, genes, drop = FALSE])
  la <- if (fam == "binomial") 0 else 2 * log(g$bin_size)
  tasks <- list()
  for (k in factors) for (f in seq_len(folds)) for (l in c(NA, ls_grid)) tasks[[length(tasks) + 1]] <- list(k = k, f = f, l = l)
  # genes used for each factor
  gsel <- lapply(stats::setNames(factors, factors), function(k) {
    o <- order(-L[, k]^2); cs <- cumsum(L[o, k]^2) / sum(L[, k]^2)
    genes[o[seq_len(max(min(5, length(o)), which(cs >= gene_share)[1]))]]
  })
  run <- function(t) {
    gk <- gsel[[t$k]]
    eta_rest <- .rff_linpred(fit, setdiff(colnames(L), t$k), genes = gk)
    off <- eta_rest - la
    bs <- binned; bs$counts <- binned$counts[, gk, drop = FALSE]; bs$genes <- gk
    W <- matrix(as.numeric(fold != t$f), N, length(gk))
    Lk <- matrix(if (is.na(t$l)) 0 else L[gk, t$k], nrow = length(gk), ncol = 1, dimnames = list(gk, t$k))
    ff <- fit_spatial_rff(bs, n_factors = 1, lengthscales = if (is.na(t$l)) ls_grid[1] else t$l,
                          density_lengthscale = NULL, offset = off, family = fam, basis = "grid",
                          learn_lengthscales = FALSE, loadings = Lk, weights = W, max_iter = max_iter,
                          factor_init = "residual_pca", seed = seed,
                          trials = if (fam == "binomial") fit$trials[, gk, drop = FALSE] else NULL)
    ho <- fold == t$f
    eta <- la + off[ho, , drop = FALSE] + matrix(ff$alpha, sum(ho), length(gk), byrow = TRUE) +
      ff$field_grid[keep, 1][ho] %*% t(ff$L[, 1])
    ll <- .rff_loglik(Yall[ho, gk, drop = FALSE], eta, fam, r = ff$dispersion,
                      n = if (fam == "binomial") fit$trials[ho, gk, drop = FALSE] else NULL)
    list(ll = rowSums(ll), dev = rowSums(.rff_deviance(Yall[ho, gk, drop = FALSE], eta, fam,
                                                    n = if (fam == "binomial") fit$trials[ho, gk, drop = FALSE] else NULL)),
         idx = which(ho))
  }
  res <- .rff_mclapply(tasks, run, n_cores)
  curves <- list(); summ <- list(); boots <- list()
  lg <- log(ls_grid)
  for (k in factors) {
    Lb <- Db <- matrix(NA_real_, N, length(ls_grid) + 1)
    for (i in seq_along(tasks)) {
      t <- tasks[[i]]; if (t$k != k) next
      col <- if (is.na(t$l)) 1 else 1 + match(t$l, ls_grid)
      Lb[res[[i]]$idx, col] <- res[[i]]$ll; Db[res[[i]]$idx, col] <- res[[i]]$dev
    }
    G <- Lb[, -1, drop = FALSE] - Lb[, 1]                 # per-bin held-out gain over no program
    gain <- colSums(G)
    ib <- which.max(gain)
    est <- exp(.rff_parabola_max(lg, gain, ib))
    # block bootstrap over bins
    Gb <- rowsum(G, bunit); nb <- nrow(Gb)
    bs_gain <- vapply(seq_len(n_boot), function(b) colSums(Gb[sample.int(nb, nb, replace = TRUE), , drop = FALSE]),
                      numeric(length(ls_grid)))
    bs_gain <- matrix(bs_gain, nrow = length(ls_grid))
    bs_est <- apply(bs_gain, 2, function(v) exp(.rff_parabola_max(lg, v, which.max(v))))
    se <- apply(bs_gain, 1, stats::sd)
    dse <- apply(sweep(bs_gain, 2, bs_gain[ib, ], "-"), 1, stats::sd)
    dd <- colSums(Db)
    dexp <- (dd[1] - dd[-1]) / max(dd[1], 1e-12)          # held-out (Poisson / binomial) deviance explained
    curves[[k]] <- data.frame(factor = k, lengthscale = ls_grid, gain = gain, se = se, delta_se = dse,
                              within_1se = gain >= gain[ib] - dse, heldout_dev_explained = dexp, row.names = NULL)
    o <- order(-abs(L[, k]))[seq_len(min(5, length(genes)))]
    summ[[k]] <- data.frame(factor = k, lengthscale = est,
                            lower = unname(stats::quantile(bs_est, 0.025)), upper = unname(stats::quantile(bs_est, 0.975)),
                            grid_best = ls_grid[ib], gain = gain[ib], gain_se = se[ib],
                            heldout_dev_explained = dexp[ib], at_boundary = ib %in% c(1, length(ls_grid)),
                            n_genes = length(gsel[[k]]),
                            reach = NA_real_, reach_lower = NA_real_, reach_upper = NA_real_,
                            identifiable = est >= 2 * g$bin_size && est <= extent / 5 && !(ib %in% c(1, length(ls_grid))),
                            top_genes = paste(sprintf("%s (%+.2f)", genes[o], L[o, k]), collapse = ", "), row.names = NULL)
    boots[[k]] <- bs_est
  }
  sm <- do.call(rbind, summ)
  if (!is.null(reach_kernel)) {
    sm$reach <- rff_reach(sm$lengthscale, g$bin_size, reach_kernel)
    sm$reach_lower <- rff_reach(sm$lower, g$bin_size, reach_kernel)
    sm$reach_upper <- rff_reach(sm$upper, g$bin_size, reach_kernel)
  }
  structure(list(summary = sm, curves = do.call(rbind, curves), boot = boots,
                 settings = list(ls_grid = ls_grid, folds = folds, holdout_size = hs, max_iter = max_iter,
                                 n_boot = n_boot, boot_block = bb, gene_share = gene_share, seed = seed)),
            class = "rff_ls_profile")
}

#' @export
print.rff_ls_profile <- function(x, ...) {
  s <- x$summary
  cat(sprintf("<rff_ls_profile> %d program(s); grid %s; %d folds, held-out units %g\n", nrow(s),
              paste(signif(x$settings$ls_grid, 3), collapse = "/"), x$settings$folds, x$settings$holdout_size))
  for (i in seq_len(nrow(s)))
    cat(sprintf("  %-8s length scale %.3g (95%% CI %.3g-%.3g)%s%s | held-out gain %.1f (se %.1f), dev. explained %.2f%% | %s\n",
                s$factor[i], s$lengthscale[i], s$lower[i], s$upper[i], if (s$at_boundary[i]) " [grid edge]" else "",
                if (!is.null(s$reach) && is.finite(s$reach[i])) sprintf(", reach %.3g%s", s$reach[i], if (isTRUE(s$identifiable[i])) "" else " (not identifiable)") else "",
                s$gain[i], s$gain_se[i], 100 * s$heldout_dev_explained[i], s$top_genes[i]))
  invisible(x)
}

# fit_spatial_rff(lengthscales = "profile"): first fit, profile, refit.
.rff_fit_profiled <- function(args, ls_grid, profile) {
  init <- if (is.null(profile$init)) 20 else profile$init
  pargs <- profile[setdiff(names(profile), "init")]
  if (isTRUE(args$learn_lengthscales)) args$learn_lengthscales <- FALSE
  base <- do.call(fit_spatial_rff, c(args[setdiff(names(args), "start")], list(lengthscales = init)))
  pr <- do.call(rff_lengthscale_profile, c(list(base, args$binned, ls_grid = ls_grid), pargs))
  K <- ncol(base$L)
  ell <- stats::setNames(rep_len(init, K), colnames(base$L))
  ell[pr$summary$factor] <- pr$summary$lengthscale
  final <- do.call(fit_spatial_rff, c(args[setdiff(names(args), "start")],
                                      list(lengthscales = unname(ell), start = base)))
  final$lengthscale_profile <- pr
  final$settings$lengthscale_method <- "profile"
  final$settings$lengthscale_init <- init
  final
}

# ---- cross-fitted (held-out) program statistics ------------------------------

# Two halves of the in-tissue bins from a checkerboard of square blocks of
# side `block` (coordinate units). Bins closer than `margin` bins to the other
# half are not used for evaluation (smoothing leaks across the border).
.rff_halves <- function(grid, keep, block, margin) {
  s <- max(2L, as.integer(round(block / grid$bin_size)))
  i <- (which(keep) - 1L) %% grid$nx; j <- (which(keep) - 1L) %/% grid$nx
  half <- ((i %/% s) + (j %/% s)) %% 2L
  unit <- as.integer(factor((i %/% s) + (j %/% s) * (grid$nx %/% s + 1L)))
  full <- matrix(-1L, grid$nx, grid$ny); full[which(keep)] <- half
  m <- max(0L, as.integer(margin))
  near <- function(h) if (m == 0) matrix(FALSE, grid$nx, grid$ny) else .dilate(full == h, m)
  nearA <- near(0L)[which(keep)]; nearB <- near(1L)[which(keep)]
  list(A = half == 0L, B = half == 1L, evalA = half == 0L & !nearB, evalB = half == 1L & !nearA, unit = unit)
}

# Cross-fitted variance shares: components found on one half (top-k right
# singular vectors of R there) and their share of the residual variance of
# the other half (evaluation bins only); both directions averaged. Returns
# the shares and per-unit numerators/denominators for a block bootstrap.
.rff_cf_shares <- function(R, hv, k) {
  one <- function(disc, ev) {
    V <- svd(R[disc, , drop = FALSE], nu = 0, nv = k)$v
    P <- (R[ev, , drop = FALSE] %*% V)^2
    list(num = rowsum(P, hv$unit[ev]), den = rowsum(rowSums(R[ev, , drop = FALSE]^2), hv$unit[ev]), V = V)
  }
  ab <- one(hv$A, hv$evalB); ba <- one(hv$B, hv$evalA)
  sh <- function(x) colSums(x$num) / sum(x$den)
  list(share = (sh(ab) + sh(ba)) / 2, ab = ab, ba = ba)
}

# Block bootstrap of the cross-fitted shares (units resampled within each half).
.rff_cf_boot <- function(cf, n_boot) {
  rs <- function(x) { n <- nrow(x$num); i <- sample.int(n, n, replace = TRUE); colSums(x$num[i, , drop = FALSE]) / sum(x$den[i]) }
  matrix(vapply(seq_len(n_boot), function(b) (rs(cf$ab) + rs(cf$ba)) / 2, numeric(ncol(cf$ab$num))), ncol = n_boot)
}

# Cross-fitted statistics of one fitted window: observed held-out shares, the
# same statistic for gene-shift surrogates, relative excess and its bootstrap
# interval, plus the residual covariance of the whole window and the mean
# surrogate covariance (for direction-specific comparisons with controls).
.rff_cf_window <- function(fit, binned, bandwidth, highpass, n_components, n_boot, block, n_boot_ci = 200,
                           keep_cov = FALSE, robust = TRUE) {
  sp <- .rff_shift_prep(fit, binned, bandwidth, highpass, robust = robust)
  k <- max(1L, min(as.integer(n_components), sp$J))
  h <- (if (is.null(bandwidth)) fit$bin_size else bandwidth) / fit$bin_size
  if (is.null(block)) block <- 12 * fit$bin_size
  hv <- .rff_halves(fit$grid, fit$in_tissue, block, margin = ceiling(1.5 * h))
  if (sum(hv$evalA) < 20 || sum(hv$evalB) < 20) stop("Too few bins for cross-fitting; use a smaller `crossfit_block`.")
  R <- sp$process(sp$Y, sp$mu, sp$v)
  cf <- .rff_cf_shares(R, hv, k)
  nul <- matrix(0, n_boot, k); Cbar <- if (keep_cov) 0 else NULL
  for (b in seq_len(n_boot)) {
    z <- sp$surrogate(sp$mu, sp$v); Rz <- sp$process(z$Y, z$mu, z$v)
    nul[b, ] <- .rff_cf_shares(Rz, hv, k)$share
    if (keep_cov) { Cz <- crossprod(Rz) / sum(Rz^2); Cbar <- Cbar + Cz / n_boot }
  }
  bs <- .rff_cf_boot(cf, n_boot_ci)
  m <- colMeans(nul)
  list(share = cf$share, null = nul, null_mean = m, excess = cf$share / m - 1,
       excess_boot = bs / m - 1, C = if (keep_cov) crossprod(R) / sum(R^2) else NULL, Cbar = Cbar,
       V = if (keep_cov) svd(R, nu = 0, nv = k)$v else NULL,
       R = R, N = sp$N, genes = fit$genes, block = block)
}

#' Negative-control reference for calibrated program tests
#'
#' **Experimental.** On real tissue, the gene-shift null of
#' [rff_program_test()] is rejected by almost any shared residual structure
#' (segmentation spill-over, mixed bins, cellularity, cell-state
#' heterogeneity), so "significant" does not separate a program of interest
#' from the structure that every region of the tissue has. A within-data null
#' makes the question relative: is the target's shared residual structure
#' stronger (rank-wise), or is a given program direction more active, than
#' in regions where no program of interest is expected (negative controls -
#' e.g. regions without TLS, a homogeneous monolayer far from the
#' inducing cells, untreated wells)?
#'
#' `rff_control_reference()` runs the cross-fitted statistic of
#' `rff_program_test(crossfit = TRUE)` in every control window: components
#' are found on one half of the bins (checkerboard of spatial blocks) and
#' their share of the residual variance is measured on the other half, for
#' the data and for gene-shift surrogates; the relative excess
#' (held-out share / surrogate mean - 1) of every component rank is stored,
#' together with each window's residual covariance and mean surrogate
#' covariance (for direction-specific comparisons). Pass the result as
#' `control` to [rff_program_test()].
#'
#' @param fits List of [fit_spatial_rff()] fits of the control windows (same
#'   genes as the target fits).
#' @param binned List of the matching `binned_transcripts` objects.
#' @param bandwidth,highpass,n_components,n_boot,crossfit_block As in
#'   [rff_program_test()]; use the same values for the target test.
#' @param seed Random seed.
#' @param n_cores Cores (windows are processed in parallel; forked).
#' @param keep_matrices Also store each control window's cross-fitting
#'   matrices (a second pass of surrogates), needed by
#'   `rff_program_test_joint(control = )`.
#'
#' @return An object of class `rff_control`: `excess` (controls x
#'   components), `share`, `null_mean`, per-window covariances `C` and
#'   `Cbar` and component directions `V`, `loo` (leave-one-out calibration:
#'   every control window tested against the others - held-out p-value,
#'   rank-wise `p_control`, direction-wise `z_dir` / `p_dir`; `called` /
#'   `called_dir` = held-out p and the control p <= 0.05, i.e. the
#'   false-positive rate of the procedure among negative controls), `genes`
#'   and `settings`.
#' @seealso [rff_program_test()]
#' @export
rff_control_reference <- function(fits, binned, bandwidth = NULL, highpass = NULL, n_components = 6, n_boot = 49,
                                  crossfit_block = NULL, seed = 1, n_cores = 1, keep_matrices = TRUE) {
  if (inherits(fits, "spatial_rff_fit")) fits <- list(fits)
  if (inherits(binned, "binned_transcripts")) binned <- list(binned)
  if (length(fits) < 3 || length(fits) != length(binned)) stop("`fits` and `binned` must be lists of >= 3 control windows of the same length.")
  genes <- fits[[1]]$genes
  if (!all(vapply(fits, function(f) identical(f$genes, genes), TRUE))) stop("All control fits must have the same genes in the same order.")
  wn <- names(fits); if (is.null(wn)) wn <- paste0("control", seq_along(fits))
  .local_seed(seed); seeds <- sample.int(.Machine$integer.max, length(fits))
  res <- .rff_mclapply(seq_along(fits), function(i) {
    set.seed(seeds[i])
    w <- .rff_cf_window(fits[[i]], binned[[i]], bandwidth, highpass, n_components, n_boot, crossfit_block, n_boot_ci = 0, keep_cov = TRUE)
    if (keep_matrices) {                                  # for joint (group) tests
      hv <- .rff_halves(fits[[i]]$grid, fits[[i]]$in_tissue, w$block,
                        margin = ceiling(1.5 * (if (is.null(bandwidth)) fits[[i]]$bin_size else bandwidth) / fits[[i]]$bin_size))
      w$M <- .rff_cf_window_mats(fits[[i]], binned[[i]], bandwidth, highpass, n_boot, w$block)$m
      w$M$C <- NULL; w$M$Cbar <- NULL
    }
    w$R <- NULL; w$excess_boot <- NULL; w
  }, n_cores)
  k <- length(res[[1]]$share)
  E <- t(vapply(res, `[[`, numeric(k), "excess")); if (k == 1) E <- matrix(E, ncol = 1)
  dimnames(E) <- list(wn, paste0("PC", seq_len(k)))
  ph <- t(vapply(res, function(w) cummax(vapply(seq_len(k), function(j) (1 + sum(w$null[, j] >= w$share[j])) / (nrow(w$null) + 1), 0)),
                 numeric(k))); if (k == 1) ph <- matrix(ph, ncol = 1)
  out <- structure(list(excess = E, share = t(vapply(res, `[[`, numeric(k), "share")),
                        null_mean = t(vapply(res, `[[`, numeric(k), "null_mean")), p_heldout = ph,
                        C = lapply(res, `[[`, "C"), Cbar = lapply(res, `[[`, "Cbar"), V = lapply(res, `[[`, "V"),
                        N = vapply(res, `[[`, 1, "N"), genes = genes,
                        M = if (keep_matrices) lapply(res, `[[`, "M") else NULL,
                        settings = list(bandwidth = bandwidth, highpass = highpass, n_components = k, n_boot = n_boot,
                                        crossfit_block = res[[1]]$block)),
                   class = "rff_control")
  out$loo <- .rff_control_loo(out)
  out
}

# Subset of the control windows of a reference (e.g. leave one out).
.rff_control_subset <- function(ref, keep) {
  r <- ref
  r$excess <- ref$excess[keep, , drop = FALSE]; r$share <- ref$share[keep, , drop = FALSE]
  r$null_mean <- ref$null_mean[keep, , drop = FALSE]; r$p_heldout <- ref$p_heldout[keep, , drop = FALSE]
  r$C <- ref$C[keep]; r$Cbar <- ref$Cbar[keep]; r$V <- ref$V[keep]; r$N <- ref$N[keep]
  if (!is.null(ref$M)) r$M <- ref$M[keep]
  r$loo <- .rff_control_loo(r)
  r
}

# Direction-specific z of a target (excess `e` along unit direction `v`)
# against control windows `idx`.
.rff_zdir <- function(ref, v, e, idx = seq_along(ref$C)) {
  ec <- vapply(idx, function(i) drop(crossprod(v, ref$C[[i]] %*% v)) / max(drop(crossprod(v, ref$Cbar[[i]] %*% v)), 1e-12) - 1, 0)
  list(z = (e - mean(ec)) / max(stats::sd(ec), 1e-12), e = ec)
}

# Leave-one-out calibration among the controls: each control window is tested
# against the others, rank-wise (p_control) and direction-wise (z_dir: its
# excess along its own component direction vs the others' excess along that
# direction; p_dir: its z_dir among the other controls' z_dir, i.e. the
# calibrated direction-specific p-value).
.rff_control_loo <- function(ref) {
  E <- ref$excess; m <- nrow(E); k <- ncol(E); wn <- rownames(E)
  pc <- matrix(vapply(seq_len(k), function(j) vapply(seq_len(m), function(i) (1 + sum(E[-i, j] >= E[i, j])) / m, 0), numeric(m)), m, k)
  Z <- matrix(NA_real_, m, k)
  if (!is.null(ref$V[[1]])) for (i in seq_len(m)) for (j in seq_len(k)) Z[i, j] <- .rff_zdir(ref, ref$V[[i]][, j], E[i, j], setdiff(seq_len(m), i))$z
  pd <- matrix(vapply(seq_len(k), function(j) vapply(seq_len(m), function(i) (1 + sum(Z[-i, j] >= Z[i, j])) / m, 0), numeric(m)), m, k)
  ph <- ref$p_heldout
  data.frame(window = rep(wn, k), component = rep(colnames(E), each = m), excess = as.vector(E),
             p_heldout = as.vector(ph), p_control = as.vector(pc), z_dir = as.vector(Z), p_dir = as.vector(pd),
             called = as.vector(ph <= 0.05 & pc <= 0.05), called_dir = as.vector(ph <= 0.05 & pd <= 0.05), row.names = NULL)
}

#' @export
print.rff_control <- function(x, ...) {
  cat(sprintf("<rff_control> %d negative-control windows, %d genes, %d components\n", nrow(x$excess), length(x$genes), ncol(x$excess)))
  q <- apply(x$excess, 2, stats::quantile, c(0.5, 0.95))
  cat("  relative excess of held-out share over gene-shift surrogates (median / 95th percentile):\n")
  cat("   ", paste(sprintf("%s %.2f/%.2f", colnames(x$excess), q[1, ], q[2, ]), collapse = "  "), "\n")
  l <- x$loo
  pc1 <- l[l$component == "PC1", ]
  cat(sprintf("  leave-one-out (first component): held-out test alone called %d/%d control windows; rank-wise control comparison %d/%d; direction-wise %d/%d\n",
              sum(pc1$p_heldout <= 0.05), nrow(pc1), sum(pc1$called), nrow(pc1), sum(pc1$called_dir), nrow(pc1)))
  invisible(x)
}

# Internal: rff_program_test(null = "shift", crossfit = TRUE [, control]).
.rff_program_test_cf <- function(fit, binned, bandwidth, n_components, n_boot, sequential, alpha, seed,
                                 highpass, crossfit_block, control, min_effect, n_boot_ci, control_stat = "direction") {
  .local_seed(seed)
  if (!is.null(control)) {
    if (!inherits(control, "rff_control")) stop("`control` must come from rff_control_reference().")
    if (!identical(control$genes, fit$genes)) stop("`control` was built for other genes.")
    cs <- control$settings
    if (!identical(cs$highpass, highpass) || !identical(cs$bandwidth, bandwidth))
      warning("`highpass` / `bandwidth` differ from those of the control reference.")
    if (is.null(crossfit_block)) crossfit_block <- cs$crossfit_block
    n_components <- min(n_components, ncol(control$excess))
  }
  w <- .rff_cf_window(fit, binned, bandwidth, highpass, n_components, n_boot, crossfit_block, n_boot_ci = n_boot_ci)
  k <- length(w$share)
  p <- vapply(seq_len(k), function(j) (1 + sum(w$null[, j] >= w$share[j])) / (n_boot + 1), 0)
  if (sequential) p <- cummax(p)
  lo <- apply(w$excess_boot, 1, stats::quantile, 0.025); hi <- apply(w$excess_boot, 1, stats::quantile, 0.975)
  # full-window components for display, matching and scores
  sv <- svd(w$R, nu = k, nv = k)
  genes <- fit$genes; Fm <- fit$field_grid[fit$in_tissue, colnames(fit$L), drop = FALSE]
  cc <- abs(stats::cor(sv$u, Fm)); cc[!is.finite(cc)] <- 0
  top <- vapply(seq_len(k), function(j) {
    o <- order(-abs(sv$v[, j]))[seq_len(min(5, length(genes)))]
    paste(sprintf("%s (%+.2f)", genes[o], sv$v[o, j]), collapse = ", ")
  }, "")
  out <- data.frame(component = paste0("PC", seq_len(k)), share = (sv$d^2 / sum(w$R^2))[seq_len(k)], p = p,
                    factor = colnames(Fm)[apply(cc, 1, which.max)], r_factor = apply(cc, 1, max), top_genes = top,
                    share_heldout = w$share, null_heldout = w$null_mean, excess = w$excess,
                    excess_lower = lo, excess_upper = hi, z_heldout = (w$share - w$null_mean) / pmax(apply(w$null, 2, stats::sd), 1e-12),
                    row.names = NULL)
  call <- out$p <= alpha
  thr <- if (is.null(min_effect)) 0 else min_effect
  if (!is.null(control)) {
    E <- control$excess[, seq_len(k), drop = FALSE]
    out$p_control <- vapply(seq_len(k), function(j) (1 + sum(E[, j] >= out$excess[j])) / (nrow(E) + 1), 0)
    out$control_q95 <- apply(E, 2, stats::quantile, 0.95)
    zd <- lapply(seq_len(k), function(j) .rff_zdir(control, sv$v[, j], out$excess[j]))
    out$z_control_dir <- vapply(zd, `[[`, 0, "z")
    Zc <- matrix(control$loo$z_dir, nrow(E))
    out$p_control_dir <- vapply(seq_len(k), function(j) (1 + sum(Zc[, j] >= out$z_control_dir[j])) / (nrow(E) + 1), 0)
    call <- call & (if (control_stat == "rank") out$p_control else out$p_control_dir) <= alpha
  }
  out$effect_threshold <- thr
  call <- call & out$excess_lower > thr
  if (sequential) call <- cumprod(call) == 1                  # stop at the first component not called
  out$call <- call
  scores <- sv$u; colnames(scores) <- out$component
  loadings <- sv$v; dimnames(loadings) <- list(genes, out$component)
  attr(out, "scores") <- scores; attr(out, "loadings") <- loadings
  attr(out, "null_shares") <- lapply(seq_len(k), function(j) w$null[, j])
  attr(out, "null") <- "shift"; attr(out, "crossfit") <- list(block = w$block, n_boot_ci = n_boot_ci,
                                                               control = !is.null(control), min_effect = min_effect)
  out
}

# ---- transfer of fixed loadings ----------------------------------------------

#' Calibrated test for a gene program transferred to new data
#'
#' **Experimental.** Tests whether gene programs learned elsewhere (fixed
#' loadings, e.g. from another data set, platform or patient group) are
#' active in new tissue, against a **permuted-loading null**: the target
#' data are left untouched, so every gene keeps its counts, its marginal
#' distribution and its own spatial structure, and only the assignment of
#' the loadings to genes is permuted (within strata of gene abundance). A
#' program counts as transferred if its genes co-vary in the target more
#' than random gene sets with the same loading values.
#'
#' Statistics:
#' * `"score"` (default): share of the processed residual variance of the
#'   target (the residuals of the gene-shift test: Pearson residuals under
#'   the fit's null mean, smoothed, optional high-pass, shared
#'   multiplicative direction removed, genes standardised) along the
#'   program direction (log-scale loadings times the square root of each
#'   gene's mean count, unit norm). It is the score statistic of a
#'   rank-one program with fixed loadings, fast enough for thousands of
#'   permutations.
#' * `"heldout"`: held-out log-likelihood gain of a fitted program field
#'   with the loadings held fixed (length scale `lengthscale`) over the
#'   model without it, by `folds`-fold cross-validation over bins
#'   (likelihood-ratio type; slow, use fewer permutations).
#'
#' The field amplitude of a fixed-loading fit (`fit_spatial_rff(loadings =
#' )$amplitude`) is not a transfer statistic: permuted loadings reached
#' similar amplitudes on real tissue, because any direction picks up some
#' shared residual structure.
#'
#' @param loadings Numeric matrix, genes (row names) x programs, on the log
#'   scale (e.g. `fit$L` or [fit_spatial_rff_joint()]`$loadings`); a named
#'   vector is one program. Genes absent from the target are dropped.
#' @param fits A [fit_spatial_rff()] fit of the target tissue (its offset
#'   and intercepts give the null mean; the target's own factors are not
#'   used) or a list of fits (windows; statistics are averaged over
#'   windows with one permutation applied to all).
#' @param binned The matching `binned_transcripts` object(s).
#' @param statistic `"score"` or `"heldout"`.
#' @param n_perm Number of permutations.
#' @param strata Number of gene-abundance strata (quantiles of the mean
#'   count per bin in the target) within which loadings are permuted.
#' @param bandwidth,highpass Residual processing for `"score"` (as in
#'   [rff_program_test()]).
#' @param lengthscale,folds,max_iter Field length scale, cross-validation
#'   folds and iterations for `"heldout"`.
#' @param seed Random seed.
#' @param n_cores Cores for `"heldout"` refits.
#' @param scale `"log"` (default): `loadings` are log-scale effects (as in
#'   `fit$L`); `"residual"`: they are directions of the processed residuals,
#'   e.g. the `loadings` attribute of [rff_program_test()] or
#'   [rff_program_test_joint()].
#'
#' @return A data frame with one row per program: `program`, `n_genes`
#'   (genes shared with the target), `statistic`, `perm_mean`, `excess`
#'   (statistic / permutation mean - 1), `z`, `p`; attribute `null`
#'   (programs x permutations).
#' @seealso [fit_spatial_rff()] (`loadings`), [rff_program_test_joint()]
#'   (confirmatory mode, gene-shift null)
#' @export
rff_transfer_test <- function(loadings, fits, binned, statistic = c("score", "heldout"), n_perm = NULL, strata = 5,
                              bandwidth = NULL, highpass = NULL, lengthscale = 20, folds = 4, max_iter = 150,
                              seed = 1, n_cores = 1, scale = c("log", "residual")) {
  statistic <- match.arg(statistic); scale <- match.arg(scale)
  if (inherits(fits, "spatial_rff_fit")) fits <- list(fits)
  if (inherits(binned, "binned_transcripts")) binned <- list(binned)
  if (length(fits) != length(binned)) stop("`fits` and `binned` must match.")
  genes <- fits[[1]]$genes
  if (!all(vapply(fits, function(f) identical(f$genes, genes), TRUE))) stop("All fits must have the same genes in the same order.")
  L <- if (is.null(dim(loadings))) matrix(loadings, ncol = 1, dimnames = list(names(loadings), "program1")) else as.matrix(loadings)
  if (is.null(rownames(L))) stop("`loadings` needs gene row names.")
  if (is.null(colnames(L))) colnames(L) <- paste0("program", seq_len(ncol(L)))
  gs <- intersect(genes, rownames(L))
  if (length(gs) < 3) stop("Fewer than 3 genes of `loadings` are in the target.")
  if (is.null(n_perm)) n_perm <- if (statistic == "score") 999 else 49
  .local_seed(seed)
  # abundance strata in the target (mean count per in-tissue bin, pooled)
  mbar <- Reduce(`+`, lapply(binned, function(b) Matrix::colSums(b$counts[b$coords$in_tissue, genes, drop = FALSE]))) /
    sum(vapply(binned, function(b) sum(b$coords$in_tissue), 1))
  J <- length(genes)
  st <- if (strata > 1) as.integer(cut(rank(mbar, ties.method = "first"), strata, labels = FALSE)) else rep(1L, J)
  Lfull <- matrix(0, J, ncol(L), dimnames = list(genes, colnames(L))); Lfull[gs, ] <- L[gs, , drop = FALSE]
  perms <- lapply(seq_len(n_perm), function(i) { o <- seq_len(J); for (s in unique(st)) { w <- which(st == s); o[w] <- w[sample.int(length(w))] }; o })
  if (statistic == "score") {
    Cs <- lapply(seq_along(fits), function(w) {
      sp <- .rff_shift_prep(fits[[w]], binned[[w]], bandwidth, highpass, robust = TRUE)
      R <- sp$process(sp$Y, sp$mu, sp$v)
      list(C = crossprod(R) / sum(R^2), s = if (scale == "log") sqrt(colMeans(sp$mu)) else rep(1, sp$J))
    })
    stat <- function(v) mean(vapply(Cs, function(x) { u <- v * x$s; u <- u / max(sqrt(sum(u^2)), 1e-12); drop(crossprod(u, x$C %*% u)) }, 0))
    obs <- apply(Lfull, 2, stat)
    nul <- vapply(perms, function(o) apply(Lfull[o, , drop = FALSE], 2, stat), numeric(ncol(L)))
  } else {
    if (scale == "residual") Lfull <- Lfull / sqrt(pmax(mbar, 1e-3))
    gain <- function(v) mean(vapply(seq_along(fits), function(w) .rff_fixed_gain(fits[[w]], binned[[w]], v, lengthscale, folds, max_iter, seed), 0))
    jobs <- c(list(seq_len(J)), perms)
    vals <- .rff_mclapply(jobs, function(o) apply(Lfull[o, , drop = FALSE], 2, gain), n_cores)
    vals <- matrix(unlist(vals), nrow = ncol(L))
    obs <- vals[, 1]; nul <- vals[, -1, drop = FALSE]
  }
  nul <- matrix(nul, nrow = ncol(L))
  pm <- rowMeans(nul); psd <- apply(nul, 1, stats::sd)
  out <- data.frame(program = colnames(L), n_genes = colSums(Lfull != 0), statistic = obs, perm_mean = pm,
                    excess = obs / pm - 1, z = (obs - pm) / pmax(psd, 1e-12),
                    p = (1 + rowSums(nul >= obs)) / (n_perm + 1), row.names = NULL)
  if (statistic == "heldout") out$excess <- NA_real_              # gains can be <= 0: ratio not meaningful
  attr(out, "null") <- nul; attr(out, "statistic") <- statistic
  out
}

# Held-out log-likelihood gain of a field along fixed loadings `v` (log scale,
# genes of the fit) over the null part of the fit (offset, intercepts,
# cellularity field; the fit's own factors are not used).
.rff_fixed_gain <- function(fit, binned, v, ell, folds, max_iter, seed) {
  keep <- fit$in_tissue; g <- fit$grid; N <- sum(keep); fam <- fit$family; if (fam == "multinomial") fam <- "poisson"
  gk <- fit$genes[v != 0]; if (length(gk) < 2) return(0)
  la <- if (fam == "binomial") 0 else 2 * log(g$bin_size)
  off <- .rff_linpred(fit, NULL, genes = gk) - la
  bs <- binned; bs$counts <- binned$counts[, gk, drop = FALSE]; bs$genes <- gk
  Y <- as.matrix(binned$counts[keep, gk, drop = FALSE])
  set.seed(seed); fold <- sample(rep_len(seq_len(folds), N))
  tr <- if (fam == "binomial") fit$trials[, gk, drop = FALSE] else NULL
  tot <- 0
  for (f in seq_len(folds)) {
    W <- matrix(as.numeric(fold != f), N, length(gk)); ho <- fold == f
    ll <- vapply(list(0, v[v != 0]), function(lv) {
      Lk <- matrix(lv, nrow = length(gk), ncol = 1, dimnames = list(gk, "p"))
      ff <- fit_spatial_rff(bs, n_factors = 1, lengthscales = ell, density_lengthscale = NULL, offset = off, family = fam,
                            basis = "grid", learn_lengthscales = FALSE, loadings = Lk, weights = W, max_iter = max_iter,
                            factor_init = "residual_pca", seed = seed, trials = tr)
      eta <- la + off[ho, , drop = FALSE] + matrix(ff$alpha, sum(ho), length(gk), byrow = TRUE) + ff$field_grid[keep, 1][ho] %*% t(ff$L[, 1])
      sum(.rff_loglik(Y[ho, , drop = FALSE], eta, fam, r = ff$dispersion, n = if (is.null(tr)) NULL else tr[ho, , drop = FALSE]))
    }, 0)
    tot <- tot + ll[2] - ll[1]
  }
  tot
}

# ---- cross-fitted statistics as J x J matrices (joint / group tests) ----------

# Per-window matrices of the cross-fitted statistic: covariances of the
# processed residuals on the discovery halves (CA, CB) and on the evaluation
# halves (EA, EB), and the surrogate means of the evaluation-half covariances
# (SA, SB). For a set of windows with weights c_w, components found on
# sum_w c_w CA_w are evaluated on sum_w c_w EB_w (and vice versa); everything
# a group statistic needs is in these matrices.
.rff_cf_mats_one <- function(R, hv) {
  cp <- function(i) crossprod(R[i, , drop = FALSE])
  list(CA = cp(hv$A), CB = cp(hv$B), EA = cp(hv$evalA), EB = cp(hv$evalB))
}

.rff_cf_group_share <- function(M, cw, k) {
  pool <- function(nm) Reduce(`+`, Map(function(m, c) m[[nm]] * c, M, cw))
  CA <- pool("CA"); CB <- pool("CB"); EA <- pool("EA"); EB <- pool("EB")
  VA <- eigen(CA, symmetric = TRUE)$vectors[, seq_len(k), drop = FALSE]
  VB <- eigen(CB, symmetric = TRUE)$vectors[, seq_len(k), drop = FALSE]
  sAB <- colSums(VA * (EB %*% VA)) / sum(diag(EB)); sBA <- colSums(VB * (EA %*% VB)) / sum(diag(EA))
  list(share = (sAB + sBA) / 2, VA = VA, VB = VB)
}

# held-out excess of each window along the group's directions (found on the
# other half), relative to its surrogate mean
.rff_cf_window_excess <- function(m, VA, VB) {
  eB <- colSums(VA * (m$EB %*% VA)) / pmax(colSums(VA * (m$SB %*% VA)), 1e-12)
  eA <- colSums(VB * (m$EA %*% VB)) / pmax(colSums(VB * (m$SA %*% VB)), 1e-12)
  (eA + eB) / 2 - 1
}

# direction excess of full control windows along directions V (columns)
.rff_ctrl_excess <- function(ref, V, idx = seq_along(ref$C))
  t(vapply(idx, function(i) colSums(V * (ref$C[[i]] %*% V)) / pmax(colSums(V * (ref$Cbar[[i]] %*% V)), 1e-12) - 1, numeric(ncol(V))))

# matrices of one window: observed and surrogate means (surrogate draws are
# also returned as pooled-ready lists when keep_draws = TRUE)
.rff_cf_window_mats <- function(fit, binned, bandwidth, highpass, n_boot, block, keep_draws = FALSE) {
  sp <- .rff_shift_prep(fit, binned, bandwidth, highpass, robust = TRUE)
  h <- (if (is.null(bandwidth)) fit$bin_size else bandwidth) / fit$bin_size
  if (is.null(block)) block <- 12 * fit$bin_size
  hv <- .rff_halves(fit$grid, fit$in_tissue, block, margin = ceiling(1.5 * h))
  if (sum(hv$evalA) < 20 || sum(hv$evalB) < 20) stop("Too few bins for cross-fitting; use a smaller `crossfit_block`.")
  R <- sp$process(sp$Y, sp$mu, sp$v)
  m <- .rff_cf_mats_one(R, hv); m$N <- sp$N
  m$C <- crossprod(R) / sum(R^2)
  SA <- SB <- 0; Cbar <- 0; draws <- vector("list", if (keep_draws) n_boot else 0)
  for (b in seq_len(n_boot)) {
    z <- sp$surrogate(sp$mu, sp$v); Rz <- sp$process(z$Y, z$mu, z$v)
    mz <- .rff_cf_mats_one(Rz, hv)
    SA <- SA + mz$EA / n_boot; SB <- SB + mz$EB / n_boot; Cbar <- Cbar + crossprod(Rz) / sum(Rz^2) / n_boot
    if (keep_draws) draws[[b]] <- mz
  }
  m$SA <- SA; m$SB <- SB; m$Cbar <- Cbar; m$block <- block
  list(m = m, draws = draws)
}

# Internal: rff_program_test_joint(crossfit = TRUE [, control]).
.rff_program_test_joint_cf <- function(fits, binned, bandwidth, highpass, n_components, n_boot, sequential, alpha,
                                       seed, cw, wn, crossfit_block, control, n_group_null, n_cores) {
  genes <- fits[[1]]$genes; J <- length(genes)
  k <- max(1L, min(as.integer(n_components), J))
  if (!is.null(control)) {
    if (!inherits(control, "rff_control")) stop("`control` must come from rff_control_reference().")
    if (!identical(control$genes, genes)) stop("`control` was built for other genes.")
    if (is.null(control$M)) stop("`control` lacks the cross-fitting matrices; rebuild it with the current rff_control_reference().")
    if (is.null(crossfit_block)) crossfit_block <- control$settings$crossfit_block
  }
  .local_seed(seed); seeds <- sample.int(.Machine$integer.max, length(fits))
  W <- .rff_mclapply(seq_along(fits), function(i) { set.seed(seeds[i])
    .rff_cf_window_mats(fits[[i]], binned[[i]], bandwidth, highpass, n_boot, crossfit_block, keep_draws = TRUE) }, n_cores)
  M <- lapply(W, `[[`, "m")
  cwn <- cw / sum(cw)
  obs <- .rff_cf_group_share(M, cwn, k)
  nul <- t(vapply(seq_len(n_boot), function(b) .rff_cf_group_share(lapply(W, function(w) w$draws[[b]]), cwn, k)$share, numeric(k)))
  if (k == 1) nul <- matrix(nul, ncol = 1)
  rm(W)
  mn <- colMeans(nul)
  p <- vapply(seq_len(k), function(j) (1 + sum(nul[, j] >= obs$share[j])) / (n_boot + 1), 0)
  if (sequential) p <- cummax(p)
  # per-window held-out excess along the group's directions; window bootstrap for the interval
  E <- t(vapply(M, .rff_cf_window_excess, numeric(k), VA = obs$VA, VB = obs$VB)); if (k == 1) E <- matrix(E, ncol = 1)
  dimnames(E) <- list(wn, paste0("PC", seq_len(k)))
  ex <- colSums(E * cwn)
  bt <- vapply(seq_len(200), function(b) { i <- sample.int(nrow(E), replace = TRUE); colSums(E[i, , drop = FALSE] * cwn[i]) / sum(cwn[i]) }, numeric(k))
  bt <- matrix(bt, nrow = k)
  # full-data pooled directions for display
  Cfull <- Reduce(`+`, Map(function(m, c) m$C * c, M, cwn)); ev <- eigen(Cfull, symmetric = TRUE)
  V <- ev$vectors[, seq_len(k), drop = FALSE]; sg <- apply(V, 2, function(v) sign(v[which.max(abs(v))])); V <- sweep(V, 2, sg, "*")
  top <- vapply(seq_len(k), function(j) { o <- order(-abs(V[, j]))[seq_len(min(5, J))]
    paste(sprintf("%s (%+.2f)", genes[o], V[o, j]), collapse = ", ") }, "")
  out <- data.frame(component = paste0("PC", seq_len(k)), share = (ev$values / sum(ev$values))[seq_len(k)], p = p, top_genes = top,
                    share_heldout = obs$share, null_heldout = mn, excess = ex,
                    excess_lower = apply(bt, 1, stats::quantile, 0.025), excess_upper = apply(bt, 1, stats::quantile, 0.975),
                    row.names = NULL)
  call <- out$p <= alpha
  if (!is.null(control)) {
    # group statistic: mean held-out excess of the target windows along the group's directions minus the mean
    # excess of the control windows along the same directions; null from pseudo-target groups of controls
    stat_of <- function(Eg, dirs, idx_ctrl) colMeans(Eg) - colMeans(.rff_ctrl_excess(control, dirs, idx_ctrl))
    Tobs <- stat_of(E, (obs$VA + obs$VB) / 2, seq_along(control$C))
    nc <- length(control$M); m <- min(length(fits), nc - 3)
    if (m < 1) stop("Too few control windows for the group null.")
    Tnull <- t(vapply(seq_len(n_group_null), function(b) {
      g <- sample.int(nc, m); Mg <- control$M[g]
      og <- .rff_cf_group_share(Mg, rep(1 / m, m), k)
      Eg <- t(vapply(Mg, .rff_cf_window_excess, numeric(k), VA = og$VA, VB = og$VB)); if (k == 1) Eg <- matrix(Eg, ncol = 1)
      stat_of(Eg, (og$VA + og$VB) / 2, setdiff(seq_len(nc), g))
    }, numeric(k))); if (k == 1) Tnull <- matrix(Tnull, ncol = 1)
    out$excess_vs_control <- Tobs
    out$p_control <- vapply(seq_len(k), function(j) (1 + sum(Tnull[, j] >= Tobs[j])) / (n_group_null + 1), 0)
    call <- call & out$p_control <= alpha
    attr(out, "control_null") <- Tnull
  }
  if (sequential) call <- cumprod(call) == 1
  out$call <- call
  dimnames(V) <- list(genes, out$component)
  attr(out, "loadings") <- V; attr(out, "window_excess") <- E
  attr(out, "scores") <- NULL
  attr(out, "null_shares") <- lapply(seq_len(k), function(j) nul[, j]); attr(out, "null") <- "shift_joint_crossfit"
  attr(out, "weights") <- stats::setNames(cw, wn)
  out
}

# ---- from length scale to reach ----------------------------------------------

# RBF length scale matching (at half correlation) the bin-averaged field of
# point sources with kernel exp(-d / lambda) (shot noise; its correlation is
# the self-convolution of the kernel).
.rff_kernel_ell <- function(lambda, bin) {
  L <- 12 * lambda + 4 * bin
  s <- max(min(lambda, bin) / 4, 2 * L / 1024); n <- 2^ceiling(log2(2 * L / s))
  x <- (0:(n - 1)) * s; x <- ifelse(x > n * s / 2, x - n * s, x)
  k <- exp(-sqrt(outer(x^2, x^2, "+")) / lambda)
  box <- outer(abs(x) < bin / 2, abs(x) < bin / 2) * 1
  if (sum(box) == 0) box[1, 1] <- 1
  A <- Re(stats::fft(abs(stats::fft(k) * stats::fft(box))^2, inverse = TRUE))
  pr <- A[1, seq_len(n / 2)] / A[1, 1]; xs <- x[seq_len(n / 2)]
  i <- which(pr < 0.5)[1]
  r12 <- stats::approx(pr[(i - 1):i], xs[(i - 1):i], 0.5)$y
  r12 / sqrt(2 * log(2))
}

#' Reach (decay length) implied by a profiled length scale
#'
#' **Experimental.** The length scale `ell` of [rff_lengthscale_profile()]
#' describes a Gaussian (RBF) correlation. If a program is a response to
#' point sources (producer cells) that decays as \eqn{e^{-d/\lambda}}, the
#' response map is shot noise whose correlation is the self-convolution of
#' that kernel, which is matched at half correlation by an RBF length scale
#' of about \eqn{1.72\lambda} (more for `lambda` near the bin size, where
#' bin averaging adds to the width). `rff_reach()` inverts this relation.
#' The result is a reach under that model only, and it is only identified
#' when the length scale is well inside the window (see
#' `rff_lengthscale_profile()`'s `identifiable` column).
#'
#' @param ell Length scale(s) (coordinate units).
#' @param bin_size Bin size used for the fit.
#' @param kernel `"exponential"` (decay \eqn{e^{-d/\lambda}}) or `"gaussian"`
#'   (decay \eqn{e^{-d^2/2\lambda^2}}: its shot noise has RBF length scale
#'   \eqn{\sqrt{2}\lambda}, bin averaging ignored).
#' @return Numeric vector of reaches (`NA` where `ell` is too small to be
#'   produced by any reach at this bin size).
#' @export
#' @examples
#' rff_reach(c(20, 70, 140, 550), bin_size = 8)
rff_reach <- function(ell, bin_size, kernel = c("exponential", "gaussian")) {
  kernel <- match.arg(kernel)
  if (kernel == "gaussian") return(ell / sqrt(2))
  lo <- .rff_kernel_ell(bin_size / 50, bin_size)
  vapply(ell, function(e) {
    if (!is.finite(e) || e <= lo) return(NA_real_)
    f <- function(ll) .rff_kernel_ell(exp(ll), bin_size) - e
    up <- log(max(e, bin_size))
    exp(stats::uniroot(f, c(log(bin_size / 50), up), tol = 1e-3)$root)
  }, 0)
}
