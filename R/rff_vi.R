# Variational (ELBO) length-scale selection for one spatial program with fixed
# loadings (Gaussian variational posterior of the field, Poisson likelihood with
# closed-form expectations; Opper & Archambeau 2009 parametrisation). The field
# prior is the grid-basis RBF Gaussian process of fit_spatial_rff(basis =
# "grid"), applied by FFT on the padded torus; all matrix functions (solves,
# log-determinants, traces, diagonals) are computed by conjugate gradients,
# Lanczos quadrature and Hutchinson probes with FFT matrix-vector products.

# covariance operator of the unit-variance RBF field (length scale ell, prior
# variance s2) restricted to the in-tissue bins
.vi_kernel <- function(grid, keep, ell, s2 = 1, nugget = 0) {
  nx <- grid$nx; ny <- grid$ny; bs <- grid$bin_size
  padn <- function(n) stats::nextn(n + min(n, ceiling(2 * ell / bs) + 2))
  px <- padn(nx); py <- padn(ny); P <- px * py
  fq <- function(p) { k <- 0:(p - 1); ifelse(k < p / 2, k, k - p) * 2 * pi / (p * bs) }
  S <- exp(-0.5 * ell^2 * outer(fq(px)^2, fq(py)^2, "+")); ev <- s2 * P * S / sum(S)   # eigenvalues on the torus
  idx <- (which(keep) - 1L) %% nx + 1L + ((which(keep) - 1L) %/% nx) * px
  mv <- function(V) {                                  # K %*% V for a matrix of columns
    V <- as.matrix(V)
    apply(V, 2, function(v) { Z <- matrix(0, px, py); Z[idx] <- v; Re(stats::fft(stats::fft(Z) * ev, inverse = TRUE))[idx] / P + nugget * v })
  }
  list(mv = mv, N = length(idx))
}

# conjugate gradients for B x = b with B = I + D K D (D = diag(d)), several right-hand sides
.vi_cg <- function(Kmv, d, Bm, tol = 1e-6, maxit = 200) {
  Bmv <- function(X) X + d * Kmv(d * X)
  X <- matrix(0, nrow(Bm), ncol(Bm)); R <- Bm; P <- R; rs <- colSums(R^2)
  for (it in seq_len(maxit)) {
    AP <- Bmv(P); al <- rs / pmax(colSums(P * AP), 1e-300)
    X <- X + sweep(P, 2, al, "*"); R <- R - sweep(AP, 2, al, "*")
    rn <- colSums(R^2); if (any(!is.finite(rn))) stop("CG diverged")
    if (all(sqrt(rn) <= tol * sqrt(colSums(Bm^2)) + 1e-12)) break
    P <- R + sweep(P, 2, rn / pmax(rs, 1e-300), "*"); rs <- rn
  }
  X
}

# log det(I + D K D) by stochastic Lanczos quadrature
.vi_logdet <- function(Kmv, d, probes, steps = 25) {
  Bmv <- function(x) x + d * Kmv(d * x)
  est <- apply(probes, 2, function(z) {
    q <- z / sqrt(sum(z^2)); Q <- matrix(0, length(z), steps); a <- b <- numeric(steps); qp <- 0; bp <- 0
    for (j in seq_len(steps)) {
      Q[, j] <- q; w <- drop(Bmv(q)); a[j] <- sum(w * q); w <- w - a[j] * q - bp * qp
      w <- w - Q[, seq_len(j), drop = FALSE] %*% crossprod(Q[, seq_len(j), drop = FALSE], w)
      b[j] <- sqrt(sum(w^2)); if (b[j] < 1e-10) { steps <- j; break }
      qp <- q; bp <- b[j]; q <- drop(w) / b[j]
    }
    Tm <- diag(a[seq_len(steps)], steps); if (steps > 1) { Tm[cbind(2:steps, 1:(steps - 1))] <- b[seq_len(steps - 1)]; Tm[cbind(1:(steps - 1), 2:steps)] <- b[seq_len(steps - 1)] }
    e <- eigen(Tm, symmetric = TRUE); sum(e$vectors[1, ]^2 * log(pmax(e$values, 1e-12))) * sum(z^2)
  })
  mean(est)
}

# Gaussian variational inference of the field for fixed loadings l (genes),
# offset eta0 (N x G), prior K (unit variance x s2); returns the ELBO and q(f)
.vi_field <- function(Y, eta0, l, Kop, n_probe = 8, iters = 12, seed = 1, lfy = NULL) {
  N <- nrow(Y); set.seed(seed)
  Zp <- matrix(sample(c(-1, 1), N * n_probe, replace = TRUE), N, n_probe)
  m <- numeric(N); v <- rep(0, N)
  if (is.null(lfy)) lfy <- sum(lgamma(Y + 1))
  lam <- function(m, v) exp(pmin(eta0 + outer(m, l) + outer(v, l^2 / 2), 30))
  for (it in seq_len(iters)) {
    Mu <- lam(m, v)
    om <- drop(Mu %*% l^2); g <- drop((Y - Mu) %*% l)
    d <- sqrt(pmax(om, 1e-10)); cvec <- om * m + g
    # Sigma c = K c - K D B^-1 D K c
    Kc <- drop(Kop$mv(cvec)); Kz <- Kop$mv(Zp)
    sol <- .vi_cg(Kop$mv, d, cbind(d * Kc, d * Kz))
    m_new <- Kc - drop(Kop$mv(d * sol[, 1]))
    SZ <- Kz - Kop$mv(d * sol[, -1, drop = FALSE])
    v <- pmax(rowMeans(Zp * SZ), 1e-10)
    step <- max(abs(m_new - m)); m <- 0.5 * m + 0.5 * m_new
    if (step < 1e-3 && it > 3) break
  }
  Mu <- lam(m, v); om <- drop(Mu %*% l^2); d <- sqrt(pmax(om, 1e-10))
  Ell <- sum(Y * (eta0 + outer(m, l))) - sum(Mu) - lfy
  # m' K^-1 m with m = Sigma c:  K^-1 Sigma c = c - D B^-1 D K c
  cvec <- om * m + drop((Y - Mu) %*% l); Kc <- drop(Kop$mv(cvec))
  sol <- .vi_cg(Kop$mv, d, cbind(d * Kc, Zp))
  Kim <- cvec - d * sol[, 1]; mKm <- sum(m * Kim)
  trBinv <- mean(colSums(Zp * sol[, -1, drop = FALSE]))
  ld <- .vi_logdet(Kop$mv, d, Zp)
  KL <- 0.5 * (trBinv - N + ld + mKm)
  list(elbo = Ell - KL, ell = Ell, kl = KL, m = m, v = v)
}

#' Length scale of each program by the variational evidence (ELBO)
#'
#' **Experimental.** Alternative to the held-out criterion of
#' [rff_lengthscale_profile()]: for each program (factor) of `fit`, with its
#' loading direction fixed and everything else in the fit as offset, the
#' program field gets a Gaussian variational posterior under the grid-basis
#' RBF prior and a Poisson likelihood (closed-form expectations;
#' Opper-Archambeau parametrisation, conjugate gradients, Lanczos
#' log-determinants and Hutchinson traces with FFT products), and the ELBO -
#' a lower bound on the marginal likelihood, which, unlike the MAP
#' objective, charges for the flexibility of short length scales - is
#' maximised over the length scale grid and the field amplitude (prior
#' variance, golden-section search in log scale). Only the program genes
#' (share `gene_share` of the squared loadings) enter.
#'
#' @inheritParams rff_lengthscale_profile
#' @param amp_range Range of the field amplitude (prior standard deviation
#'   times the loading norm) searched at each length scale.
#' @param n_probe Hutchinson / Lanczos probe vectors.
#' @param nugget If `TRUE` (default), the program field is a smooth RBF
#'   part plus an independent bin-level part (relative variance estimated
#'   with the amplitude by maximising the ELBO). Without it, on real tissue
#'   the ELBO chose the smallest length scale: responses are carried by
#'   single cells, so the dominant correlation scale of a response map is
#'   the cell, not the reach.
#' @return An object of class `rff_ls_profile` (as
#'   [rff_lengthscale_profile()], `criterion = "elbo"`): `summary` with
#'   `lengthscale` (maximum of the ELBO, parabola in log length scale),
#'   `amplitude`, `elbo_gain` (ELBO at the best length scale minus the ELBO
#'   without the program), reach and identifiability; `curves` with the
#'   ELBO per length scale (`gain` = ELBO minus no-program ELBO).
#' @seealso [rff_lengthscale_profile()]
#' @export
rff_lengthscale_vi <- function(fit, binned, factors = NULL, ls_grid = 5 * 2^(0:6), gene_share = 0.95,
                               amp_range = c(0.05, 5), n_probe = 8, seed = 1, n_cores = 1, reach_kernel = "exponential",
                               nugget = TRUE) {
  if (!inherits(fit, "spatial_rff_fit")) stop("`fit` must come from fit_spatial_rff().")
  if (!(fit$family %in% c("nb", "poisson", "multinomial"))) stop("rff_lengthscale_vi() supports count families (nb / poisson / multinomial).")
  ls_grid <- sort(unique(as.numeric(ls_grid)))
  keep <- fit$in_tissue; g <- fit$grid; L <- fit$L; genes <- fit$genes
  pstr <- sqrt(colSums(sweep(L, 2, colMeans(L))^2)); uni <- 1 - pstr^2 / pmax(colSums(L^2), 1e-12)
  if (is.null(factors)) factors <- colnames(L)[pstr >= 0.1 * max(pstr) & uni <= 0.5]
  if (is.numeric(factors)) factors <- colnames(L)[factors]
  Yall <- as.matrix(binned$counts[keep, genes, drop = FALSE])
  cxy <- binned$coords[keep, c("x", "y")]; extent <- min(diff(range(cxy$x)), diff(range(cxy$y))) + g$bin_size
  summ <- list(); curves <- list()
  for (k in factors) {
    o <- order(-L[, k]^2); cs <- cumsum(L[o, k]^2) / sum(L[, k]^2)
    gk <- genes[o[seq_len(max(min(5, length(o)), which(cs >= gene_share)[1]))]]
    eta0 <- .rff_linpred(fit, setdiff(colnames(L), k), genes = gk)
    Y <- Yall[, gk, drop = FALSE]; u <- L[gk, k] / sqrt(sum(L[gk, k]^2)); lfy <- sum(lgamma(Y + 1))
    base <- sum(Y * eta0) - sum(exp(eta0)) - lfy                          # no program (ELBO = log-likelihood)
    one <- function(ell) {
      if (!nugget) {
        Kop <- .vi_kernel(g, keep, ell, 1)
        f <- function(la) { e <- tryCatch(.vi_field(Y, eta0, exp(la) * u, Kop, n_probe = n_probe, seed = seed, lfy = lfy)$elbo, error = function(err) NA_real_)
          if (!is.finite(e)) 1e15 else -e }
        op <- stats::optimize(f, log(amp_range), tol = 0.05)
        return(c(elbo = -op$objective, amp = exp(op$minimum), nugget = 0))
      }
      # field = smooth RBF part + independent bin-level part (nugget, relative variance exp(lt)); the program
      # direction and amplitude are shared, so bin-scale granularity (single cells) does not drive the length scale
      f2 <- function(par) { Kop <- .vi_kernel(g, keep, ell, 1, nugget = exp(par[2]))
        e <- tryCatch(.vi_field(Y, eta0, exp(par[1]) * u, Kop, n_probe = n_probe, seed = seed, lfy = lfy)$elbo, error = function(err) NA_real_)
        if (!is.finite(e)) 1e15 else -e }
      op <- stats::optim(c(log(0.5), log(0.5)), f2, method = "Nelder-Mead", control = list(maxit = 40, reltol = 1e-6))
      c(elbo = -op$value, amp = exp(op$par[1]), nugget = exp(op$par[2]))
    }
    res <- .rff_mclapply(as.list(ls_grid), one, n_cores)
    E <- vapply(res, `[[`, 0, "elbo"); A <- vapply(res, `[[`, 0, "amp"); Nu <- vapply(res, `[[`, 0, "nugget")
    ib <- which.max(E); est <- exp(.rff_parabola_max(log(ls_grid), E, ib))
    curves[[k]] <- data.frame(factor = k, lengthscale = ls_grid, gain = E - base, amplitude = A, nugget = Nu, row.names = NULL)
    summ[[k]] <- data.frame(factor = k, lengthscale = est, lower = NA_real_, upper = NA_real_, grid_best = ls_grid[ib],
                            gain = E[ib] - base, gain_se = NA_real_, heldout_dev_explained = NA_real_,
                            at_boundary = ib %in% c(1, length(ls_grid)), n_genes = length(gk),
                            top_genes = paste(sprintf("%s (%+.2f)", gk[1:min(5, length(gk))], L[gk[1:min(5, length(gk))], k]), collapse = ", "),
                            reach = NA_real_, reach_lower = NA_real_, reach_upper = NA_real_,
                            identifiable = est >= 2 * g$bin_size && est <= extent / 5 && !(ib %in% c(1, length(ls_grid))),
                            amplitude = A[ib], nugget = Nu[ib], row.names = NULL)
  }
  sm <- do.call(rbind, summ)
  if (!is.null(reach_kernel)) sm$reach <- rff_reach(sm$lengthscale, g$bin_size, reach_kernel)
  structure(list(summary = sm, curves = do.call(rbind, curves), boot = NULL,
                 settings = list(ls_grid = ls_grid, criterion = "elbo", folds = NA, holdout_size = NA, n_probe = n_probe, seed = seed)),
            class = "rff_ls_profile")
}

#' Posterior of a program's length scale by MCMC (reference)
#'
#' **Experimental; small windows.** Reference sampler for checking
#' [rff_lengthscale_profile()] and [rff_lengthscale_vi()]: for one program
#' of `fit` with its loading direction fixed and everything else in the fit
#' as offset, the whitened field weights of the grid-basis RBF prior are
#' updated by elliptical slice sampling and the log length scale and log
#' amplitude by random-walk Metropolis (uniform priors on the log scale),
#' under a Poisson likelihood.
#'
#' @param fit,binned As in [rff_lengthscale_profile()].
#' @param factor The factor (name or index).
#' @param ls_range,amp_range Prior ranges (uniform on the log scale).
#' @param n_iter,burn Iterations and burn-in.
#' @param gene_share Share of the squared loadings carried by the genes used.
#' @param seed Random seed.
#' @return A list with `samples` (data frame of `lengthscale`, `amplitude`
#'   after burn-in), `summary` (posterior median and 95% interval of the
#'   length scale and amplitude, acceptance rates) and `field_mean`
#'   (posterior mean field at the in-tissue bins).
#' @export
rff_lengthscale_mcmc <- function(fit, binned, factor = 1, ls_range = c(2, 640), amp_range = c(0.02, 10),
                                 n_iter = 3000, burn = 1000, gene_share = 0.95, seed = 1) {
  if (!(fit$family %in% c("nb", "poisson", "multinomial"))) stop("Count families only.")
  .local_seed(seed)
  keep <- fit$in_tissue; g <- fit$grid; L <- fit$L; genes <- fit$genes
  k <- if (is.numeric(factor)) colnames(L)[factor] else factor
  o <- order(-L[, k]^2); cs <- cumsum(L[o, k]^2) / sum(L[, k]^2)
  gk <- genes[o[seq_len(max(min(5, length(o)), which(cs >= gene_share)[1]))]]
  eta0 <- .rff_linpred(fit, setdiff(colnames(L), k), genes = gk)
  Y <- as.matrix(binned$counts[keep, gk, drop = FALSE]); u <- L[gk, k] / sqrt(sum(L[gk, k]^2))
  nx <- g$nx; ny <- g$ny; bs <- g$bin_size
  padn <- function(n) stats::nextn(n + min(n, ceiling(2 * ls_range[2] / bs) + 2))
  px <- padn(nx); py <- padn(ny); P <- px * py
  fq <- function(p) { kk <- 0:(p - 1); ifelse(kk < p / 2, kk, kk - p) * 2 * pi / (p * bs) }
  k2 <- outer(fq(px)^2, fq(py)^2, "+")
  idx <- (which(keep) - 1L) %% nx + 1L + ((which(keep) - 1L) %/% nx) * px
  field <- function(z, ell) { S <- exp(-0.5 * ell^2 * k2); w <- sqrt(S / sum(S))
    Re(stats::fft(stats::fft(matrix(z, px, py)) * w, inverse = TRUE))[idx] / sqrt(P) }
  loglik <- function(f, a) { e <- eta0 + outer(a * f, u); sum(Y * e - exp(e)) }
  z <- stats::rnorm(P); lell <- log(20); la <- log(0.5)
  f <- field(z, exp(lell)); ll <- loglik(f, exp(la))
  keepS <- matrix(NA_real_, n_iter, 2); acc <- c(0, 0); fsum <- 0
  for (it in seq_len(n_iter)) {
    # elliptical slice sampling of z
    nu <- stats::rnorm(P); lly <- ll + log(stats::runif(1)); th <- stats::runif(1, 0, 2 * pi); lo <- th - 2 * pi; hi <- th
    repeat {
      zp <- z * cos(th) + nu * sin(th); fp <- field(zp, exp(lell)); llp <- loglik(fp, exp(la))
      if (llp > lly) { z <- zp; f <- fp; ll <- llp; break }
      if (th < 0) lo <- th else hi <- th
      th <- stats::runif(1, lo, hi)
    }
    # Metropolis on log length scale and log amplitude (whitened field kept)
    lp <- lell + stats::rnorm(1, 0, 0.15)
    if (lp > log(ls_range[1]) && lp < log(ls_range[2])) {
      fp <- field(z, exp(lp)); llp <- loglik(fp, exp(la))
      if (log(stats::runif(1)) < llp - ll) { lell <- lp; f <- fp; ll <- llp; acc[1] <- acc[1] + 1 }
    }
    lq <- la + stats::rnorm(1, 0, 0.1)
    if (lq > log(amp_range[1]) && lq < log(amp_range[2])) {
      llq <- loglik(f, exp(lq)); if (log(stats::runif(1)) < llq - ll) { la <- lq; ll <- llq; acc[2] <- acc[2] + 1 }
    }
    keepS[it, ] <- exp(c(lell, la)); if (it > burn) fsum <- fsum + exp(la) * f
  }
  sm <- as.data.frame(keepS[-seq_len(burn), , drop = FALSE]); names(sm) <- c("lengthscale", "amplitude")
  q <- function(v) stats::quantile(v, c(0.5, 0.025, 0.975))
  list(samples = sm,
       summary = data.frame(factor = k, lengthscale = q(sm$lengthscale)[1], lower = q(sm$lengthscale)[2], upper = q(sm$lengthscale)[3],
                            amplitude = q(sm$amplitude)[1], acc_lengthscale = acc[1] / n_iter, acc_amplitude = acc[2] / n_iter, row.names = NULL),
       field_mean = fsum / (n_iter - burn))
}
