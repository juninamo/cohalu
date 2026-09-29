# Length-scale profile, cross-fitted / control-calibrated program test,
# transfer test and the binomial (fraction) model.

cal_grid <- function(Y, res) {
  n <- as.integer(round(sqrt(nrow(Y))))
  gx <- rep(seq_len(n) - 1, times = n); gy <- rep(seq_len(n) - 1, each = n)
  structure(list(counts = Matrix::Matrix(Y, sparse = TRUE),
                 coords = data.frame(x = (gx + 0.5) * res, y = (gy + 0.5) * res, in_tissue = TRUE),
                 grid = list(nx = n, ny = n, bin_size = res, xmin = 0, ymin = 0), genes = colnames(Y)),
            class = "binned_transcripts")
}
cal_sim <- function(seed, ell = 10, amp = 0.8, n = 30, res = 5, J = 12, prog = 1:5, own = 0.2) {
  set.seed(seed)
  f1 <- .grf_fft(n, n, res, ell)
  eta <- log(0.3) + outer(as.vector(f1), as.numeric(seq_len(J) %in% prog) * amp) +
    vapply(seq_len(J), function(j) own * as.vector(.grf_fft(n, n, res, 15)), numeric(n * n))
  Y <- matrix(stats::rpois(length(eta), exp(eta) * res^2), n * n, J, dimnames = list(NULL, paste0("g", seq_len(J))))
  cal_grid(Y, res)
}

test_that("rff_lengthscale_profile() ranks length scales and returns curves", {
  skip_on_cran()
  est <- vapply(c(6, 24), function(ell) {
    b <- cal_sim(2, ell = ell, n = 36, own = 0.1)
    fit <- fit_spatial_rff(b, n_factors = 1, lengthscales = 8, basis = "grid", learn_lengthscales = FALSE,
                           factor_init = "residual_pca", max_iter = 60)
    pr <- rff_lengthscale_profile(fit, b, ls_grid = c(3, 6, 12, 24, 48), folds = 3, holdout_size = 5, max_iter = 60, n_boot = 30)
    expect_s3_class(pr, "rff_ls_profile")
    expect_equal(nrow(pr$curves), 5)
    expect_true(all(c("lengthscale", "lower", "upper", "gain", "heldout_dev_explained", "at_boundary") %in% names(pr$summary)))
    expect_true(pr$summary$lower <= pr$summary$upper)
    expect_gt(pr$summary$gain, 0)
    expect_output(print(pr), "length scale")
    pr$summary$lengthscale
  }, 0)
  expect_lt(est[1], est[2])
  expect_error(rff_lengthscale_profile(fit_spatial_rff(cal_sim(1), n_factors = 1, max_iter = 5), cal_sim(1), ls_grid = 5), "ls_grid")
})

test_that("fit_spatial_rff(lengthscales = \"profile\") refits with profiled length scales", {
  skip_on_cran()
  b <- cal_sim(3, ell = 16)
  f <- fit_spatial_rff(b, n_factors = 2, lengthscales = "profile", ls_grid = c(4, 8, 16, 32), basis = "grid",
                       factor_init = "residual_pca", ard = 5, max_iter = 50,
                       profile = list(init = 8, folds = 3, max_iter = 50, n_boot = 20, factors = 1))
  expect_s3_class(f$lengthscale_profile, "rff_ls_profile")
  expect_equal(unname(f$lengthscales[1]), f$lengthscale_profile$summary$lengthscale, tolerance = 1e-8)
  expect_equal(unname(f$lengthscales[2]), 8)
  expect_identical(f$settings$lengthscale_method, "profile")
  expect_false(f$settings$learn_lengthscales)
  expect_error(fit_spatial_rff(b, lengthscales = "auto"), "profile")
  # start = earlier fit keeps the factor identities
  g <- fit_spatial_rff(b, n_factors = 2, lengthscales = c(8, 8), basis = "grid", learn_lengthscales = FALSE, start = f, max_iter = 5)
  expect_gt(abs(cor(g$L[, 1], f$L[, 1])), 0.9)
  expect_error(fit_spatial_rff(b, n_factors = 3, basis = "grid", start = f, max_iter = 2), "start")
})

test_that("binomial family: gradient, fraction data, shift and factor tests", {
  set.seed(4); n <- 18; res <- 6; J <- 8
  Ntr <- matrix(stats::rpois(n * n * J, 4), n * n, J, dimnames = list(NULL, paste0("g", seq_len(J))))
  pij <- stats::plogis(stats::rnorm(J) + outer(0.8 * as.vector(.grf_fft(n, n, res, 15)), rep(1, J)))
  pij <- matrix(pij, n * n, J)
  Y <- matrix(stats::rbinom(length(Ntr), Ntr, pij), n * n, J, dimnames = dimnames(Ntr))
  b <- cal_grid(Y, res)
  expect_error(fit_spatial_rff(b, family = "binomial"), "trials")
  expect_error(fit_spatial_rff(b, family = "binomial", trials = Ntr - 5), "trials")
  f <- fit_spatial_rff(b, n_factors = 2, family = "binomial", trials = Ntr, basis = "grid", lengthscales = 10,
                       learn_lengthscales = FALSE, ard = 1, max_iter = 40)
  expect_identical(f$family, "binomial"); expect_null(f$dispersion); expect_equal(dim(f$trials), dim(Y))
  old <- options(cohalu.rff_debug = TRUE); on.exit(options(old))
  invisible(fit_spatial_rff(b, n_factors = 2, family = "binomial", trials = Ntr, basis = "grid", lengthscales = 10, max_iter = 1))
  d <- get(".rff_dbg", envir = globalenv()); rm(".rff_dbg", envir = globalenv())
  set.seed(2); th <- d$theta0 + stats::rnorm(length(d$theta0), sd = 0.05)
  ch <- c(d$idx$gamma[c(2, 50)], d$idx$L[1:3], d$idx$alpha[1:2], d$idx$log_sigma0)
  num <- vapply(ch, function(i) { a <- th; a[i] <- a[i] + 1e-5; z <- th; z[i] <- z[i] - 1e-5
    (d$objgrad(a)$value - d$objgrad(z)$value) / 2e-5 }, 0)
  expect_equal(d$objgrad(th)$grad[ch], num, tolerance = 1e-4)
  # the objective is the binomial negative log-likelihood (plus priors)
  pt <- rff_program_test(f, b, null = "shift", n_boot = 9, n_components = 2)
  expect_equal(nrow(pt), 2)
  expect_error(rff_program_test(f, b, null = "parametric", n_boot = 3), "binomial")
  ft <- rff_factor_test(f, b, n_boot = 2)
  expect_equal(nrow(ft), 2)
})

test_that("cross-fitted program test, control reference and effect-size threshold", {
  skip_on_cran()
  mk <- function(seed, amp) {
    b <- cal_sim(seed, ell = 8, amp = amp, n = 36, J = 14, prog = 1:5)
    list(b = b, fit = fit_spatial_rff(b, n_factors = 3, lengthscales = 8, basis = "grid", learn_lengthscales = FALSE,
                                      factor_init = "residual_pca", ard = 20, max_iter = 40))
  }
  ctrl <- lapply(11:14, mk, amp = 0)
  ref <- rff_control_reference(lapply(ctrl, `[[`, "fit"), lapply(ctrl, `[[`, "b"), n_boot = 9, n_components = 3)
  expect_s3_class(ref, "rff_control")
  expect_equal(dim(ref$excess), c(4, 3))
  expect_true(all(c("p_heldout", "p_control", "called") %in% names(ref$loo)))
  expect_output(print(ref), "negative-control")
  w <- mk(21, amp = 1)
  cf <- rff_program_test(w$fit, w$b, null = "shift", n_boot = 19, n_components = 3, crossfit = TRUE, n_boot_ci = 50)
  expect_true(all(c("share_heldout", "null_heldout", "excess", "excess_lower", "excess_upper", "call") %in% names(cf)))
  expect_true(all(cf$excess_lower <= cf$excess_upper))
  expect_true(cf$call[1])
  expect_gt(cf$excess[1], 0.5)
  cc <- rff_program_test(w$fit, w$b, null = "shift", n_boot = 9, n_components = 3, control = ref, n_boot_ci = 50)
  expect_true(all(c("p_control", "control_q95", "p_control_dir") %in% names(cc)))
  expect_equal(cc$p_control[1], 1 / 5)
  me <- rff_program_test(w$fit, w$b, null = "shift", n_boot = 9, n_components = 3, min_effect = 100, n_boot_ci = 50)
  expect_false(any(me$call))
  expect_error(rff_program_test(w$fit, w$b, crossfit = TRUE, n_boot = 3), "shift")
  # the robust option of the plain shift test runs; the default is unchanged
  p0 <- rff_program_test(w$fit, w$b, null = "shift", n_boot = 9, n_components = 2)
  p1 <- rff_program_test(w$fit, w$b, null = "shift", n_boot = 9, n_components = 2, robust = TRUE)
  expect_equal(nrow(p1), 2); expect_null(p0$call)
  # rff_programs uses `call` for significance and reports effect sizes
  r <- rff_programs(w$b, n_factors = 3, lengthscale = 8, ard = 20, n_boot = 19, crossfit = TRUE, fit_args = list(max_iter = 40))
  expect_true(all(c("excess", "excess_lower", "p_control", "lengthscale_lower") %in% names(r$programs)))
  expect_true(any(r$programs$status == "Confirmed"))
  f <- tempfile(fileext = ".html")
  rff_report(r, w$b, file = f, open = FALSE, export = FALSE)
  expect_true(grepl("Effect sizes (cross-fitted program test)", paste(readLines(f), collapse = "\n"), fixed = TRUE))
})

test_that("rff_transfer_test(): true program vs permuted and irrelevant loadings", {
  skip_on_cran()
  src <- cal_sim(31, ell = 8, amp = 1, n = 36, J = 14, prog = 1:5)
  tgt <- cal_sim(32, ell = 8, amp = 1, n = 36, J = 14, prog = 1:5)
  fs <- fit_spatial_rff(src, n_factors = 2, lengthscales = 8, basis = "grid", learn_lengthscales = FALSE,
                        factor_init = "residual_pca", ard = 20, max_iter = 40)
  ft <- fit_spatial_rff(tgt, n_factors = 1, lengthscales = 8, basis = "grid", learn_lengthscales = FALSE, max_iter = 20)
  k <- which.max(colSums(fs$L[1:5, , drop = FALSE]^2))
  set.seed(1)
  L <- cbind(true = fs$L[, k], irrelevant = stats::setNames(c(rep(0, 7), stats::rnorm(7)), fs$genes))
  rownames(L) <- fs$genes
  tt <- rff_transfer_test(L, ft, tgt, n_perm = 199, strata = 2)
  expect_equal(tt$program, c("true", "irrelevant"))
  expect_lt(tt$p[1], 0.05)
  expect_gt(tt$p[2], 0.05)
  expect_equal(ncol(attr(tt, "null")), 199)
  th <- rff_transfer_test(L[, 1], ft, tgt, statistic = "heldout", n_perm = 4, folds = 2, max_iter = 20)
  expect_true(is.finite(th$statistic))
  expect_error(rff_transfer_test(matrix(1, 2, 1), ft, tgt), "row names")
})

test_that("rff_programs(lengthscale = \"profile\") reports profiled length scales and the report shows them", {
  skip_on_cran()
  b <- cal_sim(5, ell = 10, n = 36, J = 14, prog = 1:5)
  r <- rff_programs(b, n_factors = 3, lengthscale = "profile", ls_grid = c(5, 10, 20, 40), ard = 20, n_boot = 19, crossfit = TRUE,
                    fit_args = list(max_iter = 50), profile = list(init = 10, folds = 3, max_iter = 50, n_boot = 30))
  p <- r$programs[r$programs$status == "Confirmed", ]
  expect_gte(nrow(p), 1)
  expect_true(all(is.finite(p$lengthscale_lower)))
  expect_true(p$lengthscale_lower[1] <= p$lengthscale[1] && p$lengthscale[1] <= p$lengthscale_upper[1])
  f <- tempfile(fileext = ".html")
  rff_report(r, b, file = f, open = FALSE, export = FALSE)
  expect_true(grepl("Length-scale profiles", paste(readLines(f), collapse = "\n"), fixed = TRUE))
})

test_that("rff_reach() inverts the exponential-kernel length scale", {
  ell <- vapply(c(10, 40, 160), function(l) .rff_kernel_ell(l, 8), 0)
  expect_equal(ell / c(10, 40, 160), c(1.755, 1.724, 1.722), tolerance = 0.01)
  expect_equal(rff_reach(ell, 8), c(10, 40, 160), tolerance = 0.01)
  expect_equal(rff_reach(10, 8, kernel = "gaussian"), 10 / sqrt(2))
  expect_true(is.na(rff_reach(1, 8)))
})

test_that("variational (ELBO) length scale and the MCMC reference", {
  skip_on_cran()
  b <- cal_sim(3, ell = 12, n = 30, own = 0.05)
  f <- fit_spatial_rff(b, n_factors = 1, lengthscales = 8, basis = "grid", learn_lengthscales = FALSE,
                       factor_init = "residual_pca", max_iter = 60)
  v <- rff_lengthscale_vi(f, b, ls_grid = c(3, 6, 12, 24, 48))
  expect_s3_class(v, "rff_ls_profile")
  expect_true(all(c("lengthscale", "amplitude", "reach", "identifiable") %in% names(v$summary)))
  expect_gt(v$summary$lengthscale, 6); expect_lt(v$summary$lengthscale, 24)
  expect_gt(v$summary$gain, 0)
  # numerical pieces: Lanczos log-determinant and conjugate gradients against dense algebra
  bs <- cal_sim(2, n = 12); Kop <- .vi_kernel(bs$grid, bs$coords$in_tissue, 10); N <- Kop$N
  K <- Kop$mv(diag(N)); d <- seq(0.2, 1.5, length.out = N); B <- diag(N) + diag(d) %*% K %*% diag(d)
  set.seed(1); Z <- matrix(sample(c(-1, 1), N * 40, TRUE), N)
  expect_equal(.vi_logdet(Kop$mv, d, Z), as.numeric(determinant(B)$modulus), tolerance = 0.1)
  x <- .vi_cg(Kop$mv, d, matrix(seq_len(N), ncol = 1)); expect_lt(max(abs(B %*% x - seq_len(N))), 1e-3)
  m <- rff_lengthscale_mcmc(f, b, n_iter = 300, burn = 100, ls_range = c(2, 100))
  expect_equal(nrow(m$samples), 200)
  expect_true(m$summary$lower <= m$summary$lengthscale && m$summary$lengthscale <= m$summary$upper)
})

test_that("multinomial family: gradient, fields and shift test", {
  b <- cal_sim(4, ell = 8, amp = 1, n = 24, J = 10, prog = 1:4)
  f <- fit_spatial_rff(b, n_factors = 2, family = "multinomial", basis = "grid", lengthscales = 8, learn_lengthscales = FALSE,
                       ard = 5, max_iter = 40, factor_init = "residual_pca")
  expect_false(f$has_density); expect_equal(length(f$bin_total), sum(b$coords$in_tissue))
  old <- options(cohalu.rff_debug = TRUE); on.exit(options(old))
  invisible(fit_spatial_rff(b, n_factors = 2, family = "multinomial", basis = "grid", lengthscales = 8, max_iter = 1))
  d <- get(".rff_dbg", envir = globalenv()); rm(".rff_dbg", envir = globalenv())
  set.seed(2); th <- d$theta0 + stats::rnorm(length(d$theta0), sd = 0.05)
  ch <- c(d$idx$gamma[c(2, 50)], d$idx$L[1:3], d$idx$alpha[1:2])
  num <- vapply(ch, function(i) { a <- th; a[i] <- a[i] + 1e-5; z <- th; z[i] <- z[i] - 1e-5
    (d$objgrad(a)$value - d$objgrad(z)$value) / 2e-5 }, 0)
  expect_equal(d$objgrad(th)$grad[ch], num, tolerance = 1e-4)
  # a bin-wide multiplicative factor leaves the multinomial fit unchanged
  b2 <- b; b2$counts <- b$counts * 2
  f2 <- fit_spatial_rff(b2, n_factors = 2, family = "multinomial", basis = "grid", lengthscales = 8, learn_lengthscales = FALSE,
                        ard = 5, max_iter = 40, factor_init = "residual_pca")
  expect_equal(unname(f2$alpha - mean(f2$alpha)), unname(f$alpha - mean(f$alpha)), tolerance = 0.05)
  pt <- rff_program_test(f, b, null = "shift", n_boot = 9, n_components = 2)
  expect_equal(nrow(pt), 2)
  expect_error(rff_factor_test(f, b, n_boot = 2), "multinomial")
})

test_that("rff_program_test_joint(): cross-fitted and control-calibrated group test", {
  skip_on_cran()
  mk <- function(seed, amp) { b <- cal_sim(seed, ell = 8, amp = amp, n = 30, J = 12, prog = 1:5)
    list(b = b, fit = fit_spatial_rff(b, n_factors = 2, lengthscales = 8, basis = "grid", learn_lengthscales = FALSE,
                                      factor_init = "residual_pca", ard = 20, max_iter = 30)) }
  ctrl <- lapply(11:16, mk, amp = 0)
  ref <- rff_control_reference(lapply(ctrl, `[[`, "fit"), lapply(ctrl, `[[`, "b"), n_boot = 9, n_components = 2)
  expect_length(ref$M, 6)
  tg <- lapply(21:23, mk, amp = 0.8)
  jj <- rff_program_test_joint(lapply(tg, `[[`, "fit"), lapply(tg, `[[`, "b"), n_boot = 19, n_components = 2, control = ref, n_group_null = 99)
  expect_true(all(c("share_heldout", "excess", "excess_vs_control", "p_control", "call") %in% names(jj)))
  expect_true(jj$call[1])
  expect_equal(dim(attr(jj, "window_excess")), c(3, 2))
  expect_error(rff_program_test_joint(lapply(tg, `[[`, "fit"), lapply(tg, `[[`, "b"), loadings = diag(12)[, 1, drop = FALSE], crossfit = TRUE), "discovery")
})
