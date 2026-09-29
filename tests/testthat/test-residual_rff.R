test_that("a zero matrix offset reproduces the area-only fit", {
  tx <- simulate_transcripts(size = 80, rate = 0.02, n_genes_per_set = 2, seed = 3)
  b <- bin_transcripts(tx, bin_size = 8)
  f1 <- fit_spatial_rff(b, n_factors = 2, n_features = 16, max_iter = 20)
  z <- matrix(0, sum(b$coords$in_tissue), length(b$genes), dimnames = list(NULL, b$genes))
  f2 <- fit_spatial_rff(b, n_factors = 2, n_features = 16, max_iter = 20, offset = z)
  expect_equal(f1$objective, f2$objective, tolerance = 1e-6)
  expect_identical(f2$offset, "matrix")
  expect_error(fit_spatial_rff(b, offset = z[-1, ]), "one row per bin")
})

test_that("rff_offset() absorbs known structure and ARD shrinks unneeded factors", {
  set.seed(4); n <- 30; res <- 6
  dom <- cohalu:::.grf_fft(n, n, res, 30) > 0
  J <- 12; prof <- matrix(rnorm(J * 2, 0, 1.2), J, 2)
  tx <- do.call(rbind, lapply(seq_len(J), function(j) {
    cnt <- rpois(n * n, exp(log(0.1) + ifelse(dom, prof[j, 1], prof[j, 2])) * res^2); idx <- rep(seq_len(n * n), cnt)
    data.frame(x = ((idx - 1) %% n + runif(length(idx))) * res, y = ((idx - 1) %/% n + runif(length(idx))) * res, gene = paste0("g", j))
  }))
  b <- bin_transcripts(tx, bin_size = res, tissue_radius = Inf)
  off <- rff_offset(b, data.frame(domain = factor(as.vector(dom))))
  expect_equal(dim(off), c(n * n, J))
  # the offset carries the domain contrast of each gene
  d1 <- colMeans(off[as.vector(dom), ]) - colMeans(off[!as.vector(dom), ])
  expect_gt(cor(d1, (prof[, 1] - prof[, 2])[match(colnames(off), paste0("g", seq_len(J)))]), 0.95)
  plain <- fit_spatial_rff(b, n_factors = 3, n_features = 32, max_iter = 60)
  resid <- fit_spatial_rff(b, n_factors = 3, n_features = 32, max_iter = 60, offset = off, ard = 30)
  expect_lt(max(resid$factor_strength), 0.5 * max(plain$factor_strength))   # nothing left to explain
  expect_length(resid$factor_strength, 3)
  expect_true(!is.null(resid$settings))
})

test_that("rff_factor_test() returns a p-value per factor", {
  tx <- simulate_transcripts(size = 70, rate = 0.02, n_genes_per_set = 2, seed = 5)
  b <- bin_transcripts(tx, bin_size = 8)
  f <- fit_spatial_rff(b, n_factors = 2, n_features = 16, max_iter = 20, ard = 2)
  expect_warning(rff_factor_test(f, b, n_boot = 2, max_iter = 10), "differs")
  tt <- rff_factor_test(f, b, n_boot = 2)
  expect_equal(nrow(tt), 2)
  expect_true(all(tt$p > 0 & tt$p <= 1))
  expect_length(attr(tt, "null_max"), 2)
  expect_true(all(c("program_strength", "uniform_share") %in% names(tt)))
  expect_true(all(tt$program_strength <= tt$strength + 1e-12))
  expect_equal(unname(f$program_strength), unname(sqrt(colSums(sweep(f$L, 2, colMeans(f$L))^2))))
  expect_length(attr(rff_factor_test(f, b, n_boot = 2, statistic = "strength"), "null_max"), 2)
})

test_that("the grid factorisation reproduces the direct random-feature evaluation", {
  tx <- simulate_transcripts(size = 70, rate = 0.02, n_genes_per_set = 2, seed = 6)
  b <- bin_transcripts(tx, bin_size = 7)
  z <- matrix(0.1, sum(b$coords$in_tissue), length(b$genes), dimnames = list(NULL, b$genes))
  f1 <- fit_spatial_rff(b, n_factors = 2, n_features = 16, max_iter = 15, offset = z, ard = 1)
  old <- options(cohalu.rff_grid = FALSE); on.exit(options(old))
  f2 <- fit_spatial_rff(b, n_factors = 2, n_features = 16, max_iter = 15, offset = z, ard = 1)
  expect_equal(f1$objective, f2$objective, tolerance = 1e-8)
  expect_equal(f1$L, f2$L, tolerance = 1e-6)
  expect_equal(f1$lengthscales, f2$lengthscales, tolerance = 1e-6)
})

test_that("rff_factor_test(): parallel refits; fit_spatial_rff(init =)", {
  tx <- simulate_transcripts(size = 70, rate = 0.02, n_genes_per_set = 2, seed = 5)
  b <- bin_transcripts(tx, bin_size = 8)
  f <- fit_spatial_rff(b, n_factors = 2, n_features = 16, max_iter = 15, ard = 2)
  t1 <- rff_factor_test(f, b, n_boot = 3)
  expect_true(is.numeric(attr(t1, "elapsed")))
  if (.Platform$OS.type == "unix") {
    t2 <- rff_factor_test(f, b, n_boot = 3, n_cores = 2)
    expect_equal(attr(t1, "null_max"), attr(t2, "null_max"))
  }
  # warm start needs the same random features
  expect_error(fit_spatial_rff(b, n_factors = 2, n_features = 16, max_iter = 5, seed = 2, init = f), "same")
  g <- fit_spatial_rff(b, n_factors = 2, n_features = 16, max_iter = 5, init = f)
  expect_s3_class(g, "spatial_rff_fit")
})

test_that("grid basis: analytic gradients, fields and settings", {
  set.seed(7); n <- 16; res <- 5
  dom <- cohalu:::.grf_fft(n, n, res, 20) > 0
  J <- 6; prof <- matrix(rnorm(J * 2), J, 2)
  lam <- exp(log(0.1) + outer(as.vector(dom), prof[, 1]) + outer(!as.vector(dom), prof[, 2])) * res^2
  Y <- matrix(rpois(length(lam), lam), n * n, J, dimnames = list(NULL, paste0("g", seq_len(J))))
  gx <- rep(seq_len(n) - 1, times = n); gy <- rep(seq_len(n) - 1, each = n)
  b <- structure(list(counts = Matrix::Matrix(Y, sparse = TRUE),
                      coords = data.frame(x = (gx + 0.5) * res, y = (gy + 0.5) * res, in_tissue = TRUE),
                      grid = list(nx = n, ny = n, bin_size = res, xmin = 0, ymin = 0), genes = colnames(Y)),
                 class = "binned_transcripts")
  off <- rff_offset(b, data.frame(domain = factor(as.vector(dom))))
  f <- fit_spatial_rff(b, n_factors = 2, lengthscales = c(8, 15), offset = off, basis = "grid", ard = 1, l1 = 0.5, max_iter = 30)
  expect_identical(f$basis, "grid"); expect_null(f$Omega)
  expect_equal(f$settings$l1, 0.5)
  expect_equal(nrow(f$field_grid), n * n)
  # rff_fields() interpolates the grid fields: exact at bin centres
  fx <- rff_fields(f, b$coords[c(1, 50, n * n), ])
  expect_equal(unname(as.matrix(fx)), unname(f$field_grid[c(1, 50, n * n), ]), tolerance = 1e-10)
  # fixed length scales stay at their initial values
  f2 <- fit_spatial_rff(b, n_factors = 2, lengthscales = c(8, 15), offset = off, basis = "grid", learn_lengthscales = FALSE, max_iter = 20)
  expect_equal(unname(f2$lengthscales), c(8, 15))
  # gene-specific L1 weights are accepted; wrong length is an error
  expect_s3_class(fit_spatial_rff(b, n_factors = 1, offset = off, basis = "grid", l1 = seq_len(J) / 10, max_iter = 5), "spatial_rff_fit")
  expect_error(fit_spatial_rff(b, n_factors = 1, offset = off, l1 = c(1, 2)), "l1")
  # analytic gradient of the grid-basis objective equals finite differences
  old <- options(cohalu.rff_debug = TRUE); on.exit(options(old))
  invisible(fit_spatial_rff(b, n_factors = 2, lengthscales = c(8, 15), offset = off, basis = "grid", ard = 1, l1 = 0.3, max_iter = 1))
  d <- get(".rff_dbg", envir = globalenv()); rm(".rff_dbg", envir = globalenv())
  set.seed(1); th <- d$theta0 + rnorm(length(d$theta0), sd = 0.05)
  ch <- c(d$idx$gamma[c(3, 400, 900)], d$idx$log_ell, d$idx$L[1:3], d$idx$alpha[1], d$idx$log_r[1])
  num <- vapply(ch, function(i) { a <- th; a[i] <- a[i] + 1e-5; z <- th; z[i] <- z[i] - 1e-5
    (d$objgrad(a)$value - d$objgrad(z)$value) / 2e-5 }, 0)
  expect_equal(d$objgrad(th)$grad[ch], num, tolerance = 1e-4)
})

test_that("sequential rff_factor_test() and rff_program_test()", {
  tx <- simulate_transcripts(size = 70, rate = 0.02, n_genes_per_set = 2, seed = 5)
  b <- bin_transcripts(tx, bin_size = 8)
  f <- fit_spatial_rff(b, n_factors = 2, n_features = 16, max_iter = 15, ard = 2)
  t1 <- rff_factor_test(f, b, n_boot = 3)
  ts <- rff_factor_test(f, b, n_boot = 3, sequential = TRUE, alpha = 1)
  expect_equal(attr(t1, "null_max"), attr(ts, "null_max"))    # first stage unchanged
  expect_true("p_sequential" %in% names(ts))
  ps <- ts$p_sequential[order(-ts$program_strength)]
  expect_true(all(diff(ps[!is.na(ps)]) >= 0))
  pt <- rff_program_test(f, b, n_boot = 5, n_components = 3)
  expect_equal(nrow(pt), 3)
  expect_true(all(pt$share > 0 & pt$share < 1))
  expect_true(pt$p[1] > 0 && pt$p[1] <= 1)
  expect_equal(dim(attr(pt, "scores")), c(sum(b$coords$in_tissue), 3))
  expect_equal(rownames(attr(pt, "loadings")), f$genes)
  pn <- rff_program_test(f, b, n_boot = 5, n_components = 2, sequential = FALSE)
  expect_false(anyNA(pn$p))
})

test_that("rff_program_test(null = \"shift\"): gene-shift surrogates", {
  set.seed(8); n <- 20; res <- 5; J <- 8
  dom <- cohalu:::.grf_fft(n, n, res, 20) > 0
  prof <- matrix(rnorm(J * 2), J, 2)
  own <- vapply(seq_len(J), function(j) as.vector(cohalu:::.grf_fft(n, n, res, 10)), numeric(n * n))  # gene-specific fields
  lam <- exp(log(0.3) + outer(as.vector(dom), prof[, 1]) + outer(!as.vector(dom), prof[, 2]) + 0.5 * own) * res^2
  Y <- matrix(rpois(length(lam), lam), n * n, J, dimnames = list(NULL, paste0("g", seq_len(J))))
  gx <- rep(seq_len(n) - 1, times = n); gy <- rep(seq_len(n) - 1, each = n)
  b <- structure(list(counts = Matrix::Matrix(Y, sparse = TRUE),
                      coords = data.frame(x = (gx + 0.5) * res, y = (gy + 0.5) * res, in_tissue = TRUE),
                      grid = list(nx = n, ny = n, bin_size = res, xmin = 0, ymin = 0), genes = colnames(Y)),
                 class = "binned_transcripts")
  off <- rff_offset(b, data.frame(domain = factor(as.vector(dom))))
  f <- fit_spatial_rff(b, n_factors = 2, lengthscales = 8, offset = off, basis = "grid", learn_lengthscales = FALSE, ard = 10, max_iter = 20)
  pt <- rff_program_test(f, b, n_boot = 19, n_components = 3, null = "shift")
  expect_identical(attr(pt, "null"), "shift")
  expect_equal(nrow(pt), 3)
  expect_true(all(pt$share[!is.na(pt$share)] > 0))
  ps <- pt$p[!is.na(pt$p)]
  expect_true(all(diff(ps) >= 0))                 # sequential p-values are cumulative maxima
  expect_length(attr(pt, "null_shares"), 3)
  pn <- rff_program_test(f, b, n_boot = 9, n_components = 2, null = "shift", sequential = FALSE)
  expect_false(anyNA(pn$p))
  expect_error(rff_program_test(f, b, null = "nope"))
  ph <- rff_program_test(f, b, n_boot = 9, n_components = 2, null = "shift", highpass = 20)
  expect_equal(nrow(ph), 2)
  expect_error(rff_program_test(f, b, n_boot = 3, highpass = 20), "highpass")
})

test_that("rff_offset(): data-driven covariates (kmeans, pca)", {
  tx <- simulate_transcripts(size = 90, rate = 0.02, n_genes_per_set = 3, seed = 9)
  b <- bin_transcripts(tx, bin_size = 6)
  nk <- sum(b$coords$in_tissue)
  o1 <- rff_offset(b, method = "kmeans", k = 4)
  expect_equal(dim(o1), c(nk, length(b$genes)))
  expect_true(all(is.finite(o1)))
  expect_equal(o1, rff_offset(b, method = "kmeans", k = 4))            # reproducible (seed)
  o2 <- rff_offset(b, method = "pca", k = 3)
  expect_equal(dim(o2), c(nk, length(b$genes)))
  # given covariates are added to the data-driven ones
  cv <- b$coords[b$coords$in_tissue, "x", drop = FALSE]
  expect_equal(dim(rff_offset(b, cv, method = "pca", k = 2)), c(nk, length(b$genes)))
  expect_error(rff_offset(b), "covariates")
  expect_error(rff_offset(b, method = "kmeans", k = 1), "k >= 2")
  # labels given as covariates still work positionally (backward compatible)
  expect_equal(dim(rff_offset(b, cv)), c(nk, length(b$genes)))
})

# small multi-window helper: 3 tissues sharing a 3-gene program
.joint_windows <- function(W = 3, n = 18, res = 5, J = 8, amp = 1.2) {
  lapply(seq_len(W), function(w) {
    set.seed(100 + w)
    dom <- cohalu:::.grf_fft(n, n, res, 20) > 0
    prog <- pmax(cohalu:::.grf_fft(n, n, res, 8) - 0.5, 0)
    prof <- matrix(rnorm(J * 2), J, 2)
    lam <- exp(log(0.3) + outer(as.vector(dom), prof[, 1]) + outer(!as.vector(dom), prof[, 2]) +
               amp * outer(as.vector(prog), as.numeric(seq_len(J) <= 3))) * res^2
    Y <- matrix(rpois(length(lam), lam), n * n, J, dimnames = list(NULL, paste0("g", seq_len(J))))
    gx <- rep(seq_len(n) - 1, times = n); gy <- rep(seq_len(n) - 1, each = n)
    b <- structure(list(counts = Matrix::Matrix(Y, sparse = TRUE),
                        coords = data.frame(x = (gx + 0.5) * res, y = (gy + 0.5) * res, in_tissue = TRUE),
                        grid = list(nx = n, ny = n, bin_size = res, xmin = 0, ymin = 0), genes = colnames(Y)),
                   class = "binned_transcripts")
    list(b = b, off = rff_offset(b, data.frame(domain = factor(as.vector(dom)))), prog = as.vector(prog))
  })
}

test_that("rff_offset(base =) adds known expectations and handles empty genes", {
  tx <- simulate_transcripts(size = 80, rate = 0.02, n_genes_per_set = 2, seed = 11)
  b <- bin_transcripts(tx, bin_size = 8)
  nk <- sum(b$coords$in_tissue); J <- length(b$genes)
  base <- matrix(rnorm(nk * J, sd = 0.2), nk, J, dimnames = list(NULL, b$genes))
  o0 <- rff_offset(b, base = base)                               # intercepts only
  expect_equal(dim(o0), c(nk, J))
  # intercept-only: the offset is base plus a constant per gene
  d <- o0 - base
  expect_lt(max(apply(d, 2, stats::sd)), 1e-8)
  cv <- b$coords[b$coords$in_tissue, "x", drop = FALSE]
  expect_equal(dim(rff_offset(b, cv, base = base)), c(nk, J))
  expect_equal(rff_offset(b, cv), rff_offset(b, cv, base = base * 0))
  expect_error(rff_offset(b, base = base[-1, ]), "one row")
  b0 <- b; b0$counts[, 1] <- 0
  expect_true(all(is.finite(rff_offset(b0, cv))))                # a gene without counts
})

test_that("rff_expected_offset() follows cell-type ownership", {
  tx <- simulate_transcripts(size = 80, rate = 0.03, n_genes_per_set = 3, seed = 12)
  tx$type <- ifelse(tx$x < 40, "left", "right")
  b <- bin_transcripts(tx, bin_size = 8)
  b2 <- b; b2$genes <- b$genes[1:3]; b2$counts <- b$counts[, 1:3]
  off <- rff_expected_offset(b2, tx)
  expect_equal(dim(off), c(nrow(b2$coords), 3))
  expect_true(all(is.finite(off)))
  expect_equal(colnames(off), b2$genes)
  expect_equal(ncol(attr(off, "type_counts")), 2)
  # usable directly as an offset and as a base
  expect_s3_class(fit_spatial_rff(b2, n_factors = 1, offset = off, n_features = 8, max_iter = 5), "spatial_rff_fit")
  expect_equal(dim(rff_offset(b2, base = off)), c(sum(b2$coords$in_tissue), 3))
  expect_error(rff_expected_offset(b2, tx, type_col = "nope"), "not found")
})

test_that("fit_spatial_rff(loadings =) keeps loadings fixed and reports amplitudes", {
  ws <- .joint_windows(W = 1)
  L0 <- matrix(c(1, 1, 1, rep(0, 5)) / sqrt(3), ncol = 1, dimnames = list(paste0("g", 1:8), "p1"))
  f <- fit_spatial_rff(ws[[1]]$b, offset = ws[[1]]$off, loadings = L0, basis = "grid", lengthscales = 8,
                       learn_lengthscales = FALSE, max_iter = 40, factor_init = "residual_pca")
  expect_equal(ncol(f$L), 1)
  expect_equal(unname(f$L[, 1] / f$amplitude), unname(L0[, 1]), tolerance = 1e-8)
  expect_gt(f$amplitude, 0)
  expect_gt(abs(cor(f$field_grid[, "factor1"], ws[[1]]$prog)), 0.3)
  expect_null(fit_spatial_rff(ws[[1]]$b, n_factors = 1, n_features = 8, max_iter = 3)$amplitude)
  # partial gene names are filled with 0
  f2 <- fit_spatial_rff(ws[[1]]$b, loadings = L0[1:3, , drop = FALSE], n_features = 8, max_iter = 3)
  expect_equal(unname(f2$L[4:8, 1]), rep(0, 5))
  expect_error(fit_spatial_rff(ws[[1]]$b, loadings = matrix(1, 3, 1)), "one row per gene")
})

test_that("rff_program_test_joint(): shared loadings, confirmatory mode and joint fit", {
  ws <- .joint_windows(W = 3)
  bs <- lapply(ws, `[[`, "b")
  fits <- lapply(ws, function(w) fit_spatial_rff(w$b, n_factors = 2, offset = w$off, basis = "grid", lengthscales = 8,
                                                 learn_lengthscales = FALSE, ard = 10, max_iter = 30))
  jt <- rff_program_test_joint(fits, bs, n_boot = 9, n_components = 3)
  expect_identical(attr(jt, "null"), "shift_joint")
  expect_equal(nrow(jt), 3)
  V <- attr(jt, "loadings")
  expect_equal(unname(colSums(V^2)), rep(1, 3), tolerance = 1e-8)
  # PC1 contrasts the program genes with the others (the shared direction is removed)
  expect_equal(length(unique(sign(V[1:3, 1]))), 1)
  expect_gt(-sign(V[1, 1]) * mean(V[4:8, 1]), 0)
  expect_length(attr(jt, "scores"), 3)
  expect_equal(dim(attr(jt, "window_share")), c(3, 3))
  ps <- jt$p[!is.na(jt$p)]; expect_true(all(diff(ps) >= 0))
  # results do not depend on n_cores; numeric weights accepted
  skip_on_os("windows")
  expect_equal(rff_program_test_joint(fits, bs, n_boot = 5, n_components = 2, n_cores = 2)$share,
               rff_program_test_joint(fits, bs, n_boot = 5, n_components = 2)$share)
  expect_equal(nrow(rff_program_test_joint(fits, bs, n_boot = 3, n_components = 2, weights = c(1, 1, 0.5))), 2)
  expect_error(rff_program_test_joint(fits, bs, weights = c(1, -1, 1)), "non-negative")
  cf <- rff_program_test_joint(fits, bs, loadings = V[, 1, drop = FALSE], n_boot = 9)
  expect_equal(nrow(cf$windows), 3)
  expect_equal(nrow(cf$combined), 1)
  expect_true(all(cf$windows$p > 0 & cf$windows$p <= 1))
  expect_gt(cf$combined$z, 0)
  expect_error(rff_program_test_joint(fits, bs[1:2]), "same length")
  jf <- fit_spatial_rff_joint(bs, lapply(ws, `[[`, "off"), loadings = V[, 1, drop = FALSE], n_iter = 1,
                              lengthscales = 8, max_iter = 20)
  expect_s3_class(jf, "spatial_rff_joint")
  expect_equal(dim(jf$amplitude), c(3, 1))
  expect_equal(unname(sqrt(colSums(jf$loadings^2))), 1, tolerance = 1e-8)
  expect_length(jf$fits, 3)
  expect_output(print(jf), "shared programs")
})

test_that("rff_report() shows one window of a joint fit with the joint test", {
  skip_on_cran()
  ws <- .joint_windows(W = 2)
  bs <- lapply(ws, `[[`, "b"); names(bs) <- c("w1", "w2")
  fits <- lapply(ws, function(w) fit_spatial_rff(w$b, n_factors = 2, offset = w$off, basis = "grid", lengthscales = 8,
                                                 learn_lengthscales = FALSE, ard = 10, max_iter = 20))
  names(fits) <- names(bs)
  jt <- rff_program_test_joint(fits, bs, n_boot = 5, n_components = 2)
  jf <- fit_spatial_rff_joint(bs, lapply(ws, `[[`, "off"), loadings = attr(jt, "loadings")[, 1, drop = FALSE],
                              n_iter = 1, lengthscales = 8, max_iter = 15)
  f <- file.path(tempdir(), "joint_report.html")
  out <- rff_report(jf, bs, test = jt, file = f, window = "w2", export = FALSE, open = FALSE)
  expect_true(file.exists(f))
  h <- paste(readLines(f, warn = FALSE), collapse = "\n")
  expect_true(grepl("shared loadings over 2 windows", h, fixed = TRUE))
  expect_error(rff_report(jf, bs, window = "nope", file = f, export = FALSE, open = FALSE), "window")
})
