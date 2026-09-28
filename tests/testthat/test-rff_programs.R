# Small tissue: two known domains (given as covariates) plus, optionally, a
# 4-gene minor program (genes g1-g4) in ~10 um spots.
rp_sim <- function(seed, amp = 1.5, n = 40, res = 5, J = 20) {
  set.seed(seed)
  dom <- .grf_fft(n, n, res, 40) > 0
  minor <- pmax(.grf_fft(n, n, res, 10) - 1, 0)
  prof <- matrix(stats::rnorm(J * 2), J, 2)
  eta <- log(0.1) + outer(as.vector(dom), prof[, 1]) + outer(!as.vector(dom), prof[, 2]) +
    amp * outer(as.vector(minor), as.numeric(seq_len(J) <= 4))
  Y <- matrix(stats::rpois(length(eta), exp(eta) * res^2), n * n, J)
  genes <- paste0("g", seq_len(J))
  gx <- rep(seq_len(n) - 1, times = n); gy <- rep(seq_len(n) - 1, each = n)
  b <- structure(list(counts = Matrix::Matrix(Y, sparse = TRUE, dimnames = list(NULL, genes)),
                      coords = data.frame(x = (gx + 0.5) * res, y = (gy + 0.5) * res, in_tissue = TRUE),
                      grid = list(nx = n, ny = n, bin_size = res, xmin = 0, ymin = 0), genes = genes),
                 class = "binned_transcripts")
  list(b = b, cv = data.frame(domain = factor(as.vector(dom))))
}

test_that("rff_programs confirms a planted program and nothing in null tissue", {
  skip_on_cran()
  s <- rp_sim(1)
  r <- rff_programs(s$b, offset = s$cv, n_factors = 4, n_boot = 19, fit_args = list(max_iter = 100))
  expect_s3_class(r, "rff_programs")
  p <- r$programs
  expect_true(any(p$status == "Confirmed"))
  k <- p$key[p$status == "Confirmed"][1]
  ld <- attr(p, "loadings")[, k]
  expect_gte(sum(names(sort(abs(ld), decreasing = TRUE))[1:4] %in% paste0("g", 1:4)), 3)
  expect_true(all(c("program", "status", "p", "detection_axis", "program_map", "r", "lengthscale",
                    "uniform_share", "top_up", "top_down", "n_genes") %in% names(p)))
  expect_equal(r$offset_info$method, "covariates")
  expect_output(print(r), "confirmed")
  expect_output(print(summary(r)), "Programs by status")
  expect_false(is.null(as.data.frame(r)$status))
  expect_null(attr(as.data.frame(r), "loadings"))
  for (seed in 1:2) {
    r0 <- rff_programs(rp_sim(seed, amp = 0)$b, offset = rp_sim(seed, amp = 0)$cv, n_factors = 4, n_boot = 19,
                       fit_args = list(max_iter = 100))
    expect_false(any(r0$programs$status == "Confirmed"), info = seed)
  }
  # reuse of a fit gives the same programs
  r2 <- rff_programs(s$b, fit = r$fit, n_boot = 19)
  expect_equal(r2$programs$status, p$status)
  expect_identical(r2$fit, r$fit)
  # report from the object
  f <- tempfile(fileext = ".html")
  rff_report(r, s$b, file = f, open = FALSE)
  h <- paste(readLines(f), collapse = "\n")
  expect_true(grepl("id=\"programs\"", h, fixed = TRUE))
  expect_true(grepl("res &lt;- rff_programs(", h, fixed = TRUE))
  pc <- utils::read.csv(file.path(sub("\\.html$", "_files", f), "programs.csv"))
  expect_true("Confirmed" %in% pc$status)
})

test_that("rff_programs offset shortcuts, factor test only, and report argument", {
  skip_on_cran()
  s <- rp_sim(2)
  fa <- list(max_iter = 30)
  for (o in list("area", "kmeans", "pca", "smoothed_total")) {
    r <- rff_programs(s$b, offset = o, n_factors = 2, null = "none", k = 3, fit_args = fa)
    expect_equal(r$offset_info$method, if (o == "smoothed_total") "covariates" else o, info = o)
    expect_true(all(r$programs$status %in% c("Exploratory", "Cellularity/technical")), info = o)
  }
  off <- rff_offset(s$b, s$cv)
  r <- rff_programs(s$b, offset = off, n_factors = 2, null = "none", fit_args = fa)
  expect_true(is.matrix(r$fit$offset_matrix))
  expect_error(rff_programs(s$b, offset = "nonsense"), "Unknown")
  expect_warning(rff_programs(s$b, fit = fit_spatial_rff(s$b, n_factors = 2, n_features = 16, max_iter = 20),
                              null = "none"), "recommended")
  f <- tempfile(fileext = ".html")
  r <- rff_programs(s$b, offset = s$cv, n_factors = 2, null = "none", factor_test = TRUE, factor_n_boot = 4,
                    fit_args = fa, report = f, open = FALSE)
  expect_true(file.exists(f))
  expect_false(is.null(r$factor_test))
  expect_true(all(r$programs$detected_by[r$programs$status == "Confirmed"] == "factor test"))
})
