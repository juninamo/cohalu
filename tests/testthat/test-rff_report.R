rr_setup <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      tx <- simulate_transcripts(size = 100, rate = 0.02, n_genes_per_set = 3, seed = 2)
      b <- bin_transcripts(tx, bin_size = 6)
      fit <- fit_spatial_rff(b, n_factors = 2, ard = 2, n_features = 24, max_iter = 40)
      cache <<- list(tx = tx, b = b, fit = fit)
    }
    cache
  }
})

rr_read <- function(f) paste(readLines(f, warn = FALSE), collapse = "\n")

# base64 decoder (test helper)
jsonlite_free_b64dec <- function(s) {
  tbl <- c(LETTERS, letters, 0:9, "+", "/")
  ch <- strsplit(s, "")[[1]]; npad <- sum(ch == "="); ch[ch == "="] <- "A"
  v <- match(ch, tbl) - 1L
  m <- matrix(v, 4)
  w <- m[1, ] * 262144 + m[2, ] * 4096 + m[3, ] * 64 + m[4, ]
  r <- as.vector(rbind(w %/% 65536, (w %/% 256) %% 256, w %% 256))
  as.raw(r[seq_len(length(r) - npad)])
}

test_that("internal base64 and JSON helpers are correct", {
  expect_equal(.rr_b64(charToRaw("Man")), "TWFu")
  expect_equal(.rr_b64(charToRaw("Ma")), "TWE=")
  expect_equal(.rr_b64(charToRaw("M")), "TQ==")
  expect_equal(.rr_b64(as.raw(0:255)), "AAECAwQFBgcICQoLDA0ODxAREhMUFRYXGBkaGxwdHh8gISIjJCUmJygpKissLS4vMDEyMzQ1Njc4OTo7PD0+P0BBQkNERUZHSElKS0xNTk9QUVJTVFVWV1hZWltcXV5fYGFiY2RlZmdoaWprbG1ub3BxcnN0dXZ3eHl6e3x9fn+AgYKDhIWGh4iJiouMjY6PkJGSk5SVlpeYmZqbnJ2en6ChoqOkpaanqKmqq6ytrq+wsbKztLW2t7i5uru8vb6/wMHCw8TFxsfIycrLzM3Oz9DR0tPU1dbX2Nna29zd3t/g4eLj5OXm5+jp6uvs7e7v8PHy8/T19vf4+fr7/P3+/w==")
  expect_equal(.rr_json(list(a = 1, b = I(c(1, NA)), c = "x\"</", d = TRUE, e = NULL)),
               "{\"a\":1,\"b\":[1,null],\"c\":\"x\\\"<\\/\",\"d\":true,\"e\":null}")
  q <- .rr_quant(c(0, 1, NA), 16L)
  r <- readBin(jsonlite_free_b64dec(q$b64), "integer", n = 3, size = 2, signed = FALSE, endian = "little")
  expect_equal(r, c(1L, 65535L, 0L))
})

test_that("rff_report writes a complete report with tests, cells and known programs", {
  skip_on_cran()
  s <- rr_setup()
  pt <- rff_program_test(s$fit, s$b, n_boot = 19, null = "shift")
  ft <- rff_factor_test(s$fit, s$b, n_boot = 4)
  set.seed(1)
  cells <- data.frame(x = runif(200, 0, 100), y = runif(200, 0, 100),
                      type = sample(c("T", "B", "Mac"), 200, TRUE))
  f <- file.path(tempfile("rr"), "report.html")
  out <- suppressWarnings(rff_report(s$fit, s$b, test = list(shift = pt, factor = ft), file = f, cells = cells,
                                     cell_label_col = "type", known_programs = list(A = c("A_1", "A_2", "A_3")),
                                     meta = list(sample = "sim"), offset_info = list(method = "none (test)"), open = FALSE))
  expect_equal(out, normalizePath(f, mustWork = FALSE))
  expect_true(file.exists(f))
  h <- rr_read(f)
  for (id in c("id=\"overview\"", "id=\"input\"", "id=\"findings\"", "id=\"outputs\"", "id=\"programs\"", "id=\"details\"",
               "id=\"maps\"", "id=\"loadings\"", "id=\"cells\"", "id=\"known\"", "id=\"testsec\"",
               "id=\"methods\"", "id=\"reproduce\"", "Key findings", "Input data", "Output files",
               "How to reproduce", "rff_program_test", "Gene-shift", "sim"))
    expect_true(grepl(id, h, fixed = TRUE), info = id)
  expect_false(grepl("{{", h, fixed = TRUE))
  for (id in c("Detection axes vs program maps", "class=\"help\"", "you normally don't need these",
               "(fit factor1)", "(program test PC1)", "program test"))
    expect_true(grepl(id, h, fixed = TRUE), info = id)
  expect_false(grepl("<script src=|<link [^>]*href=\"http", h))
  # key sections come before the detailed panels
  expect_lt(regexpr("id=\"findings\"", h), regexpr("id=\"maps\"", h))
  expect_lt(regexpr("id=\"outputs\"", h), regexpr("id=\"programs\"", h))
  expect_lt(regexpr("id=\"programs\"", h), regexpr("id=\"details\"", h))
  expect_lt(file.size(f), 3e6)
  # embedded JSON parses (checked by structure: balanced and closes the script)
  js <- sub(".*<script id=\"D\" type=\"application/json\">", "", h)
  js <- sub("</script>.*", "", js)
  expect_true(startsWith(js, "{") && endsWith(js, "}"))
  # exported files
  od <- file.path(dirname(f), "report_files")
  fl <- list.files(od)
  for (x in c("programs.csv", "display_names.csv", "factors.csv", "loadings.csv", "fields.csv", "cells.csv", "test1_program.csv", "test1_null.csv",
              "test2_factor.csv", "known_programs.csv", "fit.rds", "settings.json", "reproduce.R"))
    expect_true(x %in% fl, info = x)
  fac <- utils::read.csv(file.path(od, "factors.csv"))
  expect_equal(nrow(fac), 2)
  expect_equal(fac$display_name, c("M1", "M2"))
  expect_true("display_name" %in% names(utils::read.csv(file.path(od, "test1_program.csv"))))
  expect_true("D1" %in% utils::read.csv(file.path(od, "display_names.csv"))$display_name)
  ld <- utils::read.csv(file.path(od, "loadings.csv"))
  expect_equal(nrow(ld), length(s$fit$genes))
  expect_true(all(c("factor1", "factor2", "PC1") %in% names(ld)))
  fd <- utils::read.csv(file.path(od, "fields.csv"))
  expect_equal(nrow(fd), sum(s$b$coords$in_tissue))
  expect_equal(nrow(utils::read.csv(file.path(od, "cells.csv"))), 200)
  expect_s3_class(readRDS(file.path(od, "fit.rds")), "spatial_rff_fit")
  expect_true(grepl("report_files/factors.csv", h, fixed = TRUE))
  expect_true(any(grepl("fit_spatial_rff", readLines(file.path(od, "reproduce.R")))))
})

test_that("rff_report works without test, cells or known programs and without export", {
  skip_on_cran()
  s <- rr_setup()
  f <- tempfile(fileext = ".html")
  rff_report(s$fit, s$b, file = f, export = FALSE, open = FALSE)
  h <- rr_read(f)
  expect_true(grepl("No test", h, fixed = TRUE))
  expect_true(grepl("class=\"hidden\"><h2>Cells", h, fixed = TRUE))
  expect_false(dir.exists(sub("\\.html$", "_files", f)))
  expect_true(grepl("export = FALSE", h, fixed = TRUE))
})

test_that("rff_report handles a grid-basis residual fit, older fits and large grids", {
  skip_on_cran()
  s <- rr_setup()
  cv <- s$b$coords[s$b$coords$in_tissue, c("x", "y")]
  off <- rff_offset(s$b, cv)
  fit <- fit_spatial_rff(s$b, n_factors = 2, offset = off, basis = "grid", lengthscales = 12,
                         learn_lengthscales = FALSE, ard = 2, max_iter = 30)
  pt <- rff_program_test(fit, s$b, n_boot = 9, null = "parametric")
  f <- tempfile(fileext = ".html")
  expect_warning(rff_report(fit, s$b, test = pt, file = f, max_pixels = 50, open = FALSE,
                            offset_info = list(method = "covariates", covariates = cv)),
                 "block-averaged")
  h <- rr_read(f)
  expect_true(grepl("Deviance explained by the offset", h, fixed = TRUE))
  expect_true(grepl("Only rff_program_test() was run", h, fixed = TRUE))
  expect_true(grepl("Calibration", h, fixed = TRUE))
  old <- s$fit
  old$settings <- NULL; old$factor_strength <- NULL; old$program_strength <- NULL; old$center <- NULL
  old$basis <- NULL; old$offset <- NULL
  f2 <- tempfile(fileext = ".html")
  expect_silent(rff_report(old, s$b, file = f2, export = FALSE, open = FALSE))
  expect_true(file.exists(f2))
  expect_error(rff_report(list(), s$b, file = f2, open = FALSE), "fit_spatial_rff")
})

test_that("matching of detection axes and program maps gives Confirmed / Candidate / Not supported", {
  set.seed(3)
  n <- 20; genes <- paste0("g", 1:6)
  f1 <- rnorm(n * n); f2 <- rnorm(n * n); f3 <- rnorm(n * n)
  L <- cbind(factor1 = c(1, 1, -1, 0, 0, 0), factor2 = c(0, 0, 0, 1, -1, 0.5), factor3 = rep(0.5, 6))
  rownames(L) <- genes
  fit <- list(L = L, genes = genes, in_tissue = rep(TRUE, n * n), lengthscales = c(8, 8, 8),
              field_grid = cbind(factor1 = f1, factor2 = f2, factor3 = f3))
  pt <- data.frame(component = c("PC1", "PC2", "PC3"), share = c(0.3, 0.2, 0.1), p = c(0.01, 0.02, 0.5))
  attr(pt, "scores") <- cbind(PC1 = f1 + rnorm(n * n, sd = 0.3), PC2 = rnorm(n * n), PC3 = rnorm(n * n))
  attr(pt, "loadings") <- cbind(PC1 = c(1, 1, -1, 0, 0, 0), PC2 = c(0, 0, 0, 0, 1, 1), PC3 = 0) / 1.7
  dimnames(attr(pt, "loadings")) <- list(genes, c("PC1", "PC2", "PC3"))
  m <- .rff_match_programs(fit, pt, alpha = 0.05, min_r = 0.5)
  expect_equal(m$status[m$program %in% "P1"], "Confirmed")
  expect_equal(m$program_map[m$program %in% "P1"], "factor1")
  expect_equal(m$status[m$program %in% "P2"], "Candidate")
  expect_true("Not supported" %in% m$status)
  expect_equal(m$program_map[m$status == "Cellularity/technical"], "factor3")
  expect_equal(ncol(attr(m, "loadings")), nrow(m))
  expect_equal(nrow(attr(m, "fields")), n * n)
  m0 <- .rff_match_programs(fit)
  expect_true(all(m0$status %in% c("Exploratory", "Cellularity/technical")))
})
