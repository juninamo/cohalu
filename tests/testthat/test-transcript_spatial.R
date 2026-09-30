test_that(".grf_fft gives unit-variance fields with an RBF covariance", {
  set.seed(1)
  f <- cohalu:::.grf_fft(200, 200, delta = 1, ell = 5)
  expect_equal(var(as.vector(f)), 1, tolerance = 0.15)
  lag_cor <- cor(as.vector(f[1:195, ]), as.vector(f[6:200, ]))
  expect_equal(lag_cor, exp(-0.5), tolerance = 0.08)
})

test_that("bin_transcripts keeps every transcript", {
  tx <- simulate_transcripts(size = 60, rate = 0.02, n_genes_per_set = 2, seed = 2)
  b <- bin_transcripts(tx, bin_size = 5)
  expect_s3_class(b, "binned_transcripts")
  expect_equal(sum(b$counts), nrow(tx))
  expect_equal(nrow(b$counts), b$grid$nx * b$grid$ny)
  expect_setequal(b$genes, unique(tx$gene))
})

test_that("pcf_cross is ~0 on the log scale under complete spatial randomness", {
  set.seed(3)
  n <- 6000
  tx <- data.frame(x = runif(n, 0, 300), y = runif(n, 0, 300),
                   gene = sample(c("a", "b"), n, replace = TRUE))
  b <- bin_transcripts(tx, bin_size = 5, tissue_radius = Inf)
  ab <- pcf_cross(b, "a", "b", r_max = 40)
  aa <- pcf_cross(b, "a", r_max = 40)          # self-pairs must be removed
  expect_lt(max(abs(ab$log_g)), 0.15)
  expect_lt(max(abs(aa$log_g)), 0.15)
})

test_that("pcf_cross and the model recover planted co-localization", {
  skip_on_cran()
  tx <- simulate_transcripts(size = 200, rate = 0.03, n_genes_per_set = 3,
                             coloc = c(A = 1.2, B = 1.2), set_sd = 0.5, seed = 4)
  sets <- with(unique(tx[, c("gene", "gene_set")]), split(gene, gene_set))
  b <- bin_transcripts(tx, bin_size = 5, tissue_radius = Inf)
  ab <- pcf_cross(b, sets$A, sets$B, r_max = 20)
  ac <- pcf_cross(b, sets$A, sets$C, r_max = 20)
  expect_gt(ab$log_g[2], ac$log_g[2] + 0.3)
  fit <- fit_spatial_rff(b, n_factors = 3, n_features = 32, max_iter = 150)
  expect_s3_class(fit, "spatial_rff_fit")
  expect_equal(dim(fit$L), c(9L, 3L))
  m_ab <- rff_pair_correlation(fit, sets$A, sets$B, r = 10, type = "composition")$log_g
  m_ac <- rff_pair_correlation(fit, sets$A, sets$C, r = 10, type = "composition")$log_g
  expect_gt(m_ab, m_ac)
  expect_gt(m_ab, 0.15)
  g_ab <- rff_pair_correlation(fit, sets$A, sets$B, r = 10, type = "composition",
                               method = "gaussian")$log_g
  expect_gt(g_ab, 0.3)
})

test_that("colocalization_per_sample feeds compare_groups", {
  skip_on_cran()
  tx <- simulate_transcripts_groups(n_samples_per_group = 3, size = 100, rate = 0.03,
                                    n_genes_per_set = 2, seed = 5)
  sets <- with(unique(tx[, c("gene", "gene_set")]), split(gene, gene_set))
  res <- colocalization_per_sample(tx, "sample_id", "group", "patient", sets,
                                   pairs = list(c("A", "B")), r = c(8, 16))
  expect_equal(nrow(res), 6L * 2L)
  expect_true(all(c("cluster_i", "cluster_j", "r", "log_g", "group") %in% colnames(res)))
  cmp <- compare_groups(res, value = "log_g", method = "wilcox", ref_group = "control",
                        pair_keys = c("cluster_i", "cluster_j", "r"))
  expect_equal(nrow(cmp), 2L)
})

test_that("fit_spatial_rff supports a smoothed-total offset", {
  tx <- simulate_transcripts(size = 80, rate = 0.03, n_genes_per_set = 2, seed = 6)
  b <- bin_transcripts(tx, bin_size = 5, tissue_radius = Inf)
  fit <- fit_spatial_rff(b, n_factors = 2, n_features = 16, max_iter = 30,
                         offset = "smoothed_total", offset_genes = c("C_1", "C_2"))
  expect_false(fit$has_density)
  expect_identical(fit$offset, "smoothed_total")
  expect_error(fit_spatial_rff(b, offset = "smoothed_total", offset_genes = "nope"),
               "offset_genes")
  lg <- rff_pair_correlation(fit, c("A_1", "A_2"), c("B_1", "B_2"), r = 10, type = "composition")
  expect_true(is.finite(lg$log_g))
})

test_that("pcf_matrix agrees with pcf_cross", {
  tx <- simulate_transcripts(size = 100, rate = 0.02, n_genes_per_set = 2, seed = 12)
  sets <- split(attr(tx, "truth")$genes, attr(tx, "truth")$set_of)
  b <- bin_transcripts(tx, bin_size = 5)
  m <- pcf_matrix(b, sets, r_max = 30)
  expect_equal(nrow(m), 6L * length(unique(m$r)))
  ab <- pcf_cross(b, sets$A, sets$B, r_max = 30)
  ab_rel <- pcf_cross(b, sets$A, sets$B, r_max = 30, relative = TRUE)
  aa <- pcf_cross(b, sets$A, r_max = 30)
  mm <- m[m$cluster_i == "A" & m$cluster_j == "B", ]
  expect_equal(mm$log_g, ab$log_g, tolerance = 1e-8)
  expect_equal(mm$log_g_rel, ab_rel$log_g, tolerance = 1e-8)
  expect_equal(m$log_g[m$cluster_i == "A" & m$cluster_j == "A"], aa$log_g, tolerance = 1e-8)
})

test_that("pcf_matrix accepts any labelled points, e.g. pixel-level factors", {
  set.seed(1); n <- 3000
  cx <- runif(15, 0, 300); cy <- runif(15, 0, 300); k <- sample(15, n, TRUE)
  px <- data.frame(X = c(cx[k] + rnorm(n, 0, 8), runif(n, 0, 300)),
                   Y = c(cy[k] + rnorm(n, 0, 8), runif(n, 0, 300)),
                   K1 = c(sample(c("1", "2"), n, TRUE), rep("3", n)))
  b <- bin_transcripts(px, bin_size = 4, x_col = "X", y_col = "Y", gene_col = "K1")
  pm <- pcf_matrix(b, list(F1 = "1", F2 = "2", F3 = "3"), r_max = 30)
  short <- pm[pm$r > 0 & pm$r <= 20 & is.finite(pm$log_g_rel), ]
  expect_gt(mean(short$log_g_rel[short$cluster_i == "F1" & short$cluster_j == "F2"]), 0.5)
  expect_lt(mean(short$log_g_rel[short$cluster_i == "F1" & short$cluster_j == "F3"]), -0.3)
})

test_that("gene-level co-localization modules recover planted gene sets", {
  tx <- simulate_transcripts(size = 300, rate = 0.02, coloc = c(A = 0.8, B = 0.8), seed = 11)
  truth <- attr(tx, "truth")
  b <- bin_transcripts(tx, bin_size = 4, tissue_radius = Inf)
  M <- colocalization_gene_matrix(b, radius = 12)
  expect_true(isSymmetric(unname(M)))
  sets <- split(truth$genes, truth$set_of)
  ab <- mean(M[sets$A, sets$B]); ac <- mean(M[sets$A, sets$C])
  expect_gt(ab, ac)                                     # A and B share a niche
  mods <- colocalization_modules(M, n_modules = 3)
  tab <- table(truth$set_of[match(mods$modules$gene, truth$genes)], mods$modules$module)
  expect_equal(sum(apply(tab, 1, max)), nrow(mods$modules))   # each set in one module
  en <- module_enrichment(mods, sets)
  expect_equal(sum(en$padj < 0.05), 3)
})

test_that("profile clustering splits a co-localizing module into its sub-populations", {
  # A and B share a niche (one co-localizing group), each with its own field
  tx <- simulate_transcripts(gene_sets = c("A", "B", "C", "D"), n_genes_per_set = 6, size = 240,
                             rate = 0.03, coloc = c(A = 1.5, B = 1.5), set_sd = 0.6, seed = 1)
  truth <- attr(tx, "truth"); sets <- split(truth$genes, truth$set_of)
  b <- bin_transcripts(tx, bin_size = 4, tissue_radius = Inf)
  M <- colocalization_gene_matrix(b, radius = 12)
  set_of <- function(g) truth$set_of[match(g, truth$genes)]
  # average linkage (default) merges A and B into one module
  mods <- colocalization_modules(M, n_modules = 3)
  m1 <- mods$modules$gene[mods$modules$module %in% "M1"]
  expect_setequal(m1, c(sets$A, sets$B))
  expect_setequal(mods$modules$gene[mods$modules$module %in% "M1"],
                  colocalization_modules(M)$modules$gene[colocalization_modules(M)$modules$module %in% "M1"])
  # the profile split of M1 separates them, k chosen by silhouette
  sub <- colocalization_submodules(M, mods, module = "M1")
  expect_equal(sub$k, 2)
  expect_equal(nrow(sub$silhouette), 7)
  expect_gt(sub$silhouette$mean_silhouette[sub$silhouette$k == 2], 0.5)
  tab <- table(set_of(sub$modules$gene), sub$modules$module)
  expect_equal(sum(apply(tab, 1, max)), 12)              # each set in one sub-module
  expect_true(all(apply(tab, 2, max) == 6))
  expect_setequal(sub$summary$module, c("M1.1", "M1.2"))
  expect_true(all(sub$modules$parent == "M1"))
  expect_true(all(sub$summary$mean_oe > sub$silhouette$between_oe[sub$silhouette$k == 2]))
  expect_equal(sub$universe, rownames(M))
  # module_enrichment() works on the sub-modules (background: all genes of M)
  en <- module_enrichment(sub, sets)
  sig <- en[en$padj < 0.05, ]
  expect_equal(nrow(sig), 2)
  expect_setequal(sig$gene_set, c("A", "B"))
  expect_true(all(sig$set_size == 6))
  # fixed k and a gene vector also work
  sub2 <- colocalization_submodules(M, m1, k = 2)
  expect_equal(sort(sub2$summary$size), c(6, 6))
  expect_error(colocalization_submodules(M, mods, module = "M9"), "not found")
  # method = "profile" for the whole matrix recovers all four sets
  p4 <- colocalization_modules(M, method = "profile", n_modules = 4)
  tab4 <- table(set_of(p4$modules$gene), p4$modules$module)
  expect_true(all(apply(tab4, 1, max) == 6) && all(apply(tab4, 2, max) == 6))
  pa <- colocalization_modules(M, method = "profile")
  expect_equal(pa$method, "profile")
  expect_true(all(c("k", "mean_silhouette", "within_oe", "between_oe") %in% names(pa$silhouette)))
})

test_that("gene-level O/E is ~0 when gene labels are random", {
  set.seed(4); n <- 20000
  tx <- data.frame(x = runif(n, 0, 300), y = runif(n, 0, 300))
  tx$x <- tx$x + ifelse(tx$x > 150, 0, 0)               # uneven density: denser left half
  keep <- tx$x > 150 | runif(n) < 0.4; tx <- tx[keep, ]
  tx$gene <- sample(paste0("g", 1:6), nrow(tx), TRUE)
  b <- bin_transcripts(tx, bin_size = 4, tissue_radius = 3)
  M <- colocalization_gene_matrix(b, radius = 15, min_count = 10)
  expect_lt(max(abs(M)), 0.1)
})
