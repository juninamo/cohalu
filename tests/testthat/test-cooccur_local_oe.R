test_that("closed-form expected pair counts match the permutation mean", {
  set.seed(1)
  n <- 400
  d <- data.frame(x = runif(n, 0, 200), y = runif(n, 0, 200),
                  cell_type = sample(c("A", "B", "C"), n, TRUE, prob = c(.2, .3, .5)))
  nb <- cohalu:::.radius_neighbours(as.matrix(d[, c("x", "y")]), 20, 100)
  ex <- cohalu:::.pairs_expected(d$cell_type, nb, "A", "B")
  # permutation mean of the total number of pairs, keeping each cell's own label
  perm_tot <- replicate(300, {
    lab <- d$cell_type; others <- sample(lab)            # approximate: full shuffle
    sum(cohalu:::.pairs_expected(others, nb, "A", "B")$pairs)
  })
  expect_equal(sum(ex$expected), mean(perm_tot), tolerance = 0.03)
})

test_that("cooccur_local_oe is ~0 under randomness and positive when planted", {
  set.seed(2)
  n <- 1200
  d <- data.frame(x = runif(n, 0, 400), y = runif(n, 0, 400),
                  cell_type = sample(c("A", "B", "other"), n, TRUE, prob = c(.1, .1, .8)))
  rownames(d) <- paste0("c", seq_len(n))
  r0 <- cooccur_local_oe(d, "A", "B", radius = 25)
  expect_lt(abs(attr(r0, "section_log2_oe")), 0.25)
  expect_true(all(c("n_x", "n_y", "pairs", "expected", "local_log2_oe") %in% colnames(r0)))
  sim <- generate_sim(close_ratio = 1, n_types = 6, n_cells = 900, max_loc = 400,
                      test_type = "circle", distance_param = 15, seed = 3)
  rownames(sim) <- paste0("c", seq_len(nrow(sim)))
  r1 <- cooccur_local_oe(sim, "cell_type_1", "cell_type_2", radius = 25, n_perms = 49)
  expect_gt(attr(r1, "section_log2_oe"), 0.5)
  expect_true(all(r1$p > 0 & r1$p <= 1))
})

test_that("cooccur_local_per_sample reports an abundance-adjusted log2_oe", {
  skip_on_cran()
  df <- generate_sim_groups(n_samples_per_group = 2,
                            group_close_ratio = list(case = 0.8, control = 0.1),
                            n_types = 5, n_cells = 400, max_loc = 300,
                            test_type = "distribute", distance_param = 8, seed = 4)
  res <- cooccur_local_per_sample(df, "sample_id", "group", "cell_type",
                                  "cell_type_1", "cell_type_2", patient_key = "patient",
                                  neighbors.k = 60, radius = 20, maxnsteps = 0)
  expect_true("log2_oe" %in% colnames(res))
  expect_gt(mean(res$log2_oe[res$group == "case"]), mean(res$log2_oe[res$group == "control"]))
})
