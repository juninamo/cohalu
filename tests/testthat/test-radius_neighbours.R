# Radius-based statistics must use every cell within the radius by default.

brute_pairs <- function(xy, radius) {
  d <- as.matrix(dist(xy))
  diag(d) <- Inf
  which(d <= radius, arr.ind = TRUE)
}

brute_counts <- function(xy, cl, lv, radius) {
  w <- brute_pairs(xy, radius)
  m <- matrix(0, length(lv), length(lv),
              dimnames = list(paste0("Cluster", lv), paste0("Cluster", lv)))
  for (r in seq_len(nrow(w))) {
    a <- match(cl[w[r, 1]], lv)
    b <- match(cl[w[r, 2]], lv)
    m[a, b] <- m[a, b] + 1
  }
  m
}

dense_cells <- function(n = 600, side = 120, seed = 1) {
  set.seed(seed)
  list(xy = cbind(x = runif(n, 0, side), y = runif(n, 0, side)),
       cl = sample(c("A", "B", "C"), n, TRUE))
}

test_that("radius counts match brute force at any neighbourhood size", {
  s <- dense_cells()
  lv <- c("A", "B", "C")
  for (r in c(5, 40, 200)) {
    expect_equal(cohalu:::.radius_count_matrix(s$xy, s$cl, lv, r),
                 brute_counts(s$xy, s$cl, lv, r))
  }
})

test_that("cells exactly at the radius and duplicate positions are neighbours", {
  grid <- as.matrix(expand.grid(x = seq(0, 400, by = 100), y = seq(0, 400, by = 100)))
  p <- cohalu:::.radius_pairs(grid, 100)
  expect_equal(length(p$i), nrow(brute_pairs(grid, 100)))
  dup <- rbind(c(0, 0), c(0, 0), c(50, 50))
  p <- cohalu:::.radius_pairs(dup, 1)
  expect_setequal(paste(p$i, p$j), c("1 2", "2 1"))
  expect_equal(p$dist, c(0, 0))
})

test_that("splitting the query into chunks does not change the counts", {
  s <- dense_cells(n = 300)
  lv <- c("A", "B", "C")
  chunks <- split(seq_len(300), ceiling(seq_len(300) / 70))
  parts <- lapply(chunks, function(q) cohalu:::.radius_pairs(s$xy, 25, query = q))
  whole <- cohalu:::.radius_pairs(s$xy, 25)
  expect_setequal(paste(unlist(lapply(parts, `[[`, "i")), unlist(lapply(parts, `[[`, "j"))),
                  paste(whole$i, whole$j))
})

test_that("the ratio keeps changing with the radius beyond the former default cap", {
  s <- dense_cells(n = 1500, side = 300, seed = 2)
  # one corner enriched for A, so that the ratio depends on scale
  corner <- s$xy[, "x"] < 100 & s$xy[, "y"] < 100
  s$cl[corner] <- sample(c("A", "A", "B"), sum(corner), TRUE)
  df <- data.frame(s$xy, cell_type = s$cl, sample_id = "s1", group = "g", patient = "p")
  rownames(df) <- paste0("c", seq_len(nrow(df)))
  ratio_at <- function(r, ...) {
    x <- suppressWarnings(cooccur_ratio_per_sample(
      df, "sample_id", "group", "cell_type", radius = r, ...
    ))
    x$ratio[x$cluster_i == "A" & x$cluster_j == "C"]
  }
  expect_gt(abs(ratio_at(150) - ratio_at(60)), 0.01)
  # the cap reproduces the saturation it used to cause
  expect_lt(abs(ratio_at(150, k = 5) - ratio_at(60, k = 5)), 1e-8)
})

test_that("a numeric cap warns when it is reached and stays silent otherwise", {
  s <- dense_cells()
  lv <- c("A", "B", "C")
  expect_warning(cohalu:::.radius_count_matrix(s$xy, s$cl, lv, 40, k = 5), "k = NULL")
  expect_no_warning(cohalu:::.radius_count_matrix(s$xy, s$cl, lv, 40, k = 600))
  expect_no_warning(cohalu:::.radius_count_matrix(s$xy, s$cl, lv, 0.001, k = 5))
})

test_that("presence of both types uses cells beyond neighbors.k", {
  xy <- cbind(x = c(0, 1, 1.1, 1.2, 10), y = 0)
  labels <- c("Z", "A", "A", "A", "B")
  scores <- cohalu:::.radius_presence(xy, labels, "A", "B", radius = 15)
  expect_equal(scores[1], 1)
  capped <- suppressWarnings(cohalu:::.radius_presence(xy, labels, "A", "B", radius = 15, k = 3))
  expect_equal(capped[1], 0)
})

test_that("cooccur_local_oe counts all neighbours by default", {
  s <- dense_cells()
  d <- data.frame(s$xy, cell_type = s$cl)
  rownames(d) <- paste0("c", seq_len(nrow(d)))
  res <- cooccur_local_oe(d, "A", "B", radius = 60)
  w <- brute_pairs(s$xy, 60)
  expect_equal(res$n_neighbours, as.numeric(tabulate(w[, 1], nrow(d))))
  expect_warning(cooccur_local_oe(d, "A", "B", radius = 60, neighbors.k = 5), "neighbors.k = NULL")
})

test_that("search_interaction_spot keeps cells without neighbours", {
  df <- generate_sim(close_ratio = 0.8, n_types = 4, n_cells = 300,
                     max_loc = 300, test_type = "distribute",
                     distance_param = 10, seed = 1)
  df$sample_id <- "fov1"
  seu <- sim_to_seurat(df)
  target <- c("cell_type_1", "cell_type_2")
  spots <- search_interaction_spot(seu, fov = "fov1", radius = 1e-6, n_min = 0,
                                   cell_id = seu$cell, cluster_col = "cell_type",
                                   target_cluster = target)
  expect_equal(nrow(spots), sum(seu$cell_type %in% target))
  expect_equal(length(unique(spots$cluster_id)), nrow(spots))
})
