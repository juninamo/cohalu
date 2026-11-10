# ---- Local co-localization with a calibrated expectation ------------------
#
# The original sCLS marks a cell with 1 when both cluster_x and cluster_y
# occur within `radius`, then smooths the 0/1 indicator by graph diffusion.
# The indicator has no null: its average grows with the abundance of the two
# cell types and with local cell density, and the diffusion is
# mass-conserving, so it does not change per-sample means at all.
#
# A presence indicator also saturates: when A and B concentrate in one area,
# fewer cells have both nearby than when they are scattered, so its O/E can
# even be negative for strongly co-localized types. cooccur_local_oe()
# therefore counts A-B pairs in each cell's neighbourhood (n_A * n_B),
# divides by the exact expectation under label permutation with positions
# fixed (multivariate hypergeometric, closed form), smooths observed and
# expected counts with a Gaussian kernel of explicit physical width, and
# reports the local log2 observed / expected. Optional permutations give
# per-cell p-values for hotspots. Cost is O(n k) per pass.

# ---- Radius neighbourhoods -------------------------------------------------
#
# A radius search must return every cell within the radius. RANN::nn2(...,
# searchtype = "radius") keeps only the k nearest of them, so with a fixed k
# statistics stop changing once the radius exceeds the k-th neighbour distance.
# dbscan::frNN() has no such cap (the distance test is the same, d <= radius).
# A numeric k still goes through RANN, to reproduce capped results.

# Query cells per frNN() call: bounds the memory used by large radii.
.radius_chunk_size <- 10000L

.radius_chunks <- function(n) split(seq_len(n), ceiling(seq_len(n) / .radius_chunk_size))

# Internal: directed (cell, neighbour) pairs without self, for the rows of xy
# in `query`. `k` = NULL for no cap, otherwise RANN's number of slots, which
# includes the cell itself. n_capped counts query cells whose slots are full.
.radius_pairs <- function(xy, radius, k = NULL, query = seq_len(nrow(xy))) {
  if (is.null(k)) {
    nn <- dbscan::frNN(xy, eps = radius, query = xy[query, , drop = FALSE], sort = FALSE)
    i <- rep(query, lengths(nn$id))
    j <- as.integer(unlist(nn$id, use.names = FALSE))
    d <- as.numeric(unlist(nn$dist, use.names = FALSE))
    n_capped <- 0L
  } else {
    if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 1) {
      stop("The neighbour cap must be NULL or a single positive number.")
    }
    n <- nrow(xy)
    slots <- min(as.integer(k), n)
    res <- RANN::nn2(xy, xy[query, , drop = FALSE], k = slots,
                     searchtype = "radius", radius = radius)
    hit <- res$nn.idx > 0L
    i <- rep(query, slots)[hit]
    j <- res$nn.idx[hit]
    d <- res$nn.dists[hit]
    n_capped <- if (slots < n) sum(rowSums(hit) == slots) else 0L
  }
  keep <- j != i
  list(i = i[keep], j = j[keep], dist = d[keep], n_capped = n_capped)
}

# Internal: as .radius_pairs(), with a cap that excludes the cell itself.
.radius_neighbours <- function(xy, radius, k = NULL) {
  .radius_pairs(xy, radius, if (is.null(k)) NULL else k + 1L)
}

.warn_radius_cap <- function(n_capped, n, k, arg) {
  if (n_capped > 0L) {
    warning(sprintf(paste0("%d of %d cells reached %s = %d within the radius, so their ",
                           "neighbourhoods may be truncated; use %s = NULL for all neighbours."),
                    n_capped, n, arg, as.integer(k), arg), call. = FALSE)
  }
  invisible(NULL)
}

# Internal: cluster x cluster matrix of directed neighbour pairs within radius.
.radius_count_matrix <- function(xy, clusters, all_clusters, radius, k = NULL) {
  n <- nrow(xy)
  L <- length(all_clusters)
  cl <- match(clusters, all_clusters)
  counts <- numeric(L^2)
  n_capped <- 0L
  for (q in .radius_chunks(n)) {
    p <- .radius_pairs(xy, radius, k, q)
    ci <- cl[p$i]
    cj <- cl[p$j]
    ok <- !is.na(ci) & !is.na(cj)
    counts <- counts + tabulate(ci[ok] + (cj[ok] - 1L) * L, nbins = L^2)
    n_capped <- n_capped + p$n_capped
  }
  .warn_radius_cap(n_capped, n, k, "k")
  matrix(counts, L, L,
         dimnames = list(paste0("Cluster", all_clusters), paste0("Cluster", all_clusters)))
}

# Internal: 1 for cells with both cluster_x and cluster_y among the other cells
# within radius, else 0.
.radius_presence <- function(xy, labels, cluster_x, cluster_y, radius, k = NULL) {
  n <- nrow(xy)
  labels <- as.character(labels)
  is_x <- labels %in% cluster_x
  is_y <- labels %in% cluster_y
  has_x <- has_y <- logical(n)
  n_capped <- 0L
  for (q in .radius_chunks(n)) {
    p <- .radius_pairs(xy, radius, k, q)
    has_x[p$i[is_x[p$j]]] <- TRUE
    has_y[p$i[is_y[p$j]]] <- TRUE
    n_capped <- n_capped + p$n_capped
  }
  .warn_radius_cap(n_capped, n, k, "radius_k")
  as.numeric(has_x & has_y)
}

# Internal: numbers of cluster_x / cluster_y cells in each neighbourhood, the
# number of x-y pairs they form, and the expected number of pairs under
# label permutation given each cell's neighbourhood size N:
#   E[n_x n_y]      = N (N - 1) K_x K_y / (M (M - 1))          (x != y)
#   E[n_x (n_x-1)]  = N (N - 1) K_x (K_x - 1) / (M (M - 1))    (x == y)
# with M the other cells and K the cells of each type among them.
.pairs_expected <- function(labels, pairs, cluster_x, cluster_y) {
  n <- length(labels)
  # as.numeric: products of neighbour counts can overflow integers
  count <- function(is) as.numeric(tabulate(pairs$i[is[pairs$j]], nbins = n))
  isx <- labels == cluster_x
  N <- as.numeric(tabulate(pairs$i, nbins = n))
  M <- n - 1
  kx <- sum(isx) - isx
  nx <- count(isx)
  if (identical(cluster_x, cluster_y)) {
    ny <- nx
    obs <- nx * (nx - 1)
    e <- N * (N - 1) * kx * (kx - 1) / (M * (M - 1))
  } else {
    isy <- labels == cluster_y
    ky <- sum(isy) - isy
    ny <- count(isy)
    obs <- nx * ny
    e <- N * (N - 1) * kx * ky / (M * (M - 1))
  }
  list(n_x = nx, n_y = ny, pairs = obs, expected = e, n_neighbours = N)
}

#' Local co-localization score with a permutation-calibrated expectation
#'
#' An abundance-adjusted replacement for the diffusion-smoothed sCLS of
#' [cooccur_local()]. For every cell \eqn{v}, the number of
#' `cluster_x`-`cluster_y` pairs among the cells within `radius`,
#' \eqn{s_v = n_x(v)\, n_y(v)}, is compared with its exact expectation
#' under random relabelling of cells with positions fixed,
#' \deqn{E_v = N_v (N_v - 1) \frac{K_x K_y}{M (M - 1)},}
#' where \eqn{N_v} is the number of neighbours of \eqn{v}, \eqn{M} the
#' number of other cells and \eqn{K_x, K_y} the numbers of cells of each
#' type among them. Observed and expected values are summed with Gaussian
#' weights of width `bandwidth` (micrometres) around each cell, giving the
#' local
#' \deqn{\log_2 \mathrm{O/E}_v = \log_2 \frac{\sum_u w_{vu} s_u + c}{\sum_u w_{vu} E_u + c}.}
#' Summed over the whole section this is the section-level O/E, which, unlike
#' the mean sCLS, does not grow with the abundance of the two cell types or
#' with cell density. Pair counts are used rather than the 0/1 "both
#' present" indicator of the original sCLS, because the indicator saturates:
#' strongly co-localized types occupy fewer neighbourhoods than scattered
#' ones, which would give them a negative O/E.
#'
#' @param df A data.frame with `x`, `y` and a cell-type column; row names
#'   are used as cell IDs.
#' @param cluster_x,cluster_y The two cell types (may be identical).
#' @param radius Neighbourhood radius for the indicator.
#' @param neighbors.k Optional cap on the neighbours used per cell, excluding
#'   the cell itself. `NULL` (default) uses every cell within `radius`. A
#'   number truncates dense neighbourhoods to the nearest cells, and a warning
#'   reports how many cells reached it.
#' @param bandwidth Width of the Gaussian smoothing window (the kernel SD is
#'   `bandwidth / 2`, truncated at `bandwidth`). Defaults to `radius`.
#' @param cluster_col Column holding the cell types.
#' @param n_perms Number of label permutations for per-cell hotspot p-values
#'   (0 = none). Each permutation costs one pass over the neighbour index.
#' @param pseudocount Added to smoothed observed and expected pair counts.
#' @param seed Random seed for the permutations.
#'
#' @return A data.frame (one row per cell) with `n_x`, `n_y`, `pairs`,
#'   `expected`, `n_neighbours`, `local_log2_oe` and, if `n_perms > 0`,
#'   one-sided `p` (enrichment) and BH-adjusted `padj`. The attribute
#'   `section_log2_oe` holds log2(total pairs / total expected pairs).
#' @export
#' @examples
#' df <- generate_sim(close_ratio = 1, n_types = 6, n_cells = 600,
#'                    max_loc = 400, test_type = "circle",
#'                    distance_param = 15, seed = 3)
#' rownames(df) <- paste0("cell", seq_len(nrow(df)))
#' res <- cooccur_local_oe(df, "cell_type_1", "cell_type_2", radius = 25,
#'                         n_perms = 99)
#' attr(res, "section_log2_oe")
#' head(res)
cooccur_local_oe <- function(df, cluster_x, cluster_y, radius = 30, neighbors.k = NULL,
                             bandwidth = radius, cluster_col = "cell_type",
                             n_perms = 0, pseudocount = 0.5, seed = 1) {
  if (!all(c("x", "y", cluster_col) %in% colnames(df))) {
    stop("`df` needs columns x, y and ", cluster_col)
  }
  labels <- as.character(df[[cluster_col]])
  if (!any(labels == cluster_x) || !any(labels == cluster_y)) {
    stop("cluster_x or cluster_y not present in `df`.")
  }
  xy <- as.matrix(df[, c("x", "y")])
  n <- nrow(xy)
  nb <- .radius_neighbours(xy, radius, neighbors.k)
  ie <- .pairs_expected(labels, nb, cluster_x, cluster_y)

  # Gaussian smoothing weights (self included), sigma = bandwidth / 2
  sm <- if (isTRUE(all.equal(bandwidth, radius))) nb else .radius_neighbours(xy, bandwidth, neighbors.k)
  .warn_radius_cap(max(nb$n_capped, sm$n_capped), n, neighbors.k, "neighbors.k")
  sig <- bandwidth / 2
  w <- exp(-sm$dist^2 / (2 * sig^2))
  W <- Matrix::sparseMatrix(i = c(sm$i, seq_len(n)),
                            j = c(sm$j, seq_len(n)),
                            x = c(w, rep(1, n)), dims = c(n, n))
  obs_s <- as.numeric(W %*% ie$pairs)
  exp_s <- as.numeric(W %*% ie$expected)
  out <- data.frame(n_x = ie$n_x, n_y = ie$n_y, pairs = ie$pairs, expected = ie$expected,
                    n_neighbours = ie$n_neighbours,
                    local_log2_oe = log2((obs_s + pseudocount) / (exp_s + pseudocount)),
                    row.names = rownames(df))
  if (n_perms > 0) {
    .local_seed(seed)
    exceed <- numeric(n)
    for (r in seq_len(n_perms)) {
      s_perm <- .pairs_expected(sample(labels), nb, cluster_x, cluster_y)$pairs
      exceed <- exceed + (as.numeric(W %*% s_perm) >= obs_s - 1e-12)
    }
    out$p <- (exceed + 1) / (n_perms + 1)
    out$padj <- stats::p.adjust(out$p, method = "BH")
  }
  attr(out, "section_log2_oe") <- log2(sum(ie$pairs) / max(sum(ie$expected), .Machine$double.eps))
  out
}
