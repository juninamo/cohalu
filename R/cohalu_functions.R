# ---- Dependencies ----
# Only the functions called without a namespace prefix are imported; all
# other external calls use `pkg::fun()`. Importing whole namespaces caused
# "replacing previous import" conflicts (dplyr / igraph / tibble / Matrix).
#' @importFrom magrittr %>%
#' @importFrom dplyr n sym
#' @importFrom Matrix bdiag
#' @importFrom Seurat FindNeighbors Cells
NULL

# ---- Main Functions ----

utils::globalVariables(c(".", "cell", "cell_type", "n_all_cells", "rnorm", "runif", "sd", "x", "y"))

#' @name generate_sim
#' @title Simulate Spatial Coordinates and Cell Types
#'
#' @description Generate synthetic spatial transcriptomics data for simulation and benchmarking.
#'
#' @param close_ratio Proportion of close interactions between selected cell types.
#' @param n_types Number of distinct cell types.
#' @param max_loc Maximum coordinate value (spatial extent).
#' @param n_perm Number of permutations to simulate, for use in future analysis.
#' @param n_cells Total number of cells to simulate.
#' @param test_type Type of spatial pattern to simulate. One of "circle", "line", or "distribute".
#' @param distance_param Distance parameter controlling interaction distance.
#' @param seed Random seed for reproducibility.
#'
#' @return A data.frame with simulated spatial coordinates and cell type labels.
#' @export
#' @examples
#' df <- generate_sim(close_ratio = 0.8, n_types = 4, n_cells = 300,
#'                    max_loc = 300, test_type = "distribute",
#'                    distance_param = 10, seed = 1)
#' head(df)
#' table(df$cell_type)
generate_sim <- function(close_ratio = 0.7,
                         n_types = 10,
                         max_loc = 800,
                         n_perm = 100,
                         n_cells = 1500,
                         test_type = "circle",
                         distance_param = 50,
                         seed = 1234) {
  if(test_type=="distribute"){

    .seed_rng(seed)

    x_coords <- runif(n_cells, min = 0, max = max_loc)
    y_coords <- runif(n_cells, min = 0, max = max_loc)

    cell_types <- sample(paste0("cell_type_", seq_len(n_types)), n_cells, replace = TRUE)

    idx_type_1 <- which(cell_types == "cell_type_1")
    idx_type_2 <- which(cell_types == "cell_type_2")

    n_pairs <- round(length(idx_type_2) * close_ratio)
    idx_close_1 <- sample(idx_type_1, n_pairs, replace = TRUE)
    idx_close_2 <- sample(idx_type_2, n_pairs)

    angle_shift <- runif(n_pairs, 0, 2 * pi)
    x_coords[idx_close_2] <- x_coords[idx_close_1] + distance_param * cos(angle_shift) + rnorm(n_pairs, mean = 0, sd = distance_param/5)
    y_coords[idx_close_2] <- y_coords[idx_close_1] + distance_param * sin(angle_shift) + rnorm(n_pairs, mean = 0, sd = distance_param/5)
    # Keep relocated cells inside the tissue: a cell placed outside the square
    # sits in an empty region where its k nearest neighbours reach far away
    # (up to the partner cell), which leaks co-localization to long distances.
    # Out-of-bounds cells get a new partner and direction (rejection sampling).
    outside <- function(i) x_coords[i] < 0 | x_coords[i] > max_loc | y_coords[i] < 0 | y_coords[i] > max_loc
    for (iter in seq_len(100)) {
      bad <- which(outside(idx_close_2))
      if (!length(bad)) break
      idx_close_1[bad] <- idx_type_1[sample.int(length(idx_type_1), length(bad), replace = TRUE)]
      a <- runif(length(bad), 0, 2 * pi)
      x_coords[idx_close_2[bad]] <- x_coords[idx_close_1[bad]] + distance_param * cos(a) + rnorm(length(bad), 0, distance_param / 5)
      y_coords[idx_close_2[bad]] <- y_coords[idx_close_1[bad]] + distance_param * sin(a) + rnorm(length(bad), 0, distance_param / 5)
    }

    df <- data.frame(x = x_coords, y = y_coords, cell_type = cell_types) %>%
      dplyr::mutate(cell_type = factor(cell_type, levels = paste0("cell_type_", seq_len(n_types))))

  } else if(test_type=="circle"){

    .seed_rng(seed)

    n1_total <- n2_total <- ceiling(n_cells/n_types)
    n_other <- n_cells-n1_total-n2_total

    close_ratio_1 <- close_ratio_2 <- close_ratio

    center_x <- ceiling(max_loc/2)
    center_y <- ceiling(max_loc/2)
    radius <- ceiling(center_x/3)

    n1_close <- round(n1_total * close_ratio_1)

    theta_1_outer <- seq(0, 2*pi, length.out = round((n1_close / 4)))
    x_1_outer <- center_x + radius * cos(theta_1_outer)# + rnorm(length(theta_1_outer), mean = distance_param, sd = ceiling(distance_param/10))
    y_1_outer <- center_y + radius * sin(theta_1_outer)# + rnorm(length(theta_1_outer), mean = distance_param, sd = ceiling(distance_param/10))

    r_1_inner <- sqrt(runif(round(n1_close - round((n1_close / 4))), 0, radius^2))
    theta_1_inner <- runif(round(n1_close - round((n1_close / 4))), 0, 2*pi)
    x_1_inner <- center_x + r_1_inner * cos(theta_1_inner)
    y_1_inner <- center_y + r_1_inner * sin(theta_1_inner)

    n1_random <- n1_total - length(x_1_inner) - length(x_1_outer)
    x_1_random <- runif(n1_random, min = 0, max = max_loc)
    y_1_random <- runif(n1_random, min = 0, max = max_loc)

    x_1 <- c(x_1_inner, x_1_outer, x_1_random)
    y_1 <- c(y_1_inner, y_1_outer, y_1_random)

    n2_close <- round(n2_total * close_ratio_2)

    theta_1_outer <- seq(0, 2*pi, length.out = round(n2_close))
    x_2_close <- center_x + (radius+distance_param) * cos(theta_1_outer) + rnorm(n2_close, mean = 0, sd = 5)
    y_2_close <- center_y + (radius+distance_param) * sin(theta_1_outer) + rnorm(n2_close, mean = 0, sd = 5)

    n2_random <- n2_total - n2_close
    x_2_random <- runif(n2_random, min = 0, max = max_loc)
    y_2_random <- runif(n2_random, min = 0, max = max_loc)

    x_other <- runif(n_other, min = 0, max = max_loc)
    y_other <- runif(n_other, min = 0, max = max_loc)

    df <- data.frame(
      x = c(x_1, x_2_close, x_2_random, x_other),
      y = c(y_1, y_2_close, y_2_random, y_other),
      cell_type = c(rep("cell_type_1", n1_total),
                    rep("cell_type_2", n2_close),
                    rep("cell_type_2", n2_random),
                    sample(paste0("cell_type_", 3:n_types), n_other, replace = TRUE))
    ) %>%
      dplyr::mutate(cell_type = factor(cell_type, levels = paste0("cell_type_", seq_len(n_types))))

  } else if (test_type=="line"){

    .seed_rng(seed)

    n_cells_per_layer <- ceiling(n_cells/n_types * max(close_ratio,0.01))
    n_layers <- 2
    layer_spacing <- distance_param
    x_start <- ceiling(max_loc/2)-ceiling(max_loc/4)
    y_start <- ceiling(max_loc/2)-50
    x_end <- ceiling(max_loc/2)+ceiling(max_loc/4)
    step <- (x_end-x_start)/n_cells_per_layer

    x_1 <- c()
    y_1 <- c()
    x_2 <- c()
    y_2 <- c()

    for (i in 0:(n_layers - 1)) {
      y_pos <- y_start + i * layer_spacing

      if (i %% 2 == 0) {
        x_layer <- seq(x_start, x_start + n_cells_per_layer * step, length.out = n_cells_per_layer) + rnorm(n_cells_per_layer, mean = 0, sd = 1)
        y_layer <- rep(y_pos, n_cells_per_layer) + rnorm(n_cells_per_layer, 0, 2)

        n_close <- round(n_cells_per_layer * close_ratio)
        x_1 <- c(x_1, x_layer)
        y_1 <- c(y_1, y_layer)
      } else {
        x_layer <- seq(x_start, x_start + n_cells_per_layer * step, length.out = n_cells_per_layer) + rnorm(n_cells_per_layer, mean = 0, sd = 1)
        y_layer <- rep(y_pos, n_cells_per_layer) + rnorm(n_cells_per_layer, 0, 2)

        n_close <- round(n_cells_per_layer * close_ratio)
        idx_close <- sample(seq_len(n_cells_per_layer), n_close, replace = TRUE)
        angle_shift <- runif(n_close, 0, 2 * pi)
        x_layer[idx_close] <- x_1[idx_close]
        y_layer[idx_close] <- y_1[idx_close] + distance_param + rnorm(n_close, mean = 0, sd = 1)

        x_2 <- c(x_2, x_layer)
        y_2 <- c(y_2, y_layer)
      }
    }

    n1_random <- n2_random <- ceiling(ceiling(n_cells/n_types) * (1-close_ratio))
    x_1_random <- runif(n1_random, min = 0, max = max_loc)
    y_1_random <- runif(n1_random, min = 0, max = max_loc)
    x_2_random <- runif(n2_random, min = 0, max = max_loc)
    y_2_random <- runif(n2_random, min = 0, max = max_loc)
    n_other <- n_cells - n1_random - n2_random
    x_other <- runif(n_other, min = 0, max = max_loc)
    y_other <- runif(n_other, min = 0, max = max_loc)

    df <- data.frame(
      x = c(x_1, x_1_random, x_2, x_2_random, x_other),
      y = c(y_1, y_1_random, y_2, y_2_random, y_other),
      cell_type = c(rep("cell_type_1", sum(length(x_1),length(x_1_random))),
                    rep("cell_type_2", sum(length(x_2),length(x_2_random))),
                    sample(paste0("cell_type_", 3:n_types), n_other, replace = TRUE))
    ) %>%
      dplyr::mutate(cell_type = factor(cell_type, levels = paste0("cell_type_", seq_len(n_types))))

  }
  return(df)
}

#' Compute Co-occurrence Count Matrix
#'
#' @param adj Adjacency matrix.
#' @param int_clust_row Vector of cluster labels for rows.
#' @param int_clust_col Vector of cluster labels for columns.
#' @param n_cls Number of clusters.
#' @param cluster_data Original cluster assignments.
#' @param transformation If `TRUE`, `adj` holds link weights (e.g. 1 / (1 + d_u)) that are summed; `FALSE` (default) counts binary links.
#'
#' @return A co-occurrence count matrix.
#' @export
#' @examples
#' set.seed(1)
#' adj <- Matrix::rsparsematrix(20, 20, density = 0.2)
#' cl <- factor(sample(c("a", "b"), 20, replace = TRUE))
#' lab <- paste0("Cluster", cl)
#' compute_count(adj, lab, lab, n_cls = 2, cluster_data = cl)
compute_count <- function(adj, int_clust_row, int_clust_col, n_cls, cluster_data, transformation = FALSE) {
  lv <- paste0("Cluster", levels(cluster_data))
  # One-hot cluster indicator matrices: counts = t(M_row) %*% adj %*% M_col,
  # i.e. counts[i, j] = sum of adj over (row cells in i) x (col cells in j).
  indicator <- function(lab) {
    j <- match(lab, lv)
    keep <- !is.na(j)
    Matrix::sparseMatrix(i = which(keep), j = j[keep], x = 1,
                         dims = c(length(lab), length(lv)))
  }
  a <- if (transformation) adj else (adj == 1) * 1
  counts <- as.matrix(Matrix::t(indicator(int_clust_row)) %*% a %*% indicator(int_clust_col))
  dimnames(counts) <- list(lv, lv)
  counts
}

#' Permute Cluster Assignments and Recompute Counts
#'
#' @param adj Adjacency matrix.
#' @param int_clust Cluster labels.
#' @param n_cls Number of clusters.
#' @param cluster_data Original cluster data.
#' @param transformation Whether to apply adjacency transformation.
#'
#' Cluster labels are shuffled once and the same permutation is applied to
#' rows and columns of the adjacency matrix.
#'
#' @return Permuted co-occurrence count matrix.
#' @export
#' @examples
#' set.seed(1)
#' adj <- Matrix::rsparsematrix(20, 20, density = 0.2)
#' cl <- factor(sample(c("a", "b"), 20, replace = TRUE))
#' permute_clusters(adj, paste0("Cluster", cl), n_cls = 2, cluster_data = cl,
#'                  transformation = TRUE)
permute_clusters <- function(adj, int_clust, n_cls, cluster_data, transformation) {
  # The same shuffled labels must be used for rows and columns: each cell keeps
  # a single (random) label, so self-loops and mutual neighbours are handled
  # identically in the observed and permuted counts. Shuffling rows and
  # columns independently inflates same-type z-scores under the null.
  int_clust_perm <- sample(int_clust)
  compute_count(adj, int_clust_row = int_clust_perm, int_clust_col = int_clust_perm, n_cls, cluster_data, transformation)
}

#' Calculate Co-occurrence Matrix for a Given Radius
#'
#' @param seurat_obj Seurat object with spatial coordinates.
#' @param radius Radius to define local neighborhood.
#' @param sample_key Metadata column specifying sample identity.
#' @param cluster_key Metadata column specifying cluster labels.
#' @param k Optional cap on the cells returned per cell, counting the cell
#'   itself. `NULL` (default) uses every cell within `radius`. A number keeps
#'   only the `k` nearest cells, so the result stops changing with `radius`
#'   once the cap is reached; a warning reports how many cells reached it.
#'
#' @return A list with co-occurrence count and enrichment ratio matrices.
#' @export
#' @examples
#' df <- generate_sim(close_ratio = 0.8, n_types = 4, n_cells = 300,
#'                    max_loc = 300, test_type = "distribute",
#'                    distance_param = 10, seed = 1)
#' df$sample_id <- "fov1"
#' seu <- sim_to_seurat(df)
#' res <- calc_co_occurrence_for_radius(seu, radius = 20,
#'                                      sample_key = "sample_id",
#'                                      cluster_key = "cell_type")
#' round(res$ratio_mat, 2)
calc_co_occurrence_for_radius <- function(seurat_obj, radius, sample_key, cluster_key, k = NULL) {
  all_clusters <- levels(factor(seurat_obj@meta.data[[cluster_key]]))

  co_occur_count <- matrix(
    0,
    nrow = length(all_clusters),
    ncol = length(all_clusters),
    dimnames = list(
      paste0("Cluster", all_clusters),
      paste0("Cluster", all_clusters)
    )
  )

  for (name in names(seurat_obj@images)) {
    coords <- seurat_obj[[name]]$centroids@coords %>%
      as.data.frame() %>%
      dplyr::mutate(cell = Cells(seurat_obj[[name]])) %>%
      tibble::column_to_rownames(var = "cell") %>%
      dplyr::mutate(
        cluster = seurat_obj@meta.data[
          seurat_obj@meta.data[, sample_key] == name,
          cluster_key
        ]
      )

    co_occur_count <- co_occur_count +
      .radius_count_matrix(as.matrix(coords[, c("x", "y")]), coords$cluster,
                           all_clusters, radius, k)
  }

  ratio_mat <- compute_co_occurrence_ratio(co_occur_count)

  return(list(
    co_occur_count = co_occur_count,
    ratio_mat      = ratio_mat
  ))
}

#' Compute Enrichment Ratios from Count Matrix
#'
#' @param co_occur_count Matrix of observed co-occurrence counts.
#'
#' @return A matrix of normalized enrichment ratios.
#' @export
#' @examples
#' counts <- matrix(c(10, 2, 2, 6), 2,
#'                  dimnames = list(c("A", "B"), c("A", "B")))
#' compute_co_occurrence_ratio(counts)
compute_co_occurrence_ratio <- function(co_occur_count) {
  rn <- rownames(co_occur_count)
  cn <- colnames(co_occur_count)
  row_sums <- rowSums(co_occur_count)
  col_sums <- colSums(co_occur_count)
  total_sum <- sum(co_occur_count)

  ratio_mat <- matrix(
    0,
    nrow = nrow(co_occur_count),
    ncol = ncol(co_occur_count),
    dimnames = list(rn, cn)
  )

  for (i in seq_len(nrow(co_occur_count))) {
    for (j in seq_len(ncol(co_occur_count))) {
      p_exp_cond <- co_occur_count[i, j] / row_sums[i]   # p(exp=j | cond=i)
      p_exp      <- col_sums[j]       / total_sum        # p(exp=j)
      ratio_mat[i, j] <- p_exp_cond / p_exp
    }
  }
  return(ratio_mat)
}

#' Search for Spatial Interaction Spots
#'
#' @param seurat_object Seurat object.
#' @param fov Field of view identifier.
#' @param radius Radius threshold for neighborhood.
#' @param n_min Minimum number of cells to qualify as a spot.
#' @param neighbors.k Optional cap on the cells returned per cell, counting
#'   the cell itself. `NULL` (default) uses every cell within `radius`.
#' @param cell_id Vector of target cell IDs.
#' @param cluster_col Column name in metadata specifying cluster assignment.
#' @param target_cluster Target cluster(s) assumed to be interacting.
#'
#' @return A data.frame of detected interaction clusters and metadata.
#' @export
#' @examples
#' df <- generate_sim(close_ratio = 0.8, n_types = 4, n_cells = 300,
#'                    max_loc = 300, test_type = "distribute",
#'                    distance_param = 10, seed = 1)
#' df$sample_id <- "fov1"
#' seu <- sim_to_seurat(df)
#' spots <- search_interaction_spot(seu, fov = "fov1", radius = 15, n_min = 3,
#'                                  cell_id = seu$cell,
#'                                  cluster_col = "cell_type",
#'                                  target_cluster = c("cell_type_1", "cell_type_2"))
#' length(unique(spots$cluster_id))
search_interaction_spot <- function(seurat_object, fov, radius, n_min, neighbors.k = NULL, cell_id = cell_id, cluster_col = cluster_col, target_cluster = target_cluster) {
  coords <- seurat_object[[fov]]$centroids@coords %>%
    as.data.frame() %>%
    dplyr::mutate(cell = Cells(seurat_object[[fov]])) %>%
    dplyr::filter(cell %in% cell_id)
  cells <- coords$cell
  rownames(coords) <- cells
  coords <- as.matrix(coords[, c("x", "y")])
  dim(coords)

  pairs <- .radius_pairs(coords, radius, neighbors.k)
  .warn_radius_cap(pairs$n_capped, nrow(coords), neighbors.k, "neighbors.k")

  # (1) Create a ‘cell ID column’ in coords (use the row name if there is one)
  coords_df_ <- coords %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "cell_id")

  # (2) Link each cell to the cells within the radius
  edges_vec <- c(rbind(pairs$i, pairs$j))

  # (3) Creating an igraph object; n keeps cells without neighbours as vertices
  # Assuming an undirected graph (directed=FALSE)
  g <- igraph::graph(edges = edges_vec, n = nrow(coords), directed = FALSE)

  # (4) Search for a connected component to obtain a cluster ID
  comp <- igraph::components(g, mode = "weak")
  # component numbers (1, 2, 3, ...) to which each vertex (cell) belongs is stored in comp$membership.

  cluster_id <- comp$membership

  # (5) Add a cluster_id column to coords_df
  coords_df_$cluster_id <- cluster_id

  coords_df_ <- coords_df_ %>%
    dplyr::left_join(.,seurat_object@meta.data[,c("cell",cluster_col)] %>% dplyr::rename(cell_id = cell),by="cell_id") %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::mutate(x_min = min(x),
                  x_max = max(x),
                  y_min = min(y),
                  y_max = max(y),
                  n_all_cells = n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!!sym(cluster_col) %in% target_cluster) %>%
    dplyr::filter(n_all_cells > n_min) %>%
    dplyr::mutate(sample_id = fov,
                  cluster_id = paste0(fov,":cluster",cluster_id))

  return(coords_df_)
}

# Internal: spatial kNN graph. "nn": exact k nearest neighbours by a kd-tree
# (RANN), each cell counting as its own neighbour (k includes the cell, as in
# Seurat::FindNeighbors); A[u, v] = 1 if v is among the k nearest cells of u.
# "snn": Seurat's shared-nearest-neighbour graph. Seurat's default kNN search
# (annoy) is approximate and was ~80x slower than the kd-tree in 2-D
# (200,000 cells: 33 s vs 0.4 s, 99.995% identical links).
.knn_graph <- function(coords, k, connectivity_key = "nn") {
  coords <- as.matrix(coords)
  if (identical(connectivity_key, "snn")) {
    return(FindNeighbors(coords, k.param = k, verbose = FALSE)$snn)
  }
  n <- nrow(coords)
  k <- min(k, n)
  idx <- RANN::nn2(coords, k = k)$nn.idx
  Matrix::sparseMatrix(i = rep(seq_len(n), k), j = as.vector(idx), x = 1, dims = c(n, n),
                       dimnames = list(rownames(coords), rownames(coords)))
}

# Internal: permutation null shared by nhood_enrichment() and its Seurat
# method. Returns observed counts, permutation mean ("expected"), z-score and
# log2(observed / expected).
# Internal: directional neighbourhood statistics (row = centre cell type i,
# column = neighbour cell type j; not symmetric), from the unweighted kNN graph
# without the cell itself:
#   contact[i, j]   = share of type-i cells with >= 1 type-j neighbour
#   dominance[i, j] = share of type-i cells whose neighbours are at least half type j
.directional <- function(adj, int_clust, lv) {
  Ab <- Matrix::drop0((adj != 0) * 1)
  Matrix::diag(Ab) <- 0
  deg <- Matrix::rowSums(Ab)
  j <- match(int_clust, lv); keep <- !is.na(j)
  M <- Matrix::sparseMatrix(i = which(keep), j = j[keep], x = 1, dims = c(length(int_clust), length(lv)))
  nb <- as.matrix(Ab %*% M)                                   # neighbours of each type, per cell
  n_i <- Matrix::colSums(M)
  agg <- function(hit) { P <- as.matrix(Matrix::t(M) %*% (hit * 1)) / n_i; P[n_i == 0, ] <- NA_real_; dimnames(P) <- list(lv, lv); P }
  list(contact = agg(nb >= 1), dominance = agg(nb >= pmax(1, deg / 2)))
}

# Internal: one label shuffle -> contact counts and directional statistics from the same labels.
.permute_both <- function(adj, int_clust, n_cls, cluster_data, transformation, lv) {
  perm <- sample(int_clust)
  c(list(count = compute_count(adj, int_clust_row = perm, int_clust_col = perm, n_cls, cluster_data, transformation)),
    .directional(adj, perm, lv))
}

.nhood_permutation_core <- function(adj, cluster_data, transformation, n_perms, seed, n_jobs) {
  # Keep pre-set factor levels (e.g. harmonized across samples) so that absent
  # cell types yield NA rather than silently dropping rows / columns.
  if (!is.factor(cluster_data)) cluster_data <- factor(cluster_data)
  int_clust <- paste0("Cluster", cluster_data)
  n_cls <- length(levels(cluster_data))

  count <- compute_count(adj, int_clust_row = int_clust, int_clust_col = int_clust,
                         n_cls, cluster_data, transformation)
  lv <- paste0("Cluster", levels(cluster_data))
  dir_obs <- .directional(adj, int_clust, lv)

  .local_seed(seed)
  run_seq <- function() {
    .seed_rng(seed)
    lapply(seq_len(n_perms), function(x) .permute_both(adj, int_clust, n_cls, cluster_data, transformation, lv))
  }
  perms <- if (n_jobs <= 1L) {
    run_seq()
  } else {
    cl <- parallel::makeCluster(n_jobs)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    # Seed the worker RNG streams so results are reproducible for a given seed.
    parallel::clusterSetRNGStream(cl, iseed = seed)
    # The workers need the Matrix methods for the sparse adjacency (==, !=,
    # %*%); without them every worker fails and the run falls back to
    # sequential, silently losing the speed-up.
    parallel::clusterEvalQ(cl, suppressPackageStartupMessages(requireNamespace("Matrix", quietly = TRUE)))
    # Ship self-contained copies of the workers' functions so that the
    # workers do not need cohalu itself to be installed.
    fn_env <- new.env(parent = baseenv())
    fn_env$compute_count <- compute_count
    fn_env$.directional <- .directional
    fn_env$.permute_both <- .permute_both
    environment(fn_env$compute_count) <- fn_env
    environment(fn_env$.directional) <- fn_env
    environment(fn_env$.permute_both) <- fn_env
    worker <- function(x) fn_env$.permute_both(adj, int_clust, n_cls, cluster_data, transformation, lv)
    environment(worker) <- list2env(
      list(fn_env = fn_env, adj = adj, int_clust = int_clust, n_cls = n_cls,
           cluster_data = cluster_data, transformation = transformation, lv = lv),
      parent = baseenv()
    )
    out <- NULL
    for (attempt in seq_len(3L)) {
      out <- tryCatch(
        parallel::parLapply(cl, seq_len(n_perms), worker),
        error = function(e) { message("[cohalu] parallel attempt ", attempt, " failed: ", conditionMessage(e)); NULL })
      if (!is.null(out)) break
      Sys.sleep(1)
    }
    if (is.null(out)) {
      # bounded retries exhausted -> sequential fallback (avoids infinite hang)
      warning("[cohalu] parallel permutations failed; falling back to sequential (n_jobs = 1).", call. = FALSE)
      out <- run_seq()
    }
    out
  }

  dir_perm <- lapply(c(contact = "contact", dominance = "dominance"), function(nm) lapply(perms, `[[`, nm))
  perms <- lapply(perms, `[[`, "count")
  arr <- simplify2array(perms)
  perm_mean <- apply(arr, c(1, 2), mean)
  perm_sd <- apply(arr, c(1, 2), sd)
  zscore <- (count - perm_mean) / perm_sd
  zscore[!is.finite(zscore)] <- NA_real_
  dimnames(perm_mean) <- dimnames(zscore) <- dimnames(count)

  # log2(observed / expected) with a pseudocount of one average edge weight,
  # so the effect size is on the same scale with or without transformation.
  nz <- adj@x[adj@x != 0]
  pc <- if (length(nz)) mean(nz) else 1
  lr <- function(x) log2((x + pc) / (perm_mean + pc))
  log2_oe_raw <- lr(count)
  # The log of a ratio of small counts is biased downwards (Jensen), so rare
  # cell types get log2 O/E slightly below 0 even without any interaction.
  # Subtracting the mean of the same statistic over the label shuffles
  # centres it at exactly 0 under the null for any number of cells.
  null_bias <- apply(simplify2array(lapply(perms, lr)), c(1, 2), mean)
  log2_oe <- log2_oe_raw - null_bias
  log2_oe[perm_mean == 0 & count == 0] <- NA_real_
  log2_oe_raw[perm_mean == 0 & count == 0] <- NA_real_

  # Within-sample test per unordered pair: contacts i->j and j->i are summed
  # (the degree-normalised counts are not exactly symmetric). The observed
  # sum is standardised together with the shuffles (mean and SD over all
  # B + 1 values, so the observed tissue is treated exactly like a shuffle).
  #  * pvalue: two-sided normal p-value, for a single pre-specified pair.
  #  * padj: Westfall-Young single-step max-T over the K (K + 1) / 2 pairs,
  #    i.e. the share of shuffles whose largest |z| over all pairs reaches
  #    the observed |z|. It controls the family-wise error rate and adapts
  #    to the skewed, dependent null of rare pairs; BH on normal p-values
  #    (padj_bh) can exceed its level when there are many rare cell types.
  sym <- function(m) m + t(m) - diag(diag(m), nrow(m))
  ut <- upper.tri(count, diag = TRUE)
  S <- sym(count)[ut]
  SA <- vapply(perms, function(m) sym(m)[ut], numeric(sum(ut)))
  if (is.null(dim(SA))) SA <- matrix(SA, nrow = 1)
  ALL <- cbind(S, SA)
  mu <- rowMeans(ALL)
  sdv <- sqrt(rowSums((ALL - mu)^2) / (ncol(ALL) - 1))
  z_obs <- abs(S - mu) / sdv
  z_perm <- abs(SA - mu) / sdv
  z_obs[!is.finite(z_obs)] <- NA_real_
  z_perm[!is.finite(z_perm)] <- 0
  max_perm <- apply(z_perm, 2, max)
  p_maxt <- vapply(z_obs, function(q) if (is.na(q)) NA_real_ else (1 + sum(max_perm >= q)) / (length(max_perm) + 1), 0)
  to_mat <- function(v) {
    m <- matrix(NA_real_, nrow(count), ncol(count))
    m[ut] <- v
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    m
  }
  p_norm <- 2 * stats::pnorm(-z_obs)
  pvalue <- to_mat(p_norm)
  padj <- to_mat(p_maxt)
  padj_bh <- to_mat(stats::p.adjust(p_norm, method = "BH"))
  dimnames(log2_oe) <- dimnames(log2_oe_raw) <- dimnames(pvalue) <- dimnames(padj) <- dimnames(padj_bh) <- dimnames(count)

  # Directional statistics (row = centre cell type, column = neighbour type),
  # compared with the same label shuffles; log2 O/E centred like log2_oe; one
  # test per ordered pair, max-T over all K x K ordered pairs.
  dir_stats <- function(obs, parts, name) {
    PA <- simplify2array(parts); ex <- apply(PA, c(1, 2), mean); ppc <- 0.01
    plr <- function(x) log2((x + ppc) / (ex + ppc))
    lo <- plr(obs) - apply(simplify2array(lapply(parts, plr)), c(1, 2), mean)
    lo[obs == 0 & ex == 0] <- NA_real_                  # not attainable, observed or by chance
    PS <- as.vector(obs); PP <- matrix(PA, ncol = dim(PA)[3])
    ALLd <- cbind(PS, PP); m_ <- rowMeans(ALLd); s_ <- sqrt(rowSums((ALLd - m_)^2) / (ncol(ALLd) - 1))
    zo <- abs(PS - m_) / s_; zp <- abs(PP - m_) / s_; zo[!is.finite(zo)] <- NA_real_; zp[!is.finite(zp)] <- 0
    mx <- apply(zp, 2, max)
    pa <- vapply(zo, function(q) if (is.na(q)) NA_real_ else (1 + sum(mx >= q)) / (length(mx) + 1), 0)
    shape <- function(v) { m <- matrix(v, nrow(count), ncol(count)); dimnames(m) <- dimnames(count); m }
    dimnames(obs) <- dimnames(ex) <- dimnames(lo) <- dimnames(count)
    out <- list(obs, ex, lo, shape(2 * stats::pnorm(-zo)), shape(pa))
    names(out) <- paste0(name, c("", "_expected", "_log2_oe", "_pvalue", "_padj")); out
  }
  dir_out <- c(dir_stats(dir_obs$contact, dir_perm$contact, "contact"), dir_stats(dir_obs$dominance, dir_perm$dominance, "dominance"))

  c(list(zscore = zscore, count = count, expected = perm_mean, log2_oe = log2_oe,
       log2_oe_raw = log2_oe_raw, pvalue = pvalue, padj = padj, padj_bh = padj_bh), dir_out)
}

#' Neighborhood Enrichment (Seurat Method)
#'
#' @param seurat_obj A Seurat object with spatial coordinates.
#' @param cluster_key Metadata column for cluster IDs.
#' @param neighbors.k Number of neighbors to construct graph.
#' @param connectivity_key Which graph to use: "nn" or "snn".
#' @param transformation If `TRUE`, each link from cell u is weighted
#'   1 / (1 + d_u), with d_u the number of cells that chose u as a neighbour
#'   (a mild down-weighting of hub cells). `FALSE` (default since 0.99.3)
#'   counts every kNN link once, as squidpy does. With a kNN graph every cell
#'   sends k links at any density, so the weighting is not a density
#'   correction (density is handled by the label shuffles); in simulations it
#'   kept the calibration but lowered power.
#' @param n_perms Number of permutations for significance testing.
#' @param seed Random seed for reproducibility.
#' @param n_jobs Number of cores to use in parallel.
#'
#' @return Updated Seurat object; `misc[[paste0(cluster_key, "_nhood_enrichment")]]`
#'   holds `zscore`, `count`, `expected` (permutation mean), `log2_oe`,
#'   `log2_oe_raw`, `pvalue`, `padj` and `padj_bh` (see [nhood_enrichment()]).
#' @export
#' @examples
#' df <- generate_sim(close_ratio = 0.8, n_types = 4, n_cells = 300,
#'                    max_loc = 300, test_type = "distribute",
#'                    distance_param = 10, seed = 1)
#' df$sample_id <- "fov1"
#' seu <- sim_to_seurat(df)
#' seu <- nhood_enrichment.Seurat(seu, cluster_key = "cell_type",
#'                                neighbors.k = 10, n_perms = 20, n_jobs = 1)
#' res <- SeuratObject::Misc(seu, slot = "cell_type_nhood_enrichment")
#' round(res$zscore, 1)
nhood_enrichment.Seurat <- function(seurat_obj, cluster_key, neighbors.k = 30, connectivity_key = "nn", transformation = FALSE, n_perms = 100, seed = 1938493, n_jobs = 4) {
  if (!cluster_key %in% colnames(seurat_obj@meta.data)) {
    stop("Cluster key ", cluster_key, " not found in meta.data")
  }
  cluster_data <- seurat_obj@meta.data[[cluster_key]]

  all_nn <- list()
  cell_id <- vector()
  for (name in names(seurat_obj@images)) {
    coords <- seurat_obj[[name]]$centroids@coords %>%
      as.data.frame() %>%
      dplyr::mutate(cell = Cells(seurat_obj[[name]]))
    cells <- coords$cell
    rownames(coords) <- cells
    coords <- as.matrix(coords[, c("x", "y")])
    all_nn[[name]] <- .knn_graph(coords, neighbors.k, connectivity_key)
    cell_id <- c(cell_id, cells)
  }

  adj <- bdiag(all_nn)
  rownames(adj) <- colnames(adj) <- cell_id

  # --- normalize ---
  if (transformation) {
    degrees <- Matrix::colSums(adj) + 1
    adj <- adj / degrees
  }

  res <- .nhood_permutation_core(adj, cluster_data, transformation, n_perms, seed, n_jobs)

  seurat_obj@misc[[paste0(cluster_key, "_nhood_enrichment")]] <- res

  return(seurat_obj)
}

# Internal: iterative graph diffusion of a per-cell score, shared by
# cooccur_local() and its Seurat method. Each step computes
#   s <- (A + I) D^-1 s,   D = diag(colSums(A) + 1),
# starting from the previous step's result (so total mass is conserved).
# Diffusion stops early once the kurtosis of the normalized scores drops by
# less than 3 between consecutive steps (after more than 3 steps), or when
# the score is identically zero.
.diffuse_scores <- function(adj, local_score, maxnsteps, verbose = FALSE) {
  s <- matrix(local_score)
  if (maxnsteps < 1) return(s)
  degrees <- Matrix::colSums(adj) + 1
  prevmedkurt <- Inf
  for (i in seq_len(maxnsteps)) {
    s_norm <- s / degrees
    s <- as.matrix(adj %*% s_norm + s_norm)
    medkurt <- moments::kurtosis(prop.table(s, 2))
    if (is.nan(medkurt)) break
    if (prevmedkurt - medkurt < 3 && i > 3) {
      if (verbose) message("stopping after ", i, " steps")
      break
    }
    prevmedkurt <- medkurt
  }
  s
}

#' Local Co-occurrence Score (Seurat Method)
#'
#' @param seurat_obj A Seurat object.
#' @param cluster_x First cluster of interest.
#' @param cluster_y Second cluster of interest.
#' @param connectivity_key Graph type to use.
#' @param cluster_key Metadata column with cluster info.
#' @param sample_key Metadata column with sample ID.
#' @param neighbors.k Number of neighbors of the diffusion graph. It no longer
#'   limits the radius search; see `radius_k`.
#' @param radius Radius for proximity-based interaction.
#' @param maxnsteps Maximum number of diffusion steps (each step starts from
#'   the previous one; early stop on the kurtosis criterion, see
#'   [cooccur_local()]). `0` returns the raw indicator.
#' @param radius_k Optional cap on the cells returned per cell by the radius
#'   search, counting the cell itself. `NULL` (default) uses every cell within
#'   `radius`.
#'
#' @return A data.frame with local co-occurrence scores.
#' @export
#' @examples
#' df <- generate_sim(close_ratio = 0.8, n_types = 4, n_cells = 300,
#'                    max_loc = 300, test_type = "distribute",
#'                    distance_param = 10, seed = 1)
#' df$sample_id <- "fov1"
#' seu <- sim_to_seurat(df)
#' sc <- cooccur_local.Seurat(seu, cluster_x = "cell_type_1",
#'                            cluster_y = "cell_type_2",
#'                            cluster_key = "cell_type", sample_key = "sample_id",
#'                            neighbors.k = 10, radius = 20, maxnsteps = 1)
#' summary(sc[[1]])
cooccur_local.Seurat <- function(seurat_obj, cluster_x, cluster_y, connectivity_key = "nn", cluster_key = "seurat_clusters", sample_key = "sample_id", neighbors.k = 20, radius = 30, maxnsteps = 15, radius_k = NULL) {
  all_nn <- list()
  cell_id <- vector()

  for (name in names(seurat_obj@images)) {
    coords <- seurat_obj[[name]]$centroids@coords %>%
      as.data.frame() %>%
      dplyr::mutate(cell = Cells(seurat_obj[[name]]))
    cells <- coords$cell
    rownames(coords) <- cells
    coords <- as.matrix(coords[, c("x", "y")])

    all_nn[[name]] <- .knn_graph(coords, neighbors.k, connectivity_key)
    cell_id <- c(cell_id, cells)
  }


  adj <- bdiag(all_nn)
  rownames(adj) <- colnames(adj) <- cell_id

  local_score <- vector()

  for (name in names(seurat_obj@images)) {
    if(sample_key=="fov"){
      coords <- seurat_obj[[name]]$centroids@coords %>%
        as.data.frame() %>%
        dplyr::mutate(cell = Cells(seurat_obj[[name]])) %>%
        tibble::column_to_rownames(var = "cell") %>%
        dplyr::mutate(
          cluster = seurat_obj@meta.data[, cluster_key]
        )

    } else {
      coords <- seurat_obj[[name]]$centroids@coords %>%
        as.data.frame() %>%
        dplyr::mutate(cell = Cells(seurat_obj[[name]])) %>%
        tibble::column_to_rownames(var = "cell") %>%
        dplyr::mutate(
          cluster = seurat_obj@meta.data[
            seurat_obj@meta.data[, sample_key] == name,
            cluster_key
          ]
        )

    }

    local_score_ <- .radius_presence(as.matrix(coords[, c("x", "y")]), coords$cluster,
                                     cluster_x, cluster_y, radius, radius_k)
    names(local_score_) <- rownames(coords)
    local_score <- c(local_score, local_score_)
  }



  s <- .diffuse_scores(adj, local_score, maxnsteps, verbose = TRUE)

  df <- data.frame(cooccur_local_scores = as.numeric(s)) %>%
    magrittr::set_colnames(paste0("cooccur_local_", cluster_x, "_", cluster_y))
  rownames(df) <- names(local_score)
  return(df)
}

#' Local Co-occurrence Score (Generic method)
#'
#' @param df Data.frame of coordinates and cluster.
#' @param cluster_x First cluster of interest.
#' @param cluster_y Second cluster of interest.
#' @param connectivity_key Graph type to use.
#' @param neighbors.k Number of neighbors of the diffusion graph. It no longer
#'   limits the radius search; see `radius_k`.
#' @param radius Radius for neighborhood.
#' @param maxnsteps Maximum number of diffusion steps. Each step computes
#'   `s <- (A + I) D^-1 s` from the previous step; diffusion stops early once
#'   the kurtosis of the scores decreases by less than 3 between steps
#'   (checked after step 3). `0` returns the raw 0/1 indicator. Versions
#'   <= 0.99.1 always performed a single step regardless of `maxnsteps`.
#' @param radius_k Optional cap on the cells returned per cell by the radius
#'   search, counting the cell itself. `NULL` (default) uses every cell within
#'   `radius`.
#'
#' @return Data.frame with scores.
#' @export
#' @examples
#' df <- generate_sim(close_ratio = 0.8, n_types = 4, n_cells = 300,
#'                    max_loc = 300, test_type = "distribute",
#'                    distance_param = 10, seed = 1)
#' sc <- cooccur_local(df, cluster_x = "cell_type_1", cluster_y = "cell_type_2",
#'                     neighbors.k = 10, radius = 20)
#' summary(sc[[1]])
cooccur_local <- function(df, cluster_x, cluster_y, connectivity_key = "nn", neighbors.k = 20, radius = 30, maxnsteps = 1, radius_k = NULL) {
  #if (inherits(df, "Seurat")){
  #  return(cooccur_local.Seurat(df, cluster_x, cluster_y, connectivity_key, neighbors.k, radius, maxnsteps))
  #}
  coords <- df %>%
    dplyr::mutate(cell = rownames(.))
  cell_id <- coords$cell
  rownames(coords) <- coords$cell
  coords <- as.matrix(coords[, c("x", "y")])

  # transformation: advised for analytic p-value calculation.

  adj <- .knn_graph(coords, neighbors.k, connectivity_key)
  rownames(adj) <- colnames(adj) <- cell_id

  coords <- df %>%
    dplyr::mutate(cell = rownames(.)) %>%
    dplyr::rename( cluster = cell_type)

  local_score <- .radius_presence(as.matrix(coords[, c("x", "y")]), coords$cluster,
                                  cluster_x, cluster_y, radius, radius_k)
  names(local_score) <- rownames(coords)

  s <- .diffuse_scores(adj, local_score, maxnsteps)
  df_ <- data.frame(cooccur_local_scores = as.numeric(s)) %>%
    magrittr::set_colnames(paste0("cooccur_local_", cluster_x, "_", cluster_y))
  rownames(df_) <- names(local_score)
  return(df_)
}

#' Neighborhood Enrichment (Generic method)
#'
#' @param df Data.frame with spatial and cluster metadata.
#' @param cluster_key Column with cluster labels.
#' @param neighbors.k Number of neighbors to use.
#' @param connectivity_key Type of graph: "nn" or "snn".
#' @param transformation If `TRUE`, each link from cell u is weighted
#'   1 / (1 + d_u), with d_u the number of cells that chose u as a neighbour
#'   (a mild down-weighting of hub cells). `FALSE` (default since 0.99.3)
#'   counts every kNN link once, as squidpy does. With a kNN graph every cell
#'   sends k links at any density, so the weighting is not a density
#'   correction (density is handled by the label shuffles); in simulations it
#'   kept the calibration but lowered power.
#' @param n_perms Number of permutations.
#' @param seed Random seed.
#' @param n_jobs Number of parallel jobs. `1` runs sequentially.
#'
#' @return A list of cell-type x cell-type matrices:
#'   * `log2_oe`: effect size, log2 observed / expected, centred on the label
#'     shuffles so that it is 0 on average without interaction for any number
#'     of cells (the log of a ratio of small counts is otherwise biased
#'     slightly below 0 for rare cell types). Unlike the z-score, which grows
#'     with the number of cells, it is comparable across samples.
#'   * `log2_oe_raw`: log2((count + c) / (expected + c)) without centring,
#'     with c one mean edge weight.
#'   * `pvalue`: within-sample test per unordered pair (contacts i -> j and
#'     j -> i summed), two-sided normal p-value from the shuffles; use it for
#'     a single pre-specified pair.
#'   * `padj`: family-wise adjusted p-value over all K (K + 1) / 2 pairs by
#'     the Westfall-Young max-T permutation method (the share of shuffles
#'     whose largest |z| reaches the observed one). Calibrated for any
#'     number of cell types, including rare ones; its smallest value is
#'     1 / (n_perms + 1).
#'   * `padj_bh`: Benjamini-Hochberg on `pvalue` (can exceed its level when
#'     there are many rare cell types).
#'   * `zscore`, `count` (observed), `expected` (mean of the shuffles).
#'   * Directional statistics (row = centre cell type i, column = neighbour
#'     cell type j; not symmetric). The pair-level values above are nearly
#'     symmetric by construction (every i-j contact is also a j-i contact),
#'     so they cannot tell apart "type i is surrounded by j" and "type j is
#'     surrounded by i". Two directional questions are answered separately,
#'     on the unweighted kNN graph without the cell itself:
#'     - `contact[i, j]`: share of type-i cells with at least one type-j
#'       neighbour ("how much of population i touches j");
#'     - `dominance[i, j]`: share of type-i cells whose neighbours are at
#'       least half type j ("is the neighbourhood of i dominated by j").
#'     For each, `*_expected` is the mean over the same label shuffles,
#'     `*_log2_oe` the centred log2 ratio, `*_pvalue` a test per ordered pair
#'     and `*_padj` the max-T adjustment over all K x K ordered pairs.
#'     `contact` saturates near 1 when j is abundant, `dominance` is near 0
#'     when j is rare; read the two together.
#' @export
#' @examples
#' df <- generate_sim(close_ratio = 0.8, n_types = 4, n_cells = 300,
#'                    max_loc = 300, test_type = "distribute",
#'                    distance_param = 10, seed = 1)
#' res <- nhood_enrichment(df, cluster_key = "cell_type", neighbors.k = 10,
#'                         n_perms = 50, n_jobs = 1)
#' round(res$zscore, 1)
#' round(res$log2_oe, 2)
nhood_enrichment <- function(df, cluster_key, neighbors.k = 30, connectivity_key = "nn", transformation = FALSE, n_perms = 100, seed = 1938493, n_jobs = 4) {
  #if(inherits(df, "Seurat")){
  #  return(nhood_enrichment.Seurat(df, cluster_key, neighbors.k, connectivity_key, transformation, n_perms, seed, n_jobs))
  #}
  cluster_data <- df[,cluster_key]

  coords <- df %>%
    dplyr::mutate(cell = rownames(.))
  cell_id <- coords$cell
  rownames(coords) <- coords$cell
  coords <- as.matrix(coords[, c("x", "y")])

  adj <- .knn_graph(coords, neighbors.k, connectivity_key)
  rownames(adj) <- colnames(adj) <- cell_id

  # --- normalize ---
  if (transformation) {
    degrees <- Matrix::colSums(adj) + 1
    adj <- adj / degrees
  }

  .nhood_permutation_core(adj, cluster_data, transformation, n_perms, seed, n_jobs)
}

#' Default Manual Colors for Clusters
#'
#' @examples
#' head(manual_colors)
#' @export
manual_colors <- c(
  "0" = "#E41A1C", "1" = "#377EB8", "2" = "#4DAF4A", "3" = "#984EA3", "4" = "#FF7F00",
  "5" = "#FFFF33", "6" = "#A65628", "7" = "lightgrey", "8" = "#999999", "9" = "#66C2A5",
  "10" = "#67000D", "11" = "#8DA0CB", "12" = "#FFD92F", "13" = "#A6D854", "14" = "#E78AC3",
  "15" = "#FC8D62", "16" = "darkgrey", "17" = "#FEB24C", "18" = "#377EB8", "19" = "lightblue",
  "20" = "#FDE0EF", "21" = "#B8E186", "22" = "#66C2A5", "23" = "#A6D855", "24" = "#E78AC4",
  "25" = "#FC8D63", "26" = "brown", "27" = "#FEB25C", "28" = "#377EB9", "29" = "lightgreen",
  "30" = "#FDE1EF"
)
