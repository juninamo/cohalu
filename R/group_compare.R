# ---- Disease / condition group comparison utilities ----
#
# These functions add a sample-aware layer on top of the single-object
# analyses in cohalu_functions.R, so that scores can be computed
# per image (or aggregated per patient) and then compared between disease
# groups using standard or mixed-effects statistics, or a group-label
# permutation test.

#' @importFrom stats wilcox.test t.test p.adjust sd quantile pnorm pt aggregate as.formula lm rnorm setNames
#' @importFrom utils head
NULL

utils::globalVariables(c(
  "cluster_i", "cluster_j", "value", "group", "patient", "sample_id",
  "cell_type"
))

# ---- Phase 0: input dispatch ---------------------------------------------

#' Build a sample design table for disease-group comparisons
#'
#' Construct a data.frame mapping image / sample identifiers to a disease
#' group label (and optionally a patient ID). This is the shared metadata
#' table consumed by the `*_per_sample` helpers and by `compare_groups()`.
#'
#' @param obj A Seurat object, a list of Seurat objects, or a data.frame.
#'   For Seurat input, `sample_key` is expected to be a `meta.data` column
#'   whose values match the image names in `obj@images`. For data.frame
#'   input, `sample_key` is a column of the data.frame.
#' @param sample_key Name of the column / image identifier.
#' @param group_key Name of the column carrying the disease group (or any
#'   condition label).
#' @param patient_key Optional column name for patient identifier, used as a
#'   random effect (or permutation block) downstream. NULL to omit.
#'
#' @return A data.frame with columns `sample_id`, `group`, `patient`, and
#'   `source_index` (which list element of the input the sample came from).
#' @export
#' @examples
#' df <- generate_sim_groups(n_samples_per_group = 3,
#'                           group_close_ratio = list(case = 0.8, control = 0.2),
#'                           n_types = 4, n_cells = 200, max_loc = 250,
#'                           test_type = "distribute", distance_param = 8,
#'                           seed = 1)
#' build_sample_design(df, sample_key = "sample_id", group_key = "group",
#'                     patient_key = "patient")
build_sample_design <- function(obj, sample_key, group_key, patient_key = NULL) {
  if (inherits(obj, "Seurat")) {
    obj_list <- list(obj)
  } else if (is.list(obj) && length(obj) > 0L &&
             all(vapply(obj, inherits, logical(1), what = "Seurat"))) {
    obj_list <- obj
  } else if (is.data.frame(obj)) {
    df <- obj
    if (!sample_key %in% colnames(df)) {
      stop(sprintf("sample_key '%s' not found in data.frame columns", sample_key))
    }
    if (!group_key %in% colnames(df)) {
      stop(sprintf("group_key '%s' not found in data.frame columns", group_key))
    }
    sids <- unique(as.character(df[[sample_key]]))
    design <- data.frame(
      sample_id = sids,
      group = vapply(sids, function(s) {
        as.character(unique(df[[group_key]][df[[sample_key]] == s]))[1]
      }, character(1)),
      patient = if (is.null(patient_key)) NA_character_ else {
        if (!patient_key %in% colnames(df)) {
          stop(sprintf("patient_key '%s' not found in data.frame columns", patient_key))
        }
        vapply(sids, function(s) {
          as.character(unique(df[[patient_key]][df[[sample_key]] == s]))[1]
        }, character(1))
      },
      source_index = 1L,
      stringsAsFactors = FALSE
    )
    return(design)
  } else {
    stop("`obj` must be a Seurat object, a list of Seurat objects, or a data.frame.")
  }

  rows <- list()
  for (k in seq_along(obj_list)) {
    seu <- obj_list[[k]]
    md <- seu@meta.data
    if (!group_key %in% colnames(md)) {
      stop(sprintf("group_key '%s' not found in meta.data of object %d", group_key, k))
    }
    if (!is.null(patient_key) && !patient_key %in% colnames(md)) {
      stop(sprintf("patient_key '%s' not found in meta.data of object %d", patient_key, k))
    }
    for (img in names(seu@images)) {
      sel <- .image_cell_selector(seu, img, md, sample_key)
      if (!any(sel)) next
      gv <- unique(as.character(md[[group_key]][sel]))
      if (length(gv) > 1L) {
        warning(sprintf(
          "Image '%s' (object %d) maps to multiple groups (%s); using first.",
          img, k, paste(gv, collapse = ", ")
        ))
      }
      pv <- if (is.null(patient_key)) NA_character_ else {
        v <- unique(as.character(md[[patient_key]][sel]))
        if (length(v) > 1L) {
          warning(sprintf(
            "Image '%s' (object %d) maps to multiple patients (%s); using first.",
            img, k, paste(v, collapse = ", ")
          ))
        }
        v[1]
      }
      rows[[length(rows) + 1L]] <- data.frame(
        sample_id = img, group = gv[1], patient = pv, source_index = k,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0L) {
    stop("No samples were found. Check sample_key and the images slot of your object.")
  }
  do.call(rbind, rows)
}

# Internal: the package's single entry point to the RNG seed. Seeds are
# always explicit, user-supplied `seed` arguments.
.seed_rng <- function(seed) {
  set.seed(seed)
}

# Internal: set the RNG seed for the calling function only; the caller's
# .Random.seed is restored when that function exits.
.local_seed <- function(seed, envir = parent.frame()) {
  genv <- globalenv()
  old <- if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
    get(".Random.seed", envir = genv, inherits = FALSE)
  } else NULL
  do.call(on.exit, list(substitute(
    if (is.null(OLD)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
        rm(".Random.seed", envir = globalenv())
    } else assign(".Random.seed", OLD, envir = globalenv()),
    list(OLD = old)), add = TRUE), envir = envir)
  .seed_rng(seed)
  invisible()
}

# Internal: permutation null for paired / repeated-measures designs. Within
# each patient observed in both groups, the two group labels are swapped
# with probability 1/2 (sign flip of the patient's contribution); patients
# seen in only one group keep their labels. Exact when 2^n_patients <=
# n_perms.
.paired_perm_null <- function(v, g, pat, g1, g2, n_perms) {
  both <- names(which(tapply(g, pat, function(x) all(c(g1, g2) %in% x))))
  stat <- function(flip) {
    gg <- g
    sw <- pat %in% both[flip]
    gg[sw] <- ifelse(g[sw] == g1, g2, g1)
    mean(v[gg == g2]) - mean(v[gg == g1])
  }
  n <- length(both)
  if (2^n <= n_perms) {
    grid <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), n)))
    list(null = apply(grid, 1, stat), exact = TRUE)
  } else {
    list(null = replicate(n_perms, stat(sample(c(FALSE, TRUE), n, replace = TRUE))),
         exact = FALSE)
  }
}

# Internal: logical selector of the meta.data rows belonging to image `img`.
# Cell names are the ground truth (they are what links an image to its
# meta.data rows); the `sample_key` column is only used as a fallback when
# the image cells cannot be matched by name.
.image_cell_selector <- function(seu, img, md, sample_key) {
  sel <- rownames(md) %in% Seurat::Cells(seu[[img]])
  if (!any(sel) && sample_key %in% colnames(md)) {
    sel <- as.character(md[[sample_key]]) == img
  }
  sel
}

# Internal: warn once if one value per patient is expected but the input
# carries several rows (images) of the same patient for a single pair.
.warn_pseudoreplication <- function(df, pair_keys, patient_col, method) {
  if (is.null(patient_col) || !patient_col %in% colnames(df)) return(invisible())
  pat <- as.character(df[[patient_col]])
  if (all(is.na(pat))) return(invisible())
  key <- do.call(paste, c(df[pair_keys], list(pat), sep = "\r"))
  if (anyDuplicated(key[!is.na(pat)])) {
    warning(
      "Several rows per patient were found for the same cluster pair; ",
      "method '", method, "' treats them as independent observations ",
      "(pseudoreplication), which inflates the false-positive rate. Use ",
      "unit = \"patient\" in the *_per_sample() call, or method = \"lmm\" / ",
      "\"perm\" with patient_key.",
      call. = FALSE)
  }
  invisible()
}

# Internal: per-sample cell counts per cluster, used to annotate pair tables
# with the abundance of the two cell types (useful for filtering rare pairs
# and for diagnosing composition effects).
.add_composition <- function(tidy_df, samples, pair_cols = c("cluster_i", "cluster_j")) {
  comp <- do.call(rbind, lapply(names(samples), function(sid) {
    tab <- table(as.character(samples[[sid]]$coords$cluster))
    data.frame(sample_id = sid, cluster = names(tab), n = as.integer(tab),
               n_cells = nrow(samples[[sid]]$coords), stringsAsFactors = FALSE)
  }))
  lookup <- stats::setNames(comp$n, paste(comp$sample_id, comp$cluster, sep = "\r"))
  totals <- stats::setNames(comp$n_cells, comp$sample_id)
  totals <- totals[!duplicated(names(totals))]
  tidy_df$n_cells <- unname(totals[tidy_df$sample_id])
  for (k in seq_along(pair_cols)) {
    n <- lookup[paste(tidy_df$sample_id, tidy_df[[pair_cols[k]]], sep = "\r")]
    n[is.na(n)] <- 0L
    tidy_df[[paste0("n_", c("i", "j")[k])]] <- unname(n)
  }
  tidy_df
}

# Internal: normalize any supported input into a per-sample list with
# `coords` (data.frame with x, y, cluster, rownames = cell IDs) plus
# `group` and `patient` metadata. Returns list(samples = ..., design = ...).
.as_sample_list <- function(obj, sample_key, cluster_key,
                            group_key = NULL, patient_key = NULL) {
  design <- if (!is.null(group_key)) {
    build_sample_design(obj, sample_key = sample_key,
                        group_key = group_key, patient_key = patient_key)
  } else {
    NULL
  }
  out <- list()

  if (is.data.frame(obj)) {
    df <- obj
    if (!cluster_key %in% colnames(df)) {
      stop(sprintf("cluster_key '%s' not found in data.frame columns", cluster_key))
    }
    if (!all(c("x", "y") %in% colnames(df))) {
      stop("data.frame input must contain 'x' and 'y' columns")
    }
    sids <- unique(as.character(df[[sample_key]]))
    for (s in sids) {
      sub <- df[as.character(df[[sample_key]]) == s, , drop = FALSE]
      cells <- if (!is.null(rownames(sub)) && !all(rownames(sub) == as.character(seq_len(nrow(sub))))) {
        rownames(sub)
      } else {
        paste0(s, "_cell", seq_len(nrow(sub)))
      }
      coords <- data.frame(
        x = sub$x, y = sub$y,
        cluster = as.character(sub[[cluster_key]]),
        row.names = cells, stringsAsFactors = FALSE
      )
      d <- if (!is.null(design)) design[design$sample_id == s, , drop = FALSE] else NULL
      out[[s]] <- list(
        coords = coords,
        group = if (!is.null(d) && nrow(d)) d$group[1] else NA_character_,
        patient = if (!is.null(d) && nrow(d)) d$patient[1] else NA_character_
      )
    }
    return(list(samples = out, design = design))
  }

  obj_list <- if (inherits(obj, "Seurat")) list(obj) else obj
  if (!is.list(obj_list)) {
    stop("`obj` must be a Seurat object, list of Seurat objects, or data.frame.")
  }

  for (k in seq_along(obj_list)) {
    seu <- obj_list[[k]]
    md <- seu@meta.data
    if (!cluster_key %in% colnames(md)) {
      stop(sprintf("cluster_key '%s' not found in meta.data of object %d", cluster_key, k))
    }
    for (img in names(seu@images)) {
      coords <- as.data.frame(seu[[img]]$centroids@coords)
      cells <- Seurat::Cells(seu[[img]])
      rownames(coords) <- cells
      sel <- .image_cell_selector(seu, img, md, sample_key)
      md_sel <- md[sel, , drop = FALSE]
      common <- intersect(cells, rownames(md_sel))
      if (length(common) == 0L) {
        # fall back: assume order matches
        if (nrow(md_sel) != nrow(coords)) {
          stop(sprintf(
            "Image '%s' has %d cells but matching meta.data has %d rows.",
            img, nrow(coords), nrow(md_sel)
          ))
        }
        coords$cluster <- as.character(md_sel[[cluster_key]])
      } else {
        coords <- coords[common, , drop = FALSE]
        coords$cluster <- as.character(md_sel[common, cluster_key])
      }
      d <- if (!is.null(design)) {
        design[design$sample_id == img & design$source_index == k, , drop = FALSE]
      } else NULL
      out[[img]] <- list(
        coords = coords,
        group = if (!is.null(d) && nrow(d)) d$group[1] else NA_character_,
        patient = if (!is.null(d) && nrow(d)) d$patient[1] else NA_character_
      )
    }
  }
  list(samples = out, design = design)
}

# Internal: stack per-sample matrix results into a tidy long data.frame.
# `per_sample` is a named list of named lists; for each sample we take
# the matrices specified by `value_map`, vectorize them, and stack.
# value_map: named character vector list-element-name -> output-column-name.
.stack_pair_results <- function(per_sample, design, value_map) {
  rows <- list()
  for (sid in names(per_sample)) {
    res <- per_sample[[sid]]
    if (is.null(res)) next
    mat0 <- res[[names(value_map)[1]]]
    if (is.null(mat0)) next
    rn <- rownames(mat0); cn <- colnames(mat0)
    grid <- expand.grid(cluster_i = rn, cluster_j = cn, stringsAsFactors = FALSE)
    for (k in seq_along(value_map)) {
      src <- names(value_map)[k]
      dst <- unname(value_map[k])
      grid[[dst]] <- as.vector(res[[src]])
    }
    grid$sample_id <- sid
    rows[[sid]] <- grid
  }
  if (length(rows) == 0L) return(NULL)
  out <- do.call(rbind, rows)
  out$cluster_i <- sub("^Cluster", "", out$cluster_i)
  out$cluster_j <- sub("^Cluster", "", out$cluster_j)
  if (!is.null(design)) {
    out <- merge(out, design[, c("sample_id", "group", "patient")],
                 by = "sample_id", all.x = TRUE)
  }
  rownames(out) <- NULL
  out
}

# Internal: aggregate sample-level rows to patient-level rows by mean.
.aggregate_to_patient <- function(tidy_df, pair_keys = NULL) {
  if (!"patient" %in% colnames(tidy_df) || all(is.na(tidy_df$patient))) {
    stop("patient column is missing or all NA; cannot aggregate to patient level.")
  }
  # distance columns (e.g. `r` from colocalization_per_sample()) are keys, not values
  if (is.null(pair_keys)) pair_keys <- intersect(c("cluster_i", "cluster_j", "r"), colnames(tidy_df))
  num_cols <- vapply(tidy_df, is.numeric, logical(1))
  num_cols <- setdiff(names(num_cols)[num_cols], pair_keys)
  by_cols <- c(list(patient = tidy_df$patient),
               if (length(pair_keys)) as.list(tidy_df[, pair_keys, drop = FALSE])
               else list(cluster_i = rep("", nrow(tidy_df)), cluster_j = rep("", nrow(tidy_df))),
               list(group = tidy_df$group))
  agg <- aggregate(tidy_df[, num_cols, drop = FALSE], by = by_cols,
                   FUN = function(x) mean(x, na.rm = TRUE))
  agg$sample_id <- agg$patient
  agg
}

#' Aggregate per-image scores to one row per patient
#'
#' Average every numeric column of a `*_per_sample()` result across the
#' images of each patient (within cluster pair and group). Use this to
#' obtain patient-level values for `compare_groups(method = "wilcox" / "t")`
#' from an image-level result without recomputing the permutations
#' (equivalent to `unit = "patient"`).
#'
#' @param per_sample_df Output of a `*_per_sample()` helper with a
#'   non-missing `patient` column.
#' @param pair_keys Columns identifying a comparison (kept separate, not
#'   averaged). Defaults to `cluster_i`, `cluster_j` and, when present, the
#'   distance `r` of [colocalization_per_sample()].
#'
#' @return A data.frame with one row per `patient x cluster_i x cluster_j`
#'   (x `r`);
#'   `sample_id` is set to the patient ID.
#' @export
#' @examples
#' df <- generate_sim_groups(n_samples_per_group = 2, n_images_per_patient = 2,
#'                           group_close_ratio = list(case = 0.8, control = 0.2),
#'                           n_types = 3, n_cells = 150, max_loc = 250,
#'                           test_type = "distribute", distance_param = 8)
#' ps <- nhood_enrichment_per_sample(df, sample_key = "sample_id",
#'                                   group_key = "group",
#'                                   cluster_key = "cell_type",
#'                                   patient_key = "patient",
#'                                   neighbors.k = 8, n_perms = 20, n_jobs = 1)
#' pp <- summarize_by_patient(ps)
#' table(pp$sample_id)
summarize_by_patient <- function(per_sample_df, pair_keys = NULL) {
  out <- .aggregate_to_patient(as.data.frame(per_sample_df), pair_keys)
  attr(out, "spatial_design") <- attr(per_sample_df, "spatial_design")
  attr(out, "value_columns") <- attr(per_sample_df, "value_columns")
  class(out) <- unique(c("cohaluSample", class(out)))
  out
}

# ---- Phase 1a: nhood_enrichment_per_sample -------------------------------

#' Per-sample neighborhood enrichment
#'
#' Run [nhood_enrichment()] independently for each sample (image) in the
#' input, returning a tidy long data.frame with one row per
#' `sample_id x cluster_i x cluster_j`. The cluster factor levels are
#' harmonized across samples so the resulting table is suitable for
#' downstream group comparison with [compare_groups()].
#'
#' @inheritParams nhood_enrichment
#' @param obj A Seurat object, a list of Seurat objects, or a data.frame
#'   with x, y, cluster_key, and sample_key columns.
#' @param sample_key Column / image identifier defining a sample.
#' @param group_key Column carrying the disease / condition label.
#' @param patient_key Optional column with patient ID (random effect /
#'   permutation block).
#' @param unit "image" returns one row per image; "patient" averages across
#'   images within a patient (requires `patient_key`).
#' @param cluster_levels Optional character vector of cluster levels to use
#'   as common dimnames across samples. Defaults to the union across all
#'   samples.
#'
#' @return A data.frame (also tagged with class `cohaluSample`) with
#'   columns `sample_id`, `cluster_i`, `cluster_j`, `zscore`, `count`,
#'   `expected`, `log2_oe`, `group`, `patient`, plus `n_cells` (cells in the
#'   sample) and `n_i` / `n_j` (cells of `cluster_i` / `cluster_j`).
#'
#' @section Choosing the value to compare:
#' The permutation z-score measures *statistical evidence* within one sample
#' and grows roughly with the square root of the number of cells, so samples
#' with more cells (larger images, denser tissue) get larger |z| for the same
#' spatial pattern. For between-group comparison, `log2_oe`
#' (log2 observed / permutation-expected) is an effect size that does not
#' scale with sample size and is usually the better choice. `count` is
#' additionally confounded by cell-type composition and should not be
#' compared directly.
#' @export
#' @examples
#' df <- generate_sim_groups(n_samples_per_group = 3,
#'                           group_close_ratio = list(case = 0.8, control = 0.2),
#'                           n_types = 4, n_cells = 200, max_loc = 250,
#'                           test_type = "distribute", distance_param = 8,
#'                           seed = 1)
#' ps <- nhood_enrichment_per_sample(df, sample_key = "sample_id",
#'                                   group_key = "group",
#'                                   cluster_key = "cell_type",
#'                                   patient_key = "patient",
#'                                   neighbors.k = 8, n_perms = 20, n_jobs = 1)
#' head(ps)
nhood_enrichment_per_sample <- function(obj, sample_key, group_key, cluster_key,
                                        patient_key = NULL,
                                        unit = c("image", "patient"),
                                        cluster_levels = NULL,
                                        neighbors.k = 30,
                                        connectivity_key = "nn",
                                        transformation = FALSE,
                                        n_perms = 100,
                                        seed = 1938493,
                                        n_jobs = 1) {
  unit <- match.arg(unit)
  sl <- .as_sample_list(obj, sample_key = sample_key,
                        cluster_key = cluster_key,
                        group_key = group_key, patient_key = patient_key)
  samples <- sl$samples
  design <- sl$design

  if (is.null(cluster_levels)) {
    cluster_levels <- sort(unique(unlist(lapply(samples, function(s) unique(s$coords$cluster)))))
  }

  per_sample <- list()
  for (sid in names(samples)) {
    coords <- samples[[sid]]$coords
    if (nrow(coords) < 2L) {
      warning(sprintf("Sample '%s' has fewer than 2 cells; skipping.", sid))
      next
    }
    df_in <- data.frame(
      x = coords$x, y = coords$y,
      stringsAsFactors = FALSE,
      row.names = rownames(coords)
    )
    df_in[[cluster_key]] <- factor(coords$cluster, levels = cluster_levels)
    res <- tryCatch(
      nhood_enrichment(df_in, cluster_key = cluster_key,
                       neighbors.k = neighbors.k,
                       connectivity_key = connectivity_key,
                       transformation = transformation,
                       n_perms = n_perms,
                       seed = seed,
                       n_jobs = n_jobs),
      error = function(e) {
        warning(sprintf("nhood_enrichment failed for sample '%s': %s", sid, conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(res)) per_sample[[sid]] <- res
  }

  tidy_df <- .stack_pair_results(per_sample, design,
                                 value_map = c(zscore = "zscore", count = "count",
                                               expected = "expected",
                                               log2_oe = "log2_oe"))
  if (is.null(tidy_df) || nrow(tidy_df) == 0L) {
    stop("No per-sample results were produced.")
  }
  tidy_df <- .add_composition(tidy_df, samples)
  if (unit == "patient") tidy_df <- .aggregate_to_patient(tidy_df)
  attr(tidy_df, "spatial_design") <- design
  attr(tidy_df, "value_columns") <- c("zscore", "log2_oe", "count", "expected")
  class(tidy_df) <- c("cohaluSample", class(tidy_df))
  tidy_df
}

# ---- Phase 1b: cooccur_local_per_sample ---------------------------------

#' Per-sample local co-occurrence summary
#'
#' Run [cooccur_local()] for one (cluster_x, cluster_y) pair on each sample
#' and summarize the per-cell scores at the sample level. Output is a tidy
#' data.frame suitable for group comparison.
#'
#' @inheritParams cooccur_local
#' @param obj A Seurat object, a list of Seurat objects, or a data.frame.
#' @param sample_key,group_key,cluster_key,patient_key Sample / group /
#'   cluster / patient column names. See [build_sample_design()].
#' @param unit "image" or "patient" (averaged across images of the same
#'   patient).
#' @param summarize Character vector of summary statistics to compute:
#'   "mean", "q90" (90th percentile), "pos_rate" (fraction of cells with
#'   score > 0), and / or "log2_oe": log2 of the number of `cluster_x`-
#'   `cluster_y` pairs within `radius` of each cell, summed over cells, over
#'   its expectation under label permutation (see [cooccur_local_oe()]). `log2_oe` is adjusted for the
#'   abundance of the two cell types and for cell density and is the
#'   recommended summary for group comparison; "mean" is unchanged by the
#'   (mass-conserving) diffusion and grows with abundance.
#'
#' @return A data.frame with one row per sample carrying the requested
#'   summary statistics, plus `n_cells`, `n_i` and `n_j` (cells of
#'   `cluster_x` / `cluster_y`).
#'
#' @section Caution:
#' The local score has no permutation null, so its sample-level summaries
#' increase with the abundance of `cluster_x` and `cluster_y`. When the two
#' groups differ in cell-type composition, check `n_i` / `n_j` (or adjust
#' for them with `covariates` in [compare_groups()]) before interpreting a
#' group difference as a change in co-localization.
#' @export
#' @examples
#' df <- generate_sim_groups(n_samples_per_group = 3,
#'                           group_close_ratio = list(case = 0.8, control = 0.2),
#'                           n_types = 4, n_cells = 200, max_loc = 250,
#'                           test_type = "distribute", distance_param = 8,
#'                           seed = 1)
#' cooccur_local_per_sample(df, sample_key = "sample_id", group_key = "group",
#'                          cluster_key = "cell_type",
#'                          cluster_x = "cell_type_1", cluster_y = "cell_type_2",
#'                          patient_key = "patient", neighbors.k = 10,
#'                          radius = 20)
cooccur_local_per_sample <- function(obj, sample_key, group_key, cluster_key,
                                     cluster_x, cluster_y,
                                     patient_key = NULL,
                                     unit = c("image", "patient"),
                                     neighbors.k = 20, radius = 30, maxnsteps = 1,
                                     summarize = c("mean", "q90", "pos_rate", "log2_oe"),
                                     radius_k = NULL) {
  unit <- match.arg(unit)
  summarize <- match.arg(summarize, several.ok = TRUE)
  sl <- .as_sample_list(obj, sample_key = sample_key,
                        cluster_key = cluster_key,
                        group_key = group_key, patient_key = patient_key)
  samples <- sl$samples
  design <- sl$design

  rows <- list()
  for (sid in names(samples)) {
    coords <- samples[[sid]]$coords
    if (nrow(coords) < 2L) next
    df_in <- data.frame(
      x = coords$x, y = coords$y,
      cell_type = as.character(coords$cluster),
      row.names = rownames(coords),
      stringsAsFactors = FALSE
    )
    res <- tryCatch(
      cooccur_local(df_in, cluster_x = cluster_x, cluster_y = cluster_y,
                    neighbors.k = neighbors.k, radius = radius,
                    maxnsteps = maxnsteps, radius_k = radius_k),
      error = function(e) {
        warning(sprintf("cooccur_local failed for sample '%s': %s", sid, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(res)) next
    sc <- as.numeric(res[[1]])
    row <- data.frame(
      sample_id = sid,
      cluster_i = cluster_x,
      cluster_j = cluster_y,
      stringsAsFactors = FALSE
    )
    if ("mean" %in% summarize) row$mean <- mean(sc, na.rm = TRUE)
    if ("q90" %in% summarize) row$q90 <- as.numeric(quantile(sc, 0.9, na.rm = TRUE))
    if ("pos_rate" %in% summarize) row$pos_rate <- mean(sc > 0, na.rm = TRUE)
    if ("log2_oe" %in% summarize) {
      nb <- .radius_pairs(as.matrix(df_in[, c("x", "y")]), radius, radius_k)
      ie <- .pairs_expected(df_in$cell_type, nb, cluster_x, cluster_y)
      row$log2_oe <- log2(sum(ie$pairs) / max(sum(ie$expected), .Machine$double.eps))
    }
    row$n_cells <- nrow(coords)
    row$n_i <- sum(coords$cluster == cluster_x)
    row$n_j <- sum(coords$cluster == cluster_y)
    rows[[sid]] <- row
  }
  if (length(rows) == 0L) stop("No per-sample results were produced.")
  tidy_df <- do.call(rbind, rows)
  if (!is.null(design)) {
    tidy_df <- merge(tidy_df, design[, c("sample_id", "group", "patient")],
                     by = "sample_id", all.x = TRUE)
  }
  if (unit == "patient") tidy_df <- .aggregate_to_patient(tidy_df)
  attr(tidy_df, "spatial_design") <- design
  attr(tidy_df, "value_columns") <- intersect(summarize, colnames(tidy_df))
  class(tidy_df) <- c("cohaluSample", class(tidy_df))
  tidy_df
}

# ---- Phase 1b: cooccur_ratio_per_sample ---------------------------------

# Internal: radius-based co-occurrence count + ratio for a single sample.
.radius_count_one <- function(coords_xy, clusters, all_clusters, radius, k = NULL) {
  co_occur_count <- .radius_count_matrix(coords_xy, clusters, all_clusters, radius, k)
  ratio_mat <- compute_co_occurrence_ratio(co_occur_count)
  list(co_occur_count = co_occur_count, ratio_mat = ratio_mat)
}

#' Per-sample radius-based co-occurrence ratio
#'
#' Compute the radius-based co-occurrence count and enrichment ratio per
#' sample, with cluster dimnames harmonized across samples.
#'
#' @inheritParams calc_co_occurrence_for_radius
#' @param obj A Seurat object, list of Seurat objects, or data.frame.
#' @param sample_key,group_key,cluster_key,patient_key See
#'   [build_sample_design()].
#' @param unit "image" or "patient".
#' @param cluster_levels Optional vector of cluster levels.
#'
#' @return A data.frame with one row per `sample_id x cluster_i x cluster_j`,
#'   columns `ratio`, `count`, `group`, `patient`, `n_cells`, `n_i`, `n_j`.
#' @export
#' @examples
#' df <- generate_sim_groups(n_samples_per_group = 3,
#'                           group_close_ratio = list(case = 0.8, control = 0.2),
#'                           n_types = 4, n_cells = 200, max_loc = 250,
#'                           test_type = "distribute", distance_param = 8,
#'                           seed = 1)
#' rs <- cooccur_ratio_per_sample(df, sample_key = "sample_id",
#'                                group_key = "group", cluster_key = "cell_type",
#'                                patient_key = "patient", radius = 20)
#' head(rs)
cooccur_ratio_per_sample <- function(obj, sample_key, group_key, cluster_key,
                                     patient_key = NULL,
                                     unit = c("image", "patient"),
                                     radius = 30, k = NULL,
                                     cluster_levels = NULL) {
  unit <- match.arg(unit)
  sl <- .as_sample_list(obj, sample_key = sample_key,
                        cluster_key = cluster_key,
                        group_key = group_key, patient_key = patient_key)
  samples <- sl$samples
  design <- sl$design

  if (is.null(cluster_levels)) {
    cluster_levels <- sort(unique(unlist(lapply(samples, function(s) unique(s$coords$cluster)))))
  }

  per_sample <- list()
  for (sid in names(samples)) {
    coords <- samples[[sid]]$coords
    if (nrow(coords) < 2L) next
    res <- tryCatch(
      .radius_count_one(coords_xy = as.matrix(coords[, c("x", "y")]),
                        clusters = as.character(coords$cluster),
                        all_clusters = cluster_levels,
                        radius = radius, k = k),
      error = function(e) {
        warning(sprintf("ratio computation failed for sample '%s': %s",
                        sid, conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(res)) per_sample[[sid]] <- res
  }

  tidy_df <- .stack_pair_results(per_sample, design,
                                 value_map = c(ratio_mat = "ratio",
                                               co_occur_count = "count"))
  if (is.null(tidy_df) || nrow(tidy_df) == 0L) {
    stop("No per-sample results were produced.")
  }
  tidy_df <- .add_composition(tidy_df, samples)
  if (unit == "patient") tidy_df <- .aggregate_to_patient(tidy_df)
  attr(tidy_df, "spatial_design") <- design
  attr(tidy_df, "value_columns") <- c("ratio", "count")
  class(tidy_df) <- c("cohaluSample", class(tidy_df))
  tidy_df
}

# ---- Phase 1b: interaction_spot_per_sample ------------------------------

#' Per-sample interaction-spot summary
#'
#' Apply [search_interaction_spot()] to each image of a Seurat object and
#' summarize how many connected interaction spots were detected and their
#' mean size, ready for group comparison.
#'
#' @param seurat_object A Seurat object (lists are not supported here).
#' @param sample_key,group_key,patient_key Column names; see
#'   [build_sample_design()].
#' @param cluster_col Cluster column in meta.data.
#' @param target_cluster Target cluster(s) of interest.
#' @param cell_id Vector of cell IDs to include. Defaults to all cells.
#' @param radius Radius defining neighborhood.
#' @param n_min Minimum number of cells per spot.
#' @param neighbors.k Optional cap on the cells returned per cell, counting
#'   the cell itself. `NULL` (default) uses every cell within `radius`.
#'
#' @return A data.frame with one row per sample carrying `n_spots`,
#'   `mean_spot_size` (number of cells), `n_cells` and `spots_per_1k_cells`.
#'   Raw `n_spots` scales with image size; compare `spots_per_1k_cells`
#'   between groups when images differ in size. Samples for which the spot
#'   search failed get `NA` (with a warning), not 0.
#' @export
#' @examples
#' df <- generate_sim_groups(n_samples_per_group = 3,
#'                           group_close_ratio = list(case = 0.8, control = 0.2),
#'                           n_types = 4, n_cells = 200, max_loc = 250,
#'                           test_type = "distribute", distance_param = 8,
#'                           seed = 1)
#' seu <- sim_to_seurat(df)
#' interaction_spot_per_sample(seu, sample_key = "sample_id",
#'                             group_key = "group", cluster_col = "cell_type",
#'                             target_cluster = c("cell_type_1", "cell_type_2"),
#'                             radius = 15, n_min = 3)
interaction_spot_per_sample <- function(seurat_object, sample_key, group_key,
                                        cluster_col, target_cluster,
                                        cell_id = NULL,
                                        patient_key = NULL,
                                        radius, n_min, neighbors.k = NULL) {
  if (!inherits(seurat_object, "Seurat")) {
    stop("interaction_spot_per_sample currently supports Seurat object input only.")
  }
  design <- build_sample_design(seurat_object, sample_key = sample_key,
                                group_key = group_key, patient_key = patient_key)
  if (!"cell" %in% colnames(seurat_object@meta.data)) {
    # search_interaction_spot() joins on a `cell` column.
    seurat_object@meta.data$cell <- rownames(seurat_object@meta.data)
  }
  if (is.null(cell_id)) {
    if ("cell" %in% colnames(seurat_object@meta.data)) {
      cell_id <- as.character(seurat_object@meta.data$cell)
    } else {
      cell_id <- rownames(seurat_object@meta.data)
    }
  }
  rows <- list()
  for (fov in names(seurat_object@images)) {
    spots <- tryCatch(
      search_interaction_spot(seurat_object, fov = fov, radius = radius, n_min = n_min,
                              neighbors.k = neighbors.k, cell_id = cell_id,
                              cluster_col = cluster_col,
                              target_cluster = target_cluster),
      error = function(e) {
        warning(sprintf("search_interaction_spot failed for '%s': %s",
                        fov, conditionMessage(e)))
        NULL
      }
    )
    # A failed search is missing data (NA), not "zero spots".
    n_spots <- if (!is.null(spots)) length(unique(spots$cluster_id)) else NA_integer_
    mean_size <- if (!is.null(spots) && nrow(spots) > 0L) {
      sizes <- unique(spots[, c("cluster_id", "n_all_cells")])
      mean(sizes$n_all_cells, na.rm = TRUE)
    } else NA_real_
    n_cells_fov <- sum(Seurat::Cells(seurat_object[[fov]]) %in% cell_id)
    d <- design[design$sample_id == fov, , drop = FALSE]
    rows[[fov]] <- data.frame(
      sample_id = fov,
      group = if (nrow(d)) d$group[1] else NA_character_,
      patient = if (nrow(d)) d$patient[1] else NA_character_,
      target_cluster = paste(target_cluster, collapse = ","),
      n_spots = n_spots,
      mean_spot_size = mean_size,
      n_cells = n_cells_fov,
      spots_per_1k_cells = if (n_cells_fov > 0) 1000 * n_spots / n_cells_fov else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  attr(out, "spatial_design") <- design
  attr(out, "value_columns") <- c("spots_per_1k_cells", "n_spots", "mean_spot_size")
  class(out) <- c("cohaluSample", class(out))
  out
}

# ---- Phase 2: compare_groups --------------------------------------------

#' Compare disease groups across samples
#'
#' Given a tidy per-sample data.frame produced by one of the
#' `*_per_sample()` helpers, run a per-cluster-pair statistical test
#' between groups and return tidy results with multiple-testing adjustment.
#'
#' @section Unit of analysis:
#' In a case-control study the independent unit is the *patient*, not the
#' image. With several images per patient either aggregate first
#' (`unit = "patient"` in the `*_per_sample()` call) and use
#' `"wilcox"` / `"t"`, or keep images and use `"lmm"` or `"perm"` with
#' `patient_key`. `compare_groups()` warns when `"wilcox"` / `"t"` (or
#' `"perm"` without `patient_key`) are given several rows per patient.
#'
#' @param per_sample_df Tidy data.frame, typically the output of
#'   [nhood_enrichment_per_sample()], [cooccur_ratio_per_sample()],
#'   [cooccur_local_per_sample()], or [interaction_spot_per_sample()].
#' @param value Name of the column to test. Defaults to `"log2_oe"`, the
#'   effect size that does not grow with the number of cells (other examples:
#'   "zscore", "ratio", "mean", "spots_per_1k_cells").
#' @param group_key Name of the group column. Defaults to "group".
#' @param patient_key Patient column. If `NULL` and the data contain a
#'   `patient` column with several images per patient, that column is used
#'   automatically (a message says so). With `"wilcox"` / `"t"`,
#'   images are averaged within patient (see `unit`); required for
#'   `method = "lmm"` to use as a random effect; used as the permutation
#'   block for `method = "perm"`. If `NULL` and the data has a `patient` column, that
#'   column is only used to detect pseudoreplication.
#' @param method Statistical test:
#'   * "wilcox" — Wilcoxon rank-sum (two groups); exact p-values for small
#'     samples without ties, normal approximation otherwise.
#'   * "t" — Welch's t-test (two groups)
#'   * "lmm" — linear mixed model `value ~ group + covariates + (1|patient)`
#'     via `lme4::lmer`. p-values use Satterthwaite degrees of freedom when
#'     the `lmerTest` package is installed, and otherwise a t reference with
#'     `n_patients - n_fixed_effects` degrees of freedom (group is a
#'     patient-level factor). Falls back to `lm()` when every patient
#'     contributes a single row. Requires the `lme4` package.
#'   * "perm" — group-label permutation test on the mean difference
#'     (blocked by patient if `patient_key` is supplied). If patients
#'     appear in both groups (paired / repeated-measures designs such as
#'     pre- vs post-treatment), labels are instead swapped *within*
#'     patients. All relabelings are enumerated when there are at most
#'     `n_perms` of them (exact test).
#'   * "signrank" — paired design: values are averaged per patient and
#'     group, and the per-patient differences (test minus reference) are
#'     tested with the Wilcoxon signed-rank test. Requires `patient_key`;
#'     patients observed in only one group are dropped.
#' @param n_perms Number of permutations for `method = "perm"`.
#' @param adjust Multiple-testing adjustment method passed to
#'   [stats::p.adjust()].
#' @param pair_keys Column names that together identify a cluster pair.
#'   Defaults to `c("cluster_i", "cluster_j")`. Set to a single column name
#'   for cases like interaction_spot_per_sample (e.g. "target_cluster").
#' @param ref_group Name of the reference group (e.g. "control"). The effect
#'   is reported as `mean_test - mean_ref`. When `NULL`, groups are sorted
#'   alphabetically and the first is used as reference (a message says
#'   which); set it explicitly, because e.g. "case" sorts before "control".
#' @param covariates Optional character vector of additional columns (e.g.
#'   age, sex, batch, `n_cells`) included as fixed effects. Only used with
#'   `method = "lmm"`.
#' @param symmetric If `TRUE`, test each unordered pair once: the (i, j)
#'   and (j, i) values of each sample are averaged and reported under
#'   `cluster_i <= cluster_j`. Neighborhood enrichment scores are
#'   (nearly) symmetric, so testing both (i, j) and (j, i) doubles the
#'   multiple-testing burden without adding information.
#' @param min_n_per_group Minimum number of finite observations required in
#'   each group; pairs below this are skipped.
#' @param seed Random seed for the permutation test.
#' @param unit `"patient"` (default): with `patient_key` and `method =
#'   "wilcox"` or `"t"`, images are first averaged within patient, so the
#'   patient is the unit of analysis. `"image"` tests image-level rows as
#'   given (pseudoreplication when patients have several images; a warning
#'   is issued). `"lmm"`, `"perm"` and `"signrank"` always account for
#'   patients through `patient_key`.
#'
#' @return A data.frame with the cluster pair columns, group sizes and
#'   means, `effect` (test group mean minus reference group mean), for
#'   `method = "lmm"` also `estimate` (the adjusted model coefficient),
#'   `statistic`, raw `p`, and `padj`. Attributes `method`, `groups`,
#'   `value` and `p_method` describe the test.
#' @export
#' @examples
#' df <- generate_sim_groups(n_samples_per_group = 3,
#'                           group_close_ratio = list(case = 0.8, control = 0.2),
#'                           n_types = 4, n_cells = 200, max_loc = 250,
#'                           test_type = "distribute", distance_param = 8,
#'                           seed = 1)
#' ps <- nhood_enrichment_per_sample(df, sample_key = "sample_id",
#'                                   group_key = "group",
#'                                   cluster_key = "cell_type",
#'                                   patient_key = "patient",
#'                                   neighbors.k = 8, n_perms = 20, n_jobs = 1)
#' cmp <- compare_groups(ps, value = "log2_oe", method = "wilcox",
#'                       ref_group = "control", symmetric = TRUE)
#' head(cmp)
compare_groups <- function(per_sample_df,
                           value = "log2_oe",
                           group_key = "group",
                           patient_key = NULL,
                           method = c("wilcox", "t", "lmm", "perm", "signrank"),
                           n_perms = 1000,
                           adjust = "BH",
                           pair_keys = c("cluster_i", "cluster_j"),
                           ref_group = NULL,
                           covariates = NULL,
                           symmetric = FALSE,
                           min_n_per_group = 2,
                           seed = 1234,
                           unit = c("patient", "image")) {
  method <- match.arg(method)
  unit <- match.arg(unit)
  if (!value %in% colnames(per_sample_df)) {
    stop(sprintf("value column '%s' not in data.frame", value))
  }
  if (!group_key %in% colnames(per_sample_df)) {
    stop(sprintf("group_key column '%s' not in data.frame", group_key))
  }
  if (!all(pair_keys %in% colnames(per_sample_df))) {
    missing_keys <- setdiff(pair_keys, colnames(per_sample_df))
    stop(sprintf("pair_keys not all present: %s", paste(missing_keys, collapse = ", ")))
  }
  # Use a `patient` column automatically when patients have several images.
  if (is.null(patient_key) && "patient" %in% colnames(per_sample_df) &&
      !all(is.na(per_sample_df$patient)) &&
      any(duplicated(per_sample_df[, intersect(c("patient", group_key, pair_keys), colnames(per_sample_df)), drop = FALSE]))) {
    patient_key <- "patient"
    message("compare_groups(): using the 'patient' column as patient_key (several images per patient).")
  }
  if (!is.null(patient_key) && !patient_key %in% colnames(per_sample_df)) {
    stop(sprintf("patient_key column '%s' not in data.frame", patient_key))
  }
  if (!is.null(covariates)) {
    if (method != "lmm") stop("`covariates` is only supported with method = 'lmm'.")
    miss <- setdiff(covariates, colnames(per_sample_df))
    if (length(miss)) stop(sprintf("covariates not in data.frame: %s", paste(miss, collapse = ", ")))
  }
  if (method == "lmm" && !requireNamespace("lme4", quietly = TRUE)) {
    stop("method = 'lmm' requires the 'lme4' package; install it or choose another method.")
  }
  if (method == "signrank" && is.null(patient_key)) {
    stop("method = 'signrank' requires patient_key (the pairing variable).")
  }
  if (method == "lmm" && is.null(patient_key)) {
    warning("method = 'lmm' without patient_key is equivalent to OLS; passing patient_key is recommended.")
  }

  # Seed locally: do not clobber the caller's RNG stream (e.g. inside a
  # bootstrap / simulation loop that calls compare_groups() repeatedly).
  .local_seed(seed)
  df <- as.data.frame(per_sample_df)
  df <- df[!is.na(df[[group_key]]), , drop = FALSE]

  groups_present <- unique(as.character(df[[group_key]]))
  if (length(groups_present) < 2L) {
    stop("Need at least 2 groups to compare. Found: ", paste(groups_present, collapse = ", "))
  }
  # Determine reference (g1) and test (g2) levels.
  if (!is.null(ref_group)) {
    if (!ref_group %in% groups_present) {
      stop(sprintf("ref_group '%s' not found in '%s'. Available: %s",
                   ref_group, group_key, paste(groups_present, collapse = ", ")))
    }
    g1 <- ref_group
    other <- sort(setdiff(groups_present, ref_group))
    if (length(other) > 1L) {
      warning(sprintf(
        "compare_groups() compares two groups; using ref '%s' vs '%s'. Subset the data for other contrasts.",
        g1, other[1]
      ))
    }
    g2 <- other[1]
  } else {
    groups_sorted <- sort(groups_present)
    if (length(groups_sorted) > 2L) {
      warning(sprintf(
        "compare_groups() compares two groups; using first 2 sorted levels (%s vs %s). Pass ref_group to control this.",
        groups_sorted[1], groups_sorted[2]
      ))
    }
    g1 <- groups_sorted[1]
    g2 <- groups_sorted[2]
    message(sprintf("ref_group not set: using '%s' as reference (effect = %s - %s).", g1, g2, g1))
  }
  df <- df[as.character(df[[group_key]]) %in% c(g1, g2), , drop = FALSE]

  if (symmetric && length(pair_keys) == 2L) {
    # One test per unordered pair: average the (i, j) and (j, i) rows of each
    # sample (degree-normalised scores are directional), keyed on the sorted
    # pair; other columns are taken from the (i <= j) row.
    a_ <- as.character(df[[pair_keys[1]]]); b_ <- as.character(df[[pair_keys[2]]])
    lo_ <- pmin(a_, b_); hi_ <- pmax(a_, b_)
    sid <- if ("sample_id" %in% colnames(df)) df$sample_id else seq_len(nrow(df))
    key <- paste(sid, as.character(df[[group_key]]), lo_, hi_, sep = "\r")
    avg <- tapply(df[[value]], key, function(x) mean(x, na.rm = TRUE))
    df <- df[a_ <= b_, , drop = FALSE]
    df[[value]] <- as.numeric(avg[key[a_ <= b_]])
  }

  # Patient as the unit: for tests that assume independent rows, average the
  # images of each patient (within group and pair) first.
  if (unit == "patient" && !is.null(patient_key) && method %in% c("wilcox", "t")) {
    keys <- c(patient_key, group_key, pair_keys)
    n_before <- nrow(df)
    df <- stats::aggregate(df[, value, drop = FALSE], by = df[, keys, drop = FALSE],
                           FUN = function(x) mean(x, na.rm = TRUE))
    if (nrow(df) < n_before) {
      message(sprintf("compare_groups(): averaged images within patients (unit = \"patient\"): %d rows -> %d.", n_before, nrow(df)))
    }
  }

  # Pseudoreplication check for tests that assume one row per unit.
  pat_col <- if (!is.null(patient_key)) patient_key else if ("patient" %in% colnames(df)) "patient" else NULL
  if (method %in% c("wilcox", "t") || (method == "perm" && is.null(patient_key))) {
    .warn_pseudoreplication(df, pair_keys, pat_col, method)
  }

  p_method <- switch(method,
    wilcox = "Wilcoxon rank-sum",
    signrank = "Wilcoxon signed-rank (paired by patient)",
    t = "Welch t-test",
    perm = "permutation",
    lmm = if (requireNamespace("lmerTest", quietly = TRUE)) {
      "LMM, Satterthwaite df (lmerTest)"
    } else {
      "LMM, t with df = n_patients - n_fixed"
    }
  )

  pair_id <- do.call(paste, c(df[pair_keys], sep = "\r"))
  chunks <- split(df, factor(pair_id, levels = unique(pair_id)))

  perm_null <- function(v, unit, unit_group) {
    # unit: per-row unit index; unit_group: group label per unit.
    n_u <- length(unit_group)
    n2 <- sum(unit_group == g2)
    stat <- function(lab) {
      rg <- lab[unit]
      mean(v[rg == g2]) - mean(v[rg == g1])
    }
    if (choose(n_u, n2) <= n_perms) {
      combs <- utils::combn(n_u, n2)
      vals <- apply(combs, 2, function(idx) {
        lab <- rep(g1, n_u); lab[idx] <- g2; stat(lab)
      })
      list(null = vals, exact = TRUE)
    } else {
      list(null = replicate(n_perms, stat(sample(unit_group))), exact = FALSE)
    }
  }

  fixed_terms <- paste(c(".g", if (length(covariates)) sprintf("`%s`", covariates)), collapse = " + ")

  results <- vector("list", length(chunks))
  for (i in seq_along(chunks)) {
    sub <- chunks[[i]]
    sub <- sub[is.finite(sub[[value]]), , drop = FALSE]
    g <- as.character(sub[[group_key]])
    if (sum(g == g1) < min_n_per_group || sum(g == g2) < min_n_per_group) next
    if (method == "signrank") {
      pm <- tapply(sub[[value]], list(as.character(sub[[patient_key]]), g), mean)
      pm <- pm[stats::complete.cases(pm[, c(g1, g2), drop = FALSE]), , drop = FALSE]
      if (nrow(pm) < min_n_per_group) next
    }

    res_row <- as.list(sub[1, pair_keys, drop = FALSE])
    res_row$n_total <- nrow(sub)
    res_row[[paste0("n_", g1)]] <- sum(g == g1)
    res_row[[paste0("n_", g2)]] <- sum(g == g2)
    v <- sub[[value]]
    m1 <- mean(v[g == g1])
    m2 <- mean(v[g == g2])
    res_row[[paste0("mean_", g1)]] <- m1
    res_row[[paste0("mean_", g2)]] <- m2
    res_row$effect <- m2 - m1

    test <- list(stat = NA_real_, p = NA_real_, est = NA_real_)
    fac <- factor(g, levels = c(g1, g2))
    if (method == "signrank") {
      d <- pm[, g2] - pm[, g1]
      res_row$n_pairs <- length(d)
      tt <- tryCatch(suppressWarnings(wilcox.test(d)), error = function(e) NULL)
      if (!is.null(tt)) { test$stat <- unname(tt$statistic); test$p <- tt$p.value }
    } else if (method == "wilcox") {
      # exact = NULL lets wilcox.test use the exact distribution for small
      # samples without ties (the normal approximation is anti-conservative
      # there, e.g. 3 vs 3: 0.081 vs exact 0.10).
      tt <- tryCatch(suppressWarnings(wilcox.test(v ~ fac)), error = function(e) NULL)
      if (!is.null(tt)) { test$stat <- unname(tt$statistic); test$p <- tt$p.value }
    } else if (method == "t") {
      tt <- tryCatch(t.test(v ~ fac), error = function(e) NULL)
      if (!is.null(tt)) { test$stat <- unname(tt$statistic); test$p <- tt$p.value }
    } else if (method == "lmm") {
      sub$.v <- v
      sub$.g <- fac
      coef_name <- paste0(".g", g2)
      use_lmm <- !is.null(patient_key) && length(unique(sub[[patient_key]])) < nrow(sub)
      if (use_lmm) {
        sub$.p <- as.character(sub[[patient_key]])
        form <- as.formula(sprintf(".v ~ %s + (1 | .p)", fixed_terms))
        has_lmertest <- requireNamespace("lmerTest", quietly = TRUE)
        fit <- tryCatch(suppressMessages(suppressWarnings(
          if (has_lmertest) lmerTest::lmer(form, data = sub, REML = TRUE)
          else lme4::lmer(form, data = sub, REML = TRUE)
        )), error = function(e) NULL)
        if (!is.null(fit)) {
          cf <- summary(fit)$coefficients
          if (coef_name %in% rownames(cf)) {
            test$est <- cf[coef_name, "Estimate"]
            test$stat <- cf[coef_name, "t value"]
            if (has_lmertest && "Pr(>|t|)" %in% colnames(cf)) {
              test$p <- cf[coef_name, "Pr(>|t|)"]
            } else {
              df_resid <- max(1, length(unique(sub$.p)) - nrow(cf))
              test$p <- 2 * stats::pt(-abs(test$stat), df = df_resid)
            }
          }
        }
      } else {
        fit <- tryCatch(lm(as.formula(sprintf(".v ~ %s", fixed_terms)), data = sub),
                        error = function(e) NULL)
        if (!is.null(fit)) {
          cf <- summary(fit)$coefficients
          if (coef_name %in% rownames(cf)) {
            test$est <- cf[coef_name, "Estimate"]
            test$stat <- cf[coef_name, "t value"]
            test$p <- cf[coef_name, "Pr(>|t|)"]
          }
        }
      }
    } else if (method == "perm") {
      observed <- m2 - m1
      paired <- FALSE
      if (!is.null(patient_key)) {
        pat <- as.character(sub[[patient_key]])
        paired <- any(tapply(g, pat, function(x) length(unique(x)) > 1))
      }
      pn <- if (paired) {
        .paired_perm_null(v, g, pat, g1, g2, n_perms)
      } else {
        if (!is.null(patient_key)) {
          pat_levels <- unique(pat)
          unit <- match(pat, pat_levels)
          unit_group <- g[match(pat_levels, pat)]
        } else {
          unit <- seq_along(g)
          unit_group <- g
        }
        perm_null(v, unit, unit_group)
      }
      tol <- sqrt(.Machine$double.eps)
      test$stat <- observed
      test$p <- if (pn$exact) {
        mean(abs(pn$null) >= abs(observed) - tol)
      } else {
        (sum(abs(pn$null) >= abs(observed) - tol, na.rm = TRUE) + 1) /
          (sum(!is.na(pn$null)) + 1)
      }
    }

    if (method == "lmm") res_row$estimate <- test$est
    res_row$statistic <- test$stat
    res_row$p <- test$p
    results[[i]] <- as.data.frame(res_row, stringsAsFactors = FALSE, check.names = FALSE)
  }
  results <- results[!vapply(results, is.null, logical(1))]
  if (length(results) == 0L) {
    warning("No comparable cluster pairs found.")
    return(invisible(NULL))
  }
  out <- do.call(rbind, results)
  out$padj <- p.adjust(out$p, method = adjust)
  out <- out[order(out$p), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "method") <- method
  attr(out, "groups") <- c(g1, g2)
  attr(out, "value") <- value
  attr(out, "p_method") <- p_method
  out
}

# ---- Phase 3: multi-group simulator -------------------------------------

#' Generate multi-sample simulated data with disease-group structure
#'
#' Wraps [generate_sim()] to produce N patients per group, each with a
#' (possibly noised) group-specific `close_ratio`, and optionally several
#' images (fields of view) per patient. Returns a single tidy data.frame
#' with `x`, `y`, `cell_type`, `sample_id`, `group`, `patient` columns —
#' directly consumable by [nhood_enrichment_per_sample()] and the other
#' `*_per_sample` helpers.
#'
#' @param n_samples_per_group Integer, number of patients to generate per
#'   group.
#' @param group_close_ratio Named list of base `close_ratio` values, one
#'   entry per group, e.g. `list(disease = 0.8, control = 0.2)`.
#' @param n_types,max_loc,test_type,distance_param Passed through
#'   to [generate_sim()].
#' @param n_cells Number of cells per image; either a single number or a
#'   named list with one entry per group (e.g. to simulate groups whose
#'   images differ in size). When a group-specific value is given and
#'   `max_loc` is a single number, `max_loc` is scaled by
#'   `sqrt(n_cells / n_cells_of_first_group)` so that cell density is kept
#'   constant.
#' @param between_sample_noise SD of Gaussian noise added to the
#'   per-patient `close_ratio` around the group baseline (clipped to [0,1]).
#' @param n_images_per_patient Number of images per patient.
#' @param within_patient_noise SD of Gaussian noise added to the per-image
#'   `close_ratio` around the patient value (only used when
#'   `n_images_per_patient > 1`).
#' @param seed Random seed (controls both the noise and the per-sample
#'   seeds passed to `generate_sim`).
#'
#' @return A data.frame. `sample_id` identifies an image; `patient`
#'   identifies the patient (equal to `sample_id` when
#'   `n_images_per_patient = 1`).
#' @export
#' @examples
#' df <- generate_sim_groups(n_samples_per_group = 2, n_images_per_patient = 2,
#'                           group_close_ratio = list(case = 0.8, control = 0.2),
#'                           n_types = 4, n_cells = 150, max_loc = 250,
#'                           test_type = "distribute", distance_param = 8)
#' unique(df[, c("sample_id", "patient", "group")])
generate_sim_groups <- function(n_samples_per_group = 3,
                                group_close_ratio = list(disease = 0.8, control = 0.2),
                                n_types = 10,
                                max_loc = 800,
                                n_cells = 1500,
                                test_type = "circle",
                                distance_param = 50,
                                between_sample_noise = 0.05,
                                n_images_per_patient = 1,
                                within_patient_noise = 0.05,
                                seed = 1234) {
  .seed_rng(seed)
  groups <- names(group_close_ratio)
  if (is.null(groups) || any(groups == "")) {
    stop("group_close_ratio must be a named list, e.g. list(disease = 0.8, control = 0.2)")
  }
  if (is.list(n_cells)) {
    if (!all(groups %in% names(n_cells))) stop("n_cells list must have one entry per group.")
    n_cells_g <- n_cells[groups]
  } else {
    n_cells_g <- stats::setNames(as.list(rep(n_cells, length(groups))), groups)
  }
  base_n <- n_cells_g[[1]]
  clip01 <- function(x) max(0, min(1, x))
  out <- list()
  for (g in groups) {
    base_ratio <- group_close_ratio[[g]]
    nc <- n_cells_g[[g]]
    ml <- if (is.list(n_cells)) max_loc * sqrt(nc / base_n) else max_loc
    for (s in seq_len(n_samples_per_group)) {
      pid <- paste0(g, "_", s)
      ratio_p <- clip01(base_ratio + rnorm(1, mean = 0, sd = between_sample_noise))
      for (m in seq_len(n_images_per_patient)) {
        if (n_images_per_patient == 1L) {
          sid <- pid
          ratio_s <- ratio_p
        } else {
          sid <- paste0(pid, "_img", m)
          ratio_s <- clip01(ratio_p + rnorm(1, mean = 0, sd = within_patient_noise))
        }
        df <- generate_sim(close_ratio = ratio_s,
                           n_types = n_types,
                           max_loc = ml,
                           n_cells = nc,
                           test_type = test_type,
                           distance_param = distance_param,
                           seed = sample.int(.Machine$integer.max, 1))
        df$sample_id <- sid
        df$group <- g
        df$patient <- pid
        rownames(df) <- paste0(sid, "_cell", seq_len(nrow(df)))
        out[[sid]] <- df
      }
    }
  }
  df_all <- do.call(rbind, out)
  rownames(df_all) <- unlist(lapply(out, rownames))
  df_all
}

#' Convert simulated cells to a Seurat object with one FOV per sample
#'
#' Build a minimal Seurat object from a coordinate table such as the output
#' of [generate_sim()] or [generate_sim_groups()]: every sample becomes a
#' centroid-based FOV named after it, and all other columns are stored in
#' `meta.data`. The expression matrix is a small placeholder, so the object
#' is meant for exercising the Seurat input path of the spatial functions,
#' not for expression analysis.
#'
#' @param df A data.frame with `x`, `y` and one row per cell.
#' @param sample_key Column identifying the sample (image). If absent, all
#'   cells are put in a single FOV named `"fov"`.
#' @param n_features Number of placeholder features in the count matrix.
#'
#' @return A Seurat object with a `cell` column in `meta.data` and one image
#'   per sample in `@images`.
#' @export
#' @examples
#' df <- generate_sim_groups(n_samples_per_group = 2, n_types = 4,
#'                           n_cells = 150, max_loc = 250,
#'                           test_type = "distribute", distance_param = 10)
#' seu <- sim_to_seurat(df)
#' SeuratObject::Images(seu)
sim_to_seurat <- function(df, sample_key = "sample_id", n_features = 5) {
  if (!all(c("x", "y") %in% colnames(df))) stop("`df` must contain 'x' and 'y' columns.")
  if (is.null(rownames(df)) || anyDuplicated(rownames(df))) {
    rownames(df) <- paste0("cell", seq_len(nrow(df)))
  }
  samples <- if (sample_key %in% colnames(df)) as.character(df[[sample_key]]) else rep("fov", nrow(df))
  counts <- Matrix::sparseMatrix(
    i = rep(seq_len(n_features), length.out = nrow(df)), j = seq_len(nrow(df)), x = 1,
    dims = c(n_features, nrow(df)),
    dimnames = list(paste0("feature", seq_len(n_features)), rownames(df))
  )
  md <- df[, setdiff(colnames(df), c("x", "y")), drop = FALSE]
  md$cell <- rownames(df)
  seu <- suppressWarnings(Seurat::CreateSeuratObject(counts, meta.data = md))
  for (s in unique(samples)) {
    sel <- samples == s
    cen <- SeuratObject::CreateCentroids(
      data.frame(x = df$x[sel], y = df$y[sel], cell = rownames(df)[sel])
    )
    fov <- SeuratObject::CreateFOV(list(centroids = cen), type = "centroids",
                                   assay = "RNA", key = paste0("fov", gsub("[^A-Za-z0-9]", "", s), "_"))
    seu[[s]] <- fov
  }
  seu
}

