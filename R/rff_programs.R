# One-call pipeline for spatial gene programs: offset -> fit -> program test ->
# merged "Programs" table. The matching of detection axes (PCs of
# rff_program_test()) to program maps (factors of fit_spatial_rff()) lives in
# .rff_match_programs(), which rff_report() uses as well.

#' Find spatial gene programs (recommended entry point)
#'
#' **Experimental.** Runs the recommended residual-RFLVM pipeline in one call
#' and returns one table of **programs**:
#' 1. known structure as offset ([rff_offset()], from a shortcut, covariates
#'    or a ready matrix);
#' 2. the random-feature factor model ([fit_spatial_rff()]; grid basis, fixed
#'    length scale, residual-PCA initialisation, ARD);
#' 3. the program test ([rff_program_test()], gene-shift null by default,
#'    sequential), optionally also [rff_factor_test()];
#' 4. matching of the two kinds of components into one Programs table.
#'
#' @section Which function should I use?: Most users call `rff_programs()`
#'   and then [rff_report()] (or `report = "file.html"`). [fit_spatial_rff()],
#'   [rff_program_test()], [rff_factor_test()] and [rff_offset()] are the
#'   building blocks, for custom pipelines.
#'
#' @section Detection axes and program maps: The program test finds
#'   **detection axes** (principal components of smoothed residuals, `PC1`,
#'   ...) and decides whether and how many programs exist; the fit gives
#'   **program maps** (factors, `factor1`, ...): a smooth field and gene
#'   loadings used for interpretation. Each significant detection axis is
#'   matched to the program map whose field correlates best with its score
#'   (greedy by |r|, every map used once):
#'   * `Confirmed`: significant axis matched with |r| >= `min_r`; p from the
#'     axis, genes, length scale and map from the program map.
#'   * `Candidate`: significant axis without such a map - the RFLVM does not
#'     represent it well; genes and map come from the axis itself.
#'   * `Not supported`: a strong program map (large program strength, or
#'     significant in the factor test) without a significant matching axis.
#'   * `Cellularity/technical`: program maps whose loadings are mostly shared
#'     by all genes (uniform share > 0.5).
#'   With only the factor test (`null = "none"`, `factor_test = TRUE`),
#'   programs are the significant program maps (`detected_by = "factor
#'   test"`); without any test they are `Exploratory`.
#'
#' @param binned A `binned_transcripts` object.
#' @param offset Known structure: `NULL` / `"area"` (none), `"kmeans"` or
#'   `"pca"` ([rff_offset()] with `method`, `k`), `"smoothed_total"` (the log
#'   of the Gaussian-smoothed total transcript density, bandwidth
#'   `offset_bandwidth`, as covariate), a data frame of covariates (one row
#'   per bin or in-tissue bin; [rff_offset()] with `method = "covariates"`),
#'   or a matrix from [rff_offset()] (bins x genes).
#' @param fit Optional existing [fit_spatial_rff()] fit to reuse (then
#'   `offset`, `n_factors`, `lengthscale`, `ard` and `fit_args` are not used
#'   for fitting; `offset` only documents how the fit's offset was made).
#' @param n_factors,lengthscale,ard Fit settings (see [fit_spatial_rff()]);
#'   `lengthscale = "profile"` estimates one length scale per program by
#'   held-out likelihood ([rff_lengthscale_profile()]).
#' @param fit_args Further arguments to [fit_spatial_rff()] (they override
#'   the defaults above, e.g. `list(max_iter = 300)`).
#' @param null Null of [rff_program_test()]: `"shift"` (recommended on real
#'   data), `"parametric"`, or `"none"` (no program test).
#' @param highpass,bandwidth Passed to [rff_program_test()].
#' @param n_boot Null data sets of the program test.
#' @param alpha Significance level.
#' @param factor_test Also run [rff_factor_test()] (sequential), with
#'   `factor_n_boot` refits.
#' @param factor_n_boot Refits for the factor test.
#' @param min_r Minimum |r| between a detection axis score and a program map
#'   field for the two to be matched.
#' @param k Number of clusters / components for `offset = "kmeans"` / `"pca"`.
#' @param offset_bandwidth Bandwidth for `offset = "smoothed_total"`.
#' @param seed Random seed.
#' @param n_cores Cores for the factor test refits.
#' @param report Optional path: write [rff_report()] of the result there.
#' @param crossfit,control,min_effect Calibrated program test (see
#'   [rff_program_test()]): cross-fitted held-out effect sizes, comparison
#'   with negative-control windows ([rff_control_reference()]) and an
#'   effect-size threshold. With any of them, a detection axis counts as
#'   significant only when its `call` is `TRUE`. Recommended on real
#'   tissue, where the plain gene-shift test is rejected by almost any
#'   shared residual structure.
#' @param ls_grid,profile With `lengthscale = "profile"`: length scales to
#'   profile and settings for [rff_lengthscale_profile()] (see
#'   [fit_spatial_rff()]); the programs table then reports each program
#'   map's profiled length scale with its 95% interval.
#' @param ... Passed to [rff_report()] when `report` is given.
#'
#' @return An object of class `rff_programs`: a list with `programs` (data
#'   frame, one row per program; see Details), `fit`, `program_test`,
#'   `factor_test`, `offset_info` (how the offset was made), `settings`,
#'   `call` and `elapsed` (seconds). The programs table has columns
#'   `program` (P1, P2, ... for Confirmed / Candidate / Exploratory
#'   programs), `status`, `detected_by`, `p`, `p_factor_test`,
#'   `detection_axis` (PC of the program test), `program_map` (factor of the
#'   fit), `r` (|r| between the two), `best_map` / `best_r` (best map of an
#'   axis even if not matched), `lengthscale`, `program_strength`,
#'   `uniform_share`, `variance_share`, `n_genes` (genes with |loading| >= 25%
#'   of the largest), `top_up`, `top_down`, `key` (unique row name) and
#'   `map_label` / `axis_label` (M* / D* names of the report),
#'   `lengthscale_lower` / `lengthscale_upper` / `heldout_dev_explained`
#'   (profiled length scales), `excess` / `excess_lower` / `excess_upper` /
#'   `p_control` (calibrated test; `NA` otherwise). Attributes
#'   `loadings` (genes x rows, the loadings used for interpretation) and
#'   `fields` (all bins x rows: the program map field, or the axis score for
#'   Candidates; `NA` outside the tissue). Per-cell values: use
#'   `rff_fields(res$fit, cells)` and the `program_map` column (Confirmed
#'   programs), see Examples.
#' @seealso [rff_report()], [fit_spatial_rff()], [rff_program_test()],
#'   [rff_factor_test()], [rff_offset()]
#' @export
#' @examples
#' \donttest{
#' tx <- simulate_transcripts(size = 100, rate = 0.02, n_genes_per_set = 3)
#' b <- bin_transcripts(tx, bin_size = 5)
#' res <- rff_programs(b, n_factors = 3, lengthscale = 10, n_boot = 19,
#'                     fit_args = list(max_iter = 60))
#' res
#' as.data.frame(res)
#' # per-cell values of the confirmed programs (via their program maps)
#' cells <- data.frame(x = runif(50, 0, 100), y = runif(50, 0, 100))
#' pm <- res$programs$program_map[res$programs$status == "Confirmed"]
#' if (length(pm)) head(rff_fields(res$fit, cells)[, pm, drop = FALSE])
#' }
rff_programs <- function(binned, offset = NULL, fit = NULL, n_factors = 6, lengthscale = 8, ard = 100,
                         fit_args = list(), null = c("shift", "parametric", "none"), highpass = NULL,
                         bandwidth = NULL, n_boot = 99, alpha = 0.05, factor_test = FALSE, factor_n_boot = 19,
                         min_r = 0.5, k = 6, offset_bandwidth = 10, seed = 1, n_cores = 1, report = NULL,
                         crossfit = FALSE, control = NULL, min_effect = NULL, ls_grid = 5 * 2^(0:6), profile = list(), ...) {
  t0 <- Sys.time()
  cl <- match.call()
  null <- match.arg(null)
  if (!inherits(binned, "binned_transcripts")) stop("`binned` must come from bin_transcripts().")
  keep <- binned$coords$in_tissue
  # ---- offset ----
  off <- .rff_programs_offset(binned, offset, k, seed, offset_bandwidth, cl$offset)
  # ---- fit ----
  if (is.null(fit)) {
    args <- utils::modifyList(list(binned = binned, offset = off$offset, n_factors = n_factors, basis = "grid",
                                   lengthscales = lengthscale, learn_lengthscales = FALSE,
                                   factor_init = "residual_pca", ard = ard, seed = seed), fit_args)
    if (identical(lengthscale, "profile")) {
      args$ls_grid <- ls_grid
      args$profile <- utils::modifyList(list(n_cores = n_cores, seed = seed), profile)
    }
    fit <- do.call(fit_spatial_rff, args)
    fit_call <- as.call(c(as.name("fit_spatial_rff"), lapply(args[setdiff(names(args), "binned")], function(a) if (is.matrix(a)) as.name("off") else a)))
    fit_call$binned <- as.name("binned")
  } else {
    if (!inherits(fit, "spatial_rff_fit")) stop("`fit` must come from fit_spatial_rff().")
    st <- .rr_or(fit$settings, list())
    msg <- c(if (!identical(st$basis, "grid")) "basis = \"grid\"",
             if (!identical(st$learn_lengthscales, FALSE)) "learn_lengthscales = FALSE",
             if (!isTRUE(st$ard > 0)) "ard > 0")
    if (length(msg)) warning("The supplied fit differs from the recommended settings (", paste(msg, collapse = ", "), ").")
    if (is.null(offset)) off$info <- list(method = if (is.null(fit$offset_matrix)) .rr_or(fit$offset, "area") else NA)
    fit_call <- NULL
  }
  # ---- tests ----
  pt <- ft <- NULL
  if (null != "none")
    pt <- rff_program_test(fit, binned, bandwidth = bandwidth, n_boot = n_boot, sequential = TRUE, alpha = alpha,
                           seed = seed, null = null, highpass = if (null == "shift") highpass else NULL,
                           crossfit = crossfit, control = control, min_effect = min_effect)
  if (!is.null(pt) && !is.null(highpass)) attr(pt, "highpass") <- highpass
  if (!is.null(pt) && !is.null(bandwidth)) attr(pt, "bandwidth") <- bandwidth
  if (factor_test) ft <- rff_factor_test(fit, binned, n_boot = factor_n_boot, seed = seed, n_cores = n_cores,
                                         sequential = TRUE, alpha = alpha)
  progs <- .rff_match_programs(fit, pt, ft, alpha = alpha, min_r = min_r)
  res <- structure(list(
    programs = progs, fit = fit, program_test = pt, factor_test = ft, offset_info = off$info,
    settings = list(n_factors = ncol(fit$L), lengthscale = lengthscale, ard = ard, null = null, highpass = highpass,
                    bandwidth = bandwidth, n_boot = n_boot, alpha = alpha, factor_test = factor_test,
                    factor_n_boot = factor_n_boot, min_r = min_r, seed = seed,
                    crossfit = crossfit || !is.null(control) || !is.null(min_effect), control = !is.null(control),
                    min_effect = min_effect),
    fit_call = fit_call, call = cl,
    elapsed = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ), class = "rff_programs")
  if (!is.null(report)) res$report <- rff_report(res, binned, file = report, ...)
  res
}

# offset shortcuts -> list(offset = <"area" | matrix>, info = <offset_info>)
.rff_programs_offset <- function(binned, offset, k, seed, bw, expr) {
  keep <- binned$coords$in_tissue
  dep <- if (is.null(expr)) NULL else paste(deparse(expr, width.cutoff = 200), collapse = " ")
  if (is.null(offset) || identical(offset, "area"))
    return(list(offset = "area", info = list(method = "area")))
  if (is.character(offset) && length(offset) == 1) {
    if (offset %in% c("kmeans", "pca"))
      return(list(offset = rff_offset(binned, method = offset, k = k, seed = seed),
                  info = list(method = offset, k = k, call = sprintf("rff_offset(binned, method = \"%s\", k = %d, seed = %d)", offset, k, seed))))
    if (offset == "smoothed_total") {
      g <- binned$grid
      tot <- matrix(Matrix::rowSums(binned$counts), g$nx, g$ny); msk <- matrix(as.numeric(keep), g$nx, g$ny)
      sm <- .smooth_grid(tot * msk, bw / g$bin_size) / pmax(.smooth_grid(msk, bw / g$bin_size), 1e-8)
      cv <- data.frame(log_smoothed_total = log(as.vector(sm)[keep] + 0.01 * mean(sm[keep])))
      return(list(offset = rff_offset(binned, cv),
                  info = list(method = "covariates", covariates = cv,
                              note = sprintf("log of the total transcript density smoothed with bandwidth %g", bw),
                              call = "rff_offset(binned, data.frame(log_smoothed_total = ...))")))
    }
    stop("Unknown `offset` shortcut \"", offset, "\"; use \"area\", \"kmeans\", \"pca\" or \"smoothed_total\".")
  }
  if (is.matrix(offset) && is.numeric(offset) && !is.null(colnames(offset)) && all(binned$genes %in% colnames(offset)))
    return(list(offset = offset, info = list(method = NA, call = dep)))
  if (is.data.frame(offset) || is.matrix(offset)) {
    cv <- as.data.frame(offset)
    return(list(offset = rff_offset(binned, cv),
                info = list(method = "covariates", covariates = if (nrow(cv) == nrow(binned$coords)) cv[keep, , drop = FALSE] else cv,
                            call = if (!is.null(dep)) sprintf("rff_offset(binned, %s)", dep))))
  }
  stop("`offset` must be NULL, a shortcut string, a covariate data frame or an offset matrix from rff_offset().")
}

# Match detection axes (program test PCs) to program maps (fit factors) and
# build the Programs table. Used by rff_programs() and rff_report().
.rff_match_programs <- function(fit, program_test = NULL, factor_test = NULL, alpha = 0.05, min_r = 0.5,
                                top_n = 6) {
  L <- fit$L; genes <- .rr_or(fit$genes, rownames(L)); rownames(L) <- genes
  fn <- colnames(L); K <- length(fn)
  keep <- fit$in_tissue
  strength <- .rr_or(fit$factor_strength, stats::setNames(sqrt(colSums(L^2)), fn))[fn]
  pstr <- .rr_or(fit$program_strength, stats::setNames(sqrt(colSums(sweep(L, 2, colMeans(L))^2)), fn))[fn]
  uni <- 1 - pstr^2 / pmax(strength^2, 1e-12)
  ell <- stats::setNames(rep_len(.rr_or(fit$lengthscales, NA_real_), K), fn)
  mlab <- stats::setNames(paste0("M", seq_len(K)), fn)
  Fg <- fit$field_grid[, intersect(fn, colnames(fit$field_grid)), drop = FALSE]
  fp <- stats::setNames(rep(NA_real_, K), fn)
  if (!is.null(factor_test)) {
    pu <- if ("p_sequential" %in% names(factor_test)) factor_test$p_sequential else factor_test$p
    fp[as.character(factor_test$factor)] <- pu
  }
  strong <- strength >= 0.1 * max(strength)
  cellular <- uni > 0.5 & strong
  rows <- list(); loads <- list(); flds <- list()
  top_str <- function(v, sgn) { v <- v[is.finite(v)]; v <- if (sgn > 0) v[v > 0] else v[v < 0]
    paste(names(v)[order(-abs(v))][seq_len(min(top_n, length(v)))], collapse = ", ") }
  ngen <- function(v) { a <- abs(v[is.finite(v)]); if (!length(a) || max(a) == 0) 0L else sum(a >= 0.25 * max(a)) }
  lsp <- fit$lengthscale_profile$summary
  eff_of <- function(axis, col) {
    if (is.null(program_test) || is.na(axis) || !(col %in% names(program_test))) return(NA_real_)
    program_test[[col]][match(axis, as.character(program_test$component))]
  }
  ls_of <- function(map, col) if (is.null(lsp) || is.na(map) || !(map %in% lsp$factor)) NA_real_ else lsp[[col]][match(map, lsp$factor)]
  add <- function(status, detected_by, p, axis, axlab, map, r, best_map, best_r, share, v, field) {
    key <- paste0("row", length(rows) + 1)
    rows[[length(rows) + 1]] <<- data.frame(
      program = NA_character_, status = status, detected_by = detected_by, p = p,
      p_factor_test = if (!is.na(map)) fp[[map]] else NA_real_,
      detection_axis = axis, program_map = map, r = r, best_map = best_map, best_r = best_r,
      lengthscale = if (!is.na(map)) ell[[map]] else NA_real_,
      program_strength = if (!is.na(map)) pstr[[map]] else NA_real_,
      uniform_share = if (!is.na(map)) uni[[map]] else NA_real_,
      variance_share = share, n_genes = ngen(v), top_up = top_str(v, 1), top_down = top_str(v, -1),
      key = key, map_label = if (!is.na(map)) mlab[[map]] else NA_character_, axis_label = axlab,
      lengthscale_lower = ls_of(map, "lower"), lengthscale_upper = ls_of(map, "upper"),
      heldout_dev_explained = ls_of(map, "heldout_dev_explained"),
      excess = eff_of(axis, "excess"), excess_lower = eff_of(axis, "excess_lower"), excess_upper = eff_of(axis, "excess_upper"),
      p_control = eff_of(axis, "p_control"),
      stringsAsFactors = FALSE)
    loads[[key]] <<- v; flds[[key]] <<- field
  }
  used <- character(0)
  if (!is.null(program_test)) {
    comp <- as.character(program_test$component)
    U <- attr(program_test, "scores", exact = TRUE); V <- attr(program_test, "loadings", exact = TRUE)
    # calibrated tests (crossfit / control / min_effect) carry their decision in `call`
    sigc <- if ("call" %in% names(program_test)) comp[program_test$call %in% TRUE] else
      comp[!is.na(program_test$p) & program_test$p <= alpha]
    R <- matrix(0, length(comp), ncol(Fg), dimnames = list(comp, colnames(Fg)))
    if (!is.null(U) && nrow(U) == sum(keep)) {
      R <- suppressWarnings(abs(stats::cor(U[, comp, drop = FALSE], Fg[keep, , drop = FALSE])))
      R[!is.finite(R)] <- 0; dimnames(R) <- list(comp, colnames(Fg))
    }
    # greedy matching of significant axes to maps
    match_of <- stats::setNames(rep(NA_character_, length(sigc)), sigc)
    if (length(sigc)) {
      pr <- expand.grid(a = sigc, m = colnames(Fg), stringsAsFactors = FALSE)
      pr$r <- R[cbind(pr$a, pr$m)]
      pr <- pr[pr$r >= min_r, , drop = FALSE]; pr <- pr[order(-pr$r), , drop = FALSE]
      for (i in seq_len(nrow(pr))) if (is.na(match_of[[pr$a[i]]]) && !(pr$m[i] %in% used)) {
        match_of[[pr$a[i]]] <- pr$m[i]; used <- c(used, pr$m[i])
      }
    }
    for (a in sigc) {
      k <- match(a, comp); pa <- program_test$p[k]; bm <- colnames(R)[which.max(R[a, ])]; br <- max(R[a, ])
      axlab <- paste0("D", sub("^PC", "", a))
      full <- rep(NA_real_, length(keep)); if (!is.null(U)) full[keep] <- U[, a]
      m <- match_of[[a]]
      if (!is.na(m)) {
        v <- L[, m]
        # orient the map like the axis for display consistency is not needed: fields keep their sign
        add("Confirmed", "program test", pa, a, axlab, m, R[a, m], bm, br, program_test$share[k], v, Fg[, m])
      } else {
        v <- if (!is.null(V)) stats::setNames(V[genes, a], genes) else stats::setNames(rep(NA_real_, length(genes)), genes)
        add("Candidate", "program test", pa, a, axlab, NA_character_, NA_real_, bm, br, program_test$share[k], v, full)
      }
    }
    # maps without a significant matching axis
    cand <- if (!is.null(factor_test)) fn[!is.na(fp) & fp <= alpha] else fn[strong & pstr >= 0.5 * max(pstr)]
    for (m in setdiff(cand, c(used, fn[cellular]))) {
      bestax <- if (nrow(R)) rownames(R)[which.max(R[, m])] else NA
      add("Not supported", if (!is.null(factor_test)) "factor test" else "program strength", fp[[m]], NA_character_, NA_character_,
          m, NA_real_, bestax, if (nrow(R)) max(R[, m]) else NA_real_, NA_real_, L[, m], Fg[, m])
    }
  } else if (!is.null(factor_test)) {
    for (m in fn[!is.na(fp) & fp <= alpha & !cellular]) {
      add("Confirmed", "factor test", fp[[m]], NA_character_, NA_character_, m, NA_real_, NA_character_, NA_real_, NA_real_, L[, m], Fg[, m])
      used <- c(used, m)
    }
  } else {
    for (m in fn[strong & !cellular]) {
      add("Exploratory", "none (no test)", NA_real_, NA_character_, NA_character_, m, NA_real_, NA_character_, NA_real_, NA_real_, L[, m], Fg[, m])
      used <- c(used, m)
    }
  }
  for (m in setdiff(fn[cellular], used))
    add("Cellularity/technical", "uniform share > 0.5", fp[[m]], NA_character_, NA_character_, m, NA_real_, NA_character_, NA_real_, NA_real_, L[, m], Fg[, m])
  cols <- c("program", "status", "detected_by", "p", "p_factor_test", "detection_axis", "program_map", "r", "best_map",
            "best_r", "lengthscale", "program_strength", "uniform_share", "variance_share", "n_genes", "top_up", "top_down",
            "key", "map_label", "axis_label", "lengthscale_lower", "lengthscale_upper", "heldout_dev_explained",
            "excess", "excess_lower", "excess_upper", "p_control")
  if (!length(rows)) {
    out <- as.data.frame(stats::setNames(replicate(length(cols), logical(0), simplify = FALSE), cols))
    attr(out, "loadings") <- matrix(numeric(0), length(genes), 0, dimnames = list(genes, NULL))
    attr(out, "fields") <- matrix(numeric(0), length(keep), 0)
    return(out)
  }
  out <- do.call(rbind, rows)
  rk <- match(out$status, c("Confirmed", "Candidate", "Exploratory", "Not supported", "Cellularity/technical"))
  o <- order(rk, ifelse(is.na(out$p), Inf, out$p), -ifelse(is.na(out$r), 0, out$r), -ifelse(is.na(out$program_strength), 0, out$program_strength))
  out <- out[o, , drop = FALSE]
  main <- out$status %in% c("Confirmed", "Candidate", "Exploratory")
  out$program[main] <- paste0("P", seq_len(sum(main)))
  out$key[main] <- out$program[main]
  out$key[!main] <- ifelse(!is.na(out$map_label[!main]), out$map_label[!main], paste0("X", seq_len(sum(!main))))
  old <- vapply(rows, `[[`, "", "key")[o]
  rownames(out) <- NULL
  Lm <- do.call(cbind, loads[old]); colnames(Lm) <- out$key; rownames(Lm) <- genes
  Fm <- do.call(cbind, flds[old]); colnames(Fm) <- out$key
  attr(out, "loadings") <- Lm
  attr(out, "fields") <- Fm
  out
}

#' @export
print.rff_programs <- function(x, ...) {
  p <- x$programs
  n <- table(factor(p$status, c("Confirmed", "Candidate", "Exploratory", "Not supported", "Cellularity/technical")))
  test <- if (!is.null(x$program_test)) sprintf("program test, %s null, %d draws", .rr_or(attr(x$program_test, "null", exact = TRUE), "parametric"),
                                                 length(attr(x$program_test, "null_shares", exact = TRUE)[[1]])) else
    if (!is.null(x$factor_test)) "factor test only" else "no test (exploratory)"
  cat(sprintf("<rff_programs> %d confirmed, %d candidate%s program(s); %s; alpha %g; %.0f s\n",
              n[["Confirmed"]], n[["Candidate"]], if (n[["Exploratory"]]) sprintf(", %d exploratory", n[["Exploratory"]]) else "",
              test, x$settings$alpha, x$elapsed))
  main <- p[!is.na(p$program), , drop = FALSE]
  for (i in seq_len(nrow(main))) {
    r <- main[i, ]
    src <- paste(c(if (!is.na(r$detection_axis)) sprintf("axis %s", r$detection_axis),
                   if (!is.na(r$program_map)) sprintf("map %s%s", r$program_map, if (!is.na(r$r)) sprintf(", r %.2f", r$r) else "")), collapse = "; ")
    cat(sprintf("  %-3s %-11s p %-6s [%s]  up: %s%s\n", r$program, r$status, if (is.na(r$p)) "-" else formatC(r$p, format = "g", digits = 2),
                src, r$top_up, if (nzchar(r$top_down)) paste0("  down: ", r$top_down) else ""))
  }
  if (n[["Candidate"]]) cat("  ! Candidate: significant detection axis without a matching program map (|r| <", x$settings$min_r, "); read its own loadings.\n")
  if (n[["Not supported"]]) cat(sprintf("  ! Not supported (strong program maps without a significant axis): %s\n",
                                        paste(p$program_map[p$status == "Not supported"], collapse = ", ")))
  if (n[["Cellularity/technical"]]) cat(sprintf("  ! Cellularity-like program maps: %s\n", paste(p$program_map[p$status == "Cellularity/technical"], collapse = ", ")))
  if (is.null(x$program_test) && is.null(x$factor_test)) cat("  ! No test: exploratory maps, no significance.\n")
  invisible(x)
}

#' @export
summary.rff_programs <- function(object, ...) {
  p <- object$programs
  out <- list(counts = table(factor(p$status, c("Confirmed", "Candidate", "Exploratory", "Not supported", "Cellularity/technical"))),
              programs = p[, c("program", "status", "p", "detection_axis", "program_map", "r", "lengthscale", "n_genes", "top_up", "top_down")],
              offset = object$offset_info$method, settings = object$settings, convergence = object$fit$convergence)
  class(out) <- "summary.rff_programs"
  out
}

#' @export
print.summary.rff_programs <- function(x, ...) {
  cat("Programs by status:\n"); print(x$counts)
  cat("\nOffset:", format(.rr_or(x$offset, "not recorded")), " | fit convergence:", format(x$convergence), "\n\n")
  print(x$programs, row.names = FALSE)
  invisible(x)
}

#' @export
as.data.frame.rff_programs <- function(x, ...) {
  p <- x$programs; attr(p, "loadings") <- NULL; attr(p, "fields") <- NULL; p
}
