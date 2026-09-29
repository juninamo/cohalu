# Interactive HTML report for a (residual) random-feature LGCP analysis.
# Everything is written by base R: JSON, base64 and HTML are built here, so no
# dependency is added. The page template (CSS + JS) is inst/rff_report/template.html.

#' Interactive HTML report of a residual RFLVM analysis
#'
#' **Experimental.** Writes a single self-contained HTML file (no external
#' requests; data embedded, maps drawn on `<canvas>`) that summarises a fit of
#' [fit_spatial_rff()], optionally with the results of [rff_factor_test()]
#' and / or [rff_program_test()], per-cell factor values ([rff_fields()]) and
#' the overlap of the factors with known gene programs. Machine-readable
#' outputs (CSV tables, the fit as RDS, the settings and an R script to
#' reproduce the analysis) are written next to it and listed in the report.
#'
#' The report has, in order: an overview; **Input data** (bins, genes,
#' counts per bin, the offset / known structure and how much of each gene's
#' deviance it explains, cells, known programs, all fit and test arguments);
#' **Key findings** (plain-language bullets with ok / caution / warning
#' badges: significant programs, their genes, length scales, cellularity-like
#' factors, cell types, known-program overlap, and warnings such as a fit
#' stopped at `max_iter`, calibration caveats of the null used, few bins, an
#' irregular tissue mask, factors shrunk by ARD or dominated by one or two
#' genes); **Output files**; sortable tables of the program maps and
#' detection axes - clicking a row selects it everywhere; spatial maps with pan / zoom
#' (drag, wheel; double-click resets), hover values and a scale bar, small
#' multiples of every field, the cellularity field, total counts, the offset
#' and, for selected genes, observed counts, the null expectation (offset,
#' intercept and cellularity, without factors), their log ratio and the
#' factor contribution \eqn{\sum_k L_{jk} f_k}; a loading heatmap, ranked
#' gene bars and a gene search; per-cell summaries by label; the
#' known-program table; null distributions of the tests; methods; how to
#' reproduce; session information.
#'
#' @param fit Output of [fit_spatial_rff()] or [rff_programs()], or a
#'   multi-window fit of [fit_spatial_rff_joint()] (see `window`).
#' @param binned The `binned_transcripts` object used for the fit (a list of
#'   them for a multi-window fit).
#' @param test `NULL`, a data frame returned by [rff_factor_test()] or
#'   [rff_program_test()], or a (named) list of such data frames.
#' @param file Path of the HTML file.
#' @param title Report title (default: derived from the fit).
#' @param top_genes Number of top genes per program shown in the ranked bar
#'   lists (the heatmap shows up to 8 per program).
#' @param cells Optional data frame of cells (centroids) with coordinates in
#'   the frame of the transcripts, and optionally a label column.
#' @param cell_label_col Column of `cells` with cell-type (or other) labels.
#' @param cell_xy Coordinate columns of `cells`.
#' @param known_programs Optional named list of character vectors (gene
#'   sets); each factor / component is compared with each set (correlation of
#'   the loadings with set membership, overlap of the top genes).
#' @param offset_info Optional list describing the offset, e.g.
#'   `list(method = "kmeans", k = 6)` or `list(method = "covariates",
#'   covariates = cv, call = "rff_offset(b, cv)")`. The fit does not record
#'   how a matrix offset was made; without this the report says so.
#' @param meta Optional named list of sample / section identifiers and other
#'   notes (e.g. `list(sample = "S1", section = "A", window = "...")`), shown
#'   in the input section.
#' @param map_genes Genes whose maps are embedded (in addition to the top
#'   genes of each program, up to `n_map_genes` in total).
#' @param n_map_genes Maximum number of embedded gene maps.
#' @param max_pixels Maximum number of grid cells per embedded map; larger
#'   grids are block-averaged for display (with a warning).
#' @param max_cells Maximum number of cells drawn as points (summaries use
#'   all cells).
#' @param alpha Significance level used for flags and findings.
#' @param export Write the machine-readable outputs to `out_dir`.
#' @param out_dir Directory for the outputs (default: `<file without
#'   extension>_files` next to the HTML).
#' @param open Open the report in the browser ([utils::browseURL()]).
#' @param window Only for a multi-window fit ([fit_spatial_rff_joint()]
#'   given as `fit`, with `binned` the list of its windows): the window
#'   (name or index) whose fields are shown; default the window with the
#'   largest total program amplitude. `test` may then be the discovery-mode
#'   output of [rff_program_test_joint()] (p-values of the joint test,
#'   variance shares and scores of this window); the shared loadings and the
#'   per-window amplitudes are summarised in the input section.
#'
#' @section Names used in the report: The factors of the fit are shown as
#'   **program maps** M1, M2, ... (the smooth field and gene loadings used for
#'   interpretation) and the PCs of [rff_program_test()] as **detection axes**
#'   D1, D2, ... (which decide whether and how many programs exist). Object
#'   and column names are unchanged (`factor1`, `PC1`, ...); exported tables
#'   carry a `display_name` column and `display_names.csv` maps the two.
#'
#' @return The path of the HTML file, invisibly.
#' @seealso [fit_spatial_rff()], [rff_factor_test()], [rff_program_test()],
#'   [rff_fields()], [rff_offset()]
#' @export
#' @examples
#' \donttest{
#' tx <- simulate_transcripts(size = 100, rate = 0.02, n_genes_per_set = 3)
#' b <- bin_transcripts(tx, bin_size = 6)
#' fit <- fit_spatial_rff(b, n_factors = 2, ard = 2, n_features = 24, max_iter = 40)
#' pt <- rff_program_test(fit, b, n_boot = 19, null = "shift")
#' f <- rff_report(fit, b, test = pt, file = file.path(tempdir(), "rflvm_report.html"),
#'                 known_programs = list(A = c("A_1", "A_2", "A_3")), open = FALSE)
#' }
rff_report <- function(fit, binned, test = NULL, file = "rflvm_report.html", title = NULL,
                       top_genes = 15, cells = NULL, cell_label_col = NULL, cell_xy = c("x", "y"),
                       known_programs = NULL, offset_info = NULL, meta = NULL,
                       map_genes = NULL, n_map_genes = 30, max_pixels = 40000, max_cells = 20000,
                       alpha = 0.05, min_r = 0.5, export = TRUE, out_dir = NULL, open = interactive(),
                       window = NULL) {
  t_start <- Sys.time()
  cl <- match.call()
  res_obj <- NULL
  if (inherits(fit, "spatial_rff_joint")) {
    # multi-window fit with shared loadings: report one window (fields of that window, shared loadings)
    jo <- .rr_joint_window(fit, binned, test, window)
    fit <- jo$fit; binned <- jo$binned; test <- jo$test
    meta <- c(meta, jo$meta)
  }
  if (inherits(fit, "rff_programs")) {
    res_obj <- fit; fit <- res_obj$fit
    if (is.null(test)) {
      test <- list(`program test` = res_obj$program_test, `factor test` = res_obj$factor_test)
      test <- test[!vapply(test, is.null, TRUE)]
      if (!length(test)) test <- NULL
    }
    if (is.null(offset_info)) offset_info <- res_obj$offset_info
    if (missing(alpha)) alpha <- .rr_or(res_obj$settings$alpha, alpha)
    if (missing(min_r)) min_r <- .rr_or(res_obj$settings$min_r, min_r)
  }
  if (!inherits(fit, "spatial_rff_fit")) stop("`fit` must come from fit_spatial_rff() or rff_programs().")
  if (!inherits(binned, "binned_transcripts")) stop("`binned` must come from bin_transcripts().")
  g <- .rr_or(fit$grid, binned$grid)
  if (nrow(binned$coords) != g$nx * g$ny) stop("`binned` does not match the grid of `fit`.")
  if (!is.null(known_programs) && (!is.list(known_programs) || is.null(names(known_programs))))
    stop("`known_programs` must be a named list of gene sets.")
  keep <- .rr_or(fit$in_tissue, binned$coords$in_tissue)
  genes <- .rr_or(fit$genes, rownames(fit$L))
  L <- fit$L
  K <- ncol(L); fnames <- colnames(L); if (is.null(fnames)) fnames <- colnames(L) <- paste0("factor", seq_len(K))
  if (is.null(rownames(L))) rownames(L) <- genes
  st <- .rr_or(fit$settings, list())
  bs <- g$bin_size
  N <- sum(keep); J <- length(genes)
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  file <- file.path(normalizePath(dirname(file)), basename(file))
  if (is.null(out_dir)) out_dir <- paste0(sub("\\.html?$", "", file, ignore.case = TRUE), "_files")
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  if (is.null(title)) title <- sprintf("Spatial gene programs: %d genes, %s in-tissue bins", J, format(sum(keep), big.mark = ","))

  # ---- per-gene statistics ------------------------------------------------
  cnt <- binned$counts[keep, genes, drop = FALSE]
  gene_tot <- Matrix::colSums(cnt)
  gene_mean <- gene_tot / N
  gene_frac <- Matrix::colSums(cnt > 0) / N
  bin_tot <- Matrix::rowSums(binned$counts[keep, , drop = FALSE])
  dev_expl <- .rr_offset_dev(fit, binned, keep, genes)

  # ---- programs: factors + test components ----------------------------------
  strength <- .rr_or(fit$factor_strength, stats::setNames(sqrt(colSums(L^2)), fnames))
  prog_str <- .rr_or(fit$program_strength, stats::setNames(sqrt(colSums(sweep(L, 2, colMeans(L))^2)), fnames))
  uniform <- 1 - prog_str^2 / pmax(strength^2, 1e-12)
  ell <- .rr_or(fit$lengthscales, rep(NA_real_, K))
  tests <- .rr_parse_tests(test, fit, alpha, prog_str, strength)
  multi_pt <- sum(vapply(tests, function(t) t$kind == "program", TRUE)) > 1
  programs <- list()
  for (k in seq_len(K)) {
    programs[[fnames[k]]] <- list(id = fnames[k], name = fnames[k], label = paste0("M", k), pkg = paste("fit", fnames[k]), kind = "factor", loading = stats::setNames(L[, k], genes),
                                  strength = unname(strength[k]), program_strength = unname(prog_str[k]),
                                  uniform_share = unname(uniform[k]), lengthscale = unname(ell[k]))
  }
  comp_scores <- list()
  for (ti in seq_along(tests)) {
    t <- tests[[ti]]
    if (t$kind != "program") next
    V <- t$loadings; U <- t$scores
    for (k in seq_len(nrow(t$table))) {
      cn <- t$table$component[k]
      id <- if (multi_pt) paste0("T", ti, ".", cn) else cn
      dl <- paste0(if (multi_pt) paste0("T", ti, ":") else "", "D", sub("^PC", "", cn))
      ld <- if (!is.null(V) && cn %in% colnames(V)) V[genes, cn] else rep(NA_real_, J)
      programs[[id]] <- list(id = id, name = id, label = dl, pkg = paste("program test", cn), kind = "component", test = ti, component = cn,
                             loading = stats::setNames(ld, genes), share = t$table$share[k],
                             matched_factor = .rr_or(t$table$factor[k], NA),
                             matched_label = if (!is.na(.rr_or(t$table$factor[k], NA)) && t$table$factor[k] %in% fnames) paste0("M", match(t$table$factor[k], fnames)) else NA,
                             r_factor = .rr_or(t$table$r_factor[k], NA))
      if (!is.null(U) && cn %in% colnames(U) && nrow(U) == N) {
        full <- rep(NA_real_, g$nx * g$ny); full[keep] <- U[, cn]
        comp_scores[[id]] <- full
      }
      tests[[ti]]$table$id[k] <- id
    }
  }
  # ---- merged Programs table (detection axes matched to program maps) ----------
  raw <- if (is.data.frame(test)) list(test) else if (is.null(test)) list() else test
  pt_i <- which(vapply(raw, function(t) is.data.frame(t) && "component" %in% names(t), TRUE))[1]
  ft_i <- which(vapply(raw, function(t) is.data.frame(t) && !("component" %in% names(t)) && all(c("factor", "p") %in% names(t)), TRUE))[1]
  progtab <- .rff_match_programs(fit, if (!is.na(pt_i)) raw[[pt_i]] else NULL, if (!is.na(ft_i)) raw[[ft_i]] else NULL,
                                 alpha = alpha, min_r = min_r)
  axis_id <- function(cn) if (is.na(cn)) NA_character_ else if (multi_pt) paste0("T", pt_i, ".", cn) else cn
  PL <- attr(progtab, "loadings"); PF <- attr(progtab, "fields")
  alias <- list()
  for (i in which(!is.na(progtab$program))) {
    r <- progtab[i, ]; id <- r$program
    src <- if (!is.na(r$program_map)) r$program_map else axis_id(r$detection_axis)
    alias[[id]] <- src
    programs[[id]] <- list(id = id, name = id, label = id, kind = "program", status = r$status,
                           pkg = paste(c(if (!is.na(r$axis_label)) sprintf("%s (program test %s)", r$axis_label, r$detection_axis),
                                         if (!is.na(r$map_label)) sprintf("%s (fit %s)", r$map_label, r$program_map)), collapse = " + "),
                           loading = stats::setNames(PL[, r$key], genes), p = r$p, r = r$r,
                           axis = axis_id(r$detection_axis), map = r$program_map,
                           lengthscale = r$lengthscale, uniform_share = r$uniform_share, share = r$variance_share)
    comp_scores[[id]] <- PF[, r$key]
  }
  pids <- names(programs)
  # significance per program (any test)
  sig <- stats::setNames(rep(NA, length(pids)), pids)
  pmin_of <- stats::setNames(rep(NA_real_, length(pids)), pids)
  for (t in tests) {
    ids <- t$table$id; pv <- t$table$p_use
    for (i in seq_along(ids)) if (ids[i] %in% pids) {
      s <- !is.na(pv[i]) && pv[i] <= alpha
      sig[ids[i]] <- isTRUE(sig[ids[i]]) || s
      if (!is.na(pv[i])) pmin_of[ids[i]] <- min(pmin_of[ids[i]], pv[i], na.rm = TRUE)
    }
  }
  pmin_of[!is.finite(pmin_of)] <- NA
  for (id in names(alias)) {
    pr <- programs[[id]]
    sig[[id]] <- pr$status %in% c("Confirmed", "Candidate") && isTRUE(pr$p <= alpha)
    pmin_of[[id]] <- pr$p
  }

  # ---- known programs --------------------------------------------------------
  known <- NULL
  if (!is.null(known_programs)) known <- .rr_known(programs, known_programs, genes, top_genes)

  # ---- maps --------------------------------------------------------------------
  maps <- .rr_maps(fit, binned, keep, genes, programs, comp_scores, max_pixels, map_genes, n_map_genes, gene_tot)

  # ---- cells -------------------------------------------------------------------
  cellres <- NULL
  if (!is.null(cells)) cellres <- .rr_cells(fit, cells, cell_label_col, cell_xy, comp_scores[setdiff(names(comp_scores), names(alias))],
                                            max_cells, maps$f, alias)

  # ---- input summary, findings -------------------------------------------------
  off_desc <- .rr_offset_desc(fit, offset_info)
  mask_irreg <- .rr_mask_boundary(keep, g$nx, g$ny)
  findings <- .rr_findings(fit, tests, programs, sig, alpha, N, J, uniform, strength, st, gene_mean,
                           cellres, known, off_desc, mask_irreg, maps, offset_info, top_genes, list(), prog_str, progtab, min_r)

  # ---- reproduce code ------------------------------------------------------------
  repro <- .rr_reproduce(cl, fit, tests, off_desc, offset_info, out_dir, file, st, g, res_obj)

  # ---- export ----------------------------------------------------------------------
  files <- list()
  if (export) {
    files <- .rr_export(out_dir, fit, binned, keep, genes, programs, tests, comp_scores, cellres, known,
                        gene_tot, gene_mean, dev_expl, sig, pmin_of, repro, st, cl, alpha, g, progtab)
  }

  # ---- assemble ----------------------------------------------------------------------
  pkg_version <- tryCatch(as.character(utils::packageVersion("cohalu")), error = function(e) "unknown")
  hist_bins <- .rr_hist(bin_tot, 30)
  pt_json <- progtab; attr(pt_json, "loadings") <- NULL; attr(pt_json, "fields") <- NULL
  pt_json$id <- ifelse(!is.na(pt_json$program), pt_json$program, pt_json$program_map)
  prog_json <- lapply(programs, function(p) {
    p$loading <- I(unname(signif(p$loading, 4))); p$sig <- isTRUE(sig[[p$id]]); p$pmin <- pmin_of[[p$id]]; p
  })
  data <- list(
    title = title, date = format(Sys.time(), "%Y-%m-%d %H:%M"), version = pkg_version,
    alpha = alpha, top_genes = top_genes,
    genes = I(genes), gene_total = I(unname(gene_tot)), gene_mean = I(signif(unname(gene_mean), 4)),
    gene_frac = I(signif(unname(gene_frac), 3)),
    dev_expl = if (is.null(dev_expl)) NULL else I(signif(unname(dev_expl), 3)),
    programs = unname(prog_json), factor_ids = I(fnames), progtab = pt_json, min_r = min_r,
    has_program_test = !is.na(pt_i), has_factor_test = !is.na(ft_i),
    tests = lapply(tests, function(t) {
      t$loadings <- NULL; t$scores <- NULL
      t$table <- t$table[, setdiff(names(t$table), "top_genes"), drop = FALSE]
      t$nulls <- lapply(t$nulls, function(v) I(signif(v, 4))); t
    }),
    grid = list(nx = g$nx, ny = g$ny, bin_size = bs, xmin = g$xmin, ymin = g$ymin,
                f = maps$f, dnx = maps$dnx, dny = maps$dny),
    layers = maps$layers, gene_maps = maps$genes,
    hist = hist_bins, cells = if (is.null(cellres)) NULL else cellres$json, known = known$json
  )
  static <- list(
    input = .rr_html_input(fit, binned, keep, genes, gene_tot, bin_tot, dev_expl, off_desc, offset_info,
                           cells, cellres, known_programs, meta, st, tests, g, cl, mask_irreg),
    findings = .rr_html_findings(findings),
    outputs = .rr_html_outputs(files, out_dir, file, export),
    methods = .rr_html_methods(fit, tests),
    reproduce = paste0("<pre class=\"code\">", .rr_h(repro), "</pre>"),
    session = .rr_html_session(pkg_version, t_start)
  )
  overview <- .rr_html_overview(N, J, K, g, fit, tests, if (length(alias)) sig[names(alias)] else sig, alpha, maps)
  tpl <- .rr_template()
  html <- .rr_fill(tpl, list(
    TITLE = .rr_h(title),
    SUBTITLE = .rr_h(sprintf("%s  -  cohalu %s  -  %s in-tissue bins of %g, %d genes, %d program maps",
                             data$date, pkg_version, format(N, big.mark = ","), bs, J, K)),
    OVERVIEW = overview, INPUT = static$input, FINDINGS = static$findings, OUTPUTS = static$outputs,
    METHODS = static$methods, FVP = .rr_html_fvp(), REPRODUCE = static$reproduce, SESSION = static$session,
    DATA = .rr_json(data)
  ))
  writeLines(html, file, useBytes = TRUE)
  sz <- file.size(file)
  if (sz > 30e6) warning(sprintf("The report is large (%.0f MB); lower `max_pixels`, `n_map_genes` or `max_cells`.", sz / 1e6))
  if (isTRUE(open)) utils::browseURL(file)
  invisible(file)
}

# ---- small helpers -----------------------------------------------------------------

.rr_or <- function(a, b) if (is.null(a)) b else a

.rr_h <- function(s) {
  s <- as.character(s)
  s <- gsub("&", "&amp;", s, fixed = TRUE); s <- gsub("<", "&lt;", s, fixed = TRUE)
  s <- gsub(">", "&gt;", s, fixed = TRUE); gsub("\"", "&quot;", s, fixed = TRUE)
}

# base64 of a raw vector (RFC 4648), vectorised base R
.rr_b64 <- function(r) {
  n <- length(r)
  if (!n) return("")
  tbl <- c(LETTERS, letters, 0:9, "+", "/")
  pad <- (3L - n %% 3L) %% 3L
  v <- c(as.integer(r), integer(pad))
  m <- matrix(v, nrow = 3L)
  w <- m[1, ] * 65536L + m[2, ] * 256L + m[3, ]
  idx <- rbind(w %/% 262144L, (w %/% 4096L) %% 64L, (w %/% 64L) %% 64L, w %% 64L)
  ch <- tbl[as.vector(idx) + 1L]
  if (pad) ch[(length(ch) - pad + 1L):length(ch)] <- "="
  paste(ch, collapse = "")
}

# linear quantisation to 8 or 16 bits (code 0 = missing), base64 encoded
.rr_quant <- function(v, bits = 16L) {
  ok <- is.finite(v)
  if (!any(ok)) return(list(b64 = "", min = 0, max = 0, bits = bits))
  lo <- min(v[ok]); hi <- max(v[ok]); if (hi <= lo) hi <- lo + 1e-9
  top <- if (bits == 8L) 254 else 65534
  code <- integer(length(v))
  code[ok] <- as.integer(round((v[ok] - lo) / (hi - lo) * top)) + 1L
  r <- if (bits == 8L) as.raw(code) else writeBin(code, raw(), size = 2L, endian = "little")
  list(b64 = .rr_b64(r), min = lo, max = hi, bits = bits)
}

.rr_jstr <- function(s) {
  s <- enc2utf8(as.character(s)); na <- is.na(s)
  s <- gsub("\\", "\\\\", s, fixed = TRUE); s <- gsub("\"", "\\\"", s, fixed = TRUE)
  s <- gsub("\n", "\\n", s, fixed = TRUE); s <- gsub("\r", "\\r", s, fixed = TRUE)
  s <- gsub("\t", "\\t", s, fixed = TRUE); s <- gsub("[\001-\037]", " ", s)
  s <- gsub("</", "<\\/", s, fixed = TRUE)
  out <- paste0("\"", s, "\""); out[na] <- "null"; out
}

# Minimal JSON writer: length-1 atomics are scalars unless wrapped in I();
# named lists / data frames are objects (data frames by column).
.rr_json <- function(x, digits = 6) {
  if (is.null(x)) return("null")
  if (is.data.frame(x)) x <- lapply(as.list(x), I)
  if (is.list(x)) {
    if (!length(x)) return(if (is.null(names(x))) "[]" else "{}")
    vals <- vapply(x, .rr_json, "", digits = digits)
    if (!is.null(names(x)) && all(nzchar(names(x))))
      return(paste0("{", paste0(.rr_jstr(names(x)), ":", vals, collapse = ","), "}"))
    return(paste0("[", paste(vals, collapse = ","), "]"))
  }
  arr <- inherits(x, "AsIs") || length(x) != 1
  if (is.factor(x)) x <- as.character(x)
  x <- unclass(x); attributes(x) <- NULL
  el <- if (is.character(x)) .rr_jstr(x) else if (is.logical(x)) ifelse(is.na(x), "null", ifelse(x, "true", "false")) else
    if (is.numeric(x)) { o <- rep("null", length(x)); f <- is.finite(x)
      o[f] <- if (is.integer(x)) as.character(x[f]) else formatC(x[f], digits = digits, format = "g"); gsub(" ", "", o) } else
      .rr_jstr(as.character(x))
  if (arr) paste0("[", paste(el, collapse = ","), "]") else el
}

.rr_fill <- function(tpl, vals) {
  for (nm in names(vals)) {
    parts <- strsplit(tpl, paste0("{{", nm, "}}"), fixed = TRUE)[[1]]
    if (length(parts) > 1) tpl <- paste(parts, collapse = vals[[nm]])
  }
  tpl
}

.rr_template <- function() {
  f <- system.file("rff_report", "template.html", package = "cohalu")
  if (!nzchar(f)) stop("Report template not found (inst/rff_report/template.html).")
  paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

.rr_hist <- function(v, n = 30) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NULL)
  br <- pretty(range(v), n)
  cnt <- tabulate(findInterval(v, br, rightmost.closed = TRUE, all.inside = TRUE), length(br) - 1L)
  list(breaks = I(br), counts = I(as.integer(cnt)),
       q = as.list(stats::setNames(stats::quantile(v, c(0, .05, .25, .5, .75, .95, 1), names = FALSE),
                                   c("min", "q05", "q25", "median", "q75", "q95", "max"))),
       mean = mean(v), zero = mean(v == 0))
}

.rr_fmt <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(signif(x, d), format = "fg", digits = d, flag = "#"))
.rr_num <- function(x, d = 3) {
  o <- ifelse(is.na(x), "NA", format(signif(x, d), big.mark = ",", scientific = FALSE, trim = TRUE, drop0trailing = TRUE))
  o
}
.rr_p <- function(p) ifelse(is.na(p), "-", ifelse(p < 0.001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p)))

.rr_top <- function(v, n, sign = 0) {
  v <- v[is.finite(v)]
  if (sign > 0) v <- v[v > 0] else if (sign < 0) v <- v[v < 0]
  o <- order(-abs(v)); names(v)[o[seq_len(min(n, length(o)))]]
}

# ---- tests ---------------------------------------------------------------------------

.rr_parse_tests <- function(test, fit, alpha, prog_str, strength) {
  if (is.null(test)) return(list())
  if (is.data.frame(test)) test <- list(test)
  if (!is.list(test)) stop("`test` must be a data frame from rff_factor_test() / rff_program_test() or a list of them.")
  nm <- names(test); if (is.null(nm)) nm <- rep("", length(test))
  out <- list()
  for (i in seq_along(test)) {
    t <- test[[i]]
    if (!is.data.frame(t)) stop("`test` must be a data frame from rff_factor_test() / rff_program_test() or a list of them.")
    if ("component" %in% names(t) && "share" %in% names(t)) {
      kind <- "program"
      null <- .rr_or(attr(t, "null", exact = TRUE), "parametric")
      ns <- attr(t, "null_shares", exact = TRUE)
      tab <- data.frame(id = as.character(t$component), component = as.character(t$component),
                        share = t$share, p = t$p, p_use = t$p,
                        factor = if ("factor" %in% names(t)) as.character(t$factor) else NA,
                        r_factor = if ("r_factor" %in% names(t)) t$r_factor else NA,
                        top_genes = if ("top_genes" %in% names(t)) as.character(t$top_genes) else NA,
                        stringsAsFactors = FALSE)
      sequential <- .rr_or(attr(t, "sequential", exact = TRUE), anyNA(t$p) || (length(ns) > 1 && !identical(ns[[1]], ns[[2]])))
      nulls <- lapply(seq_len(nrow(tab)), function(k) {
        v <- if (length(ns) >= k && !is.null(ns[[k]])) ns[[k]] else if (length(ns)) ns[[1]] else numeric(0)
        as.numeric(v)
      })
      obs <- tab$share
      n_boot <- if (length(ns)) length(ns[[1]]) else NA_integer_
      stat_name <- "variance share of the component"
    } else if (all(c("factor", "p") %in% names(t))) {
      kind <- "factor"; null <- "parametric"
      statistic <- .rr_or(attr(t, "statistic", exact = TRUE), "program")
      tab <- data.frame(id = as.character(t$factor), factor = as.character(t$factor), p = t$p,
                        p_sequential = if ("p_sequential" %in% names(t)) t$p_sequential else NA_real_,
                        stringsAsFactors = FALSE)
      sequential <- "p_sequential" %in% names(t)
      tab$p_use <- if (sequential) tab$p_sequential else tab$p
      nm_ok <- tab$factor %in% names(prog_str)
      obs <- ifelse(nm_ok, if (statistic == "program") prog_str[tab$factor] else strength[tab$factor], NA)
      if ("program_strength" %in% names(t) && statistic == "program") obs <- t$program_strength
      nm0 <- attr(t, "null_max", exact = TRUE); nseq <- attr(t, "null_sequential", exact = TRUE)
      rk <- rank(-obs, ties.method = "first")
      nulls <- lapply(seq_len(nrow(tab)), function(k) {
        v <- if (sequential && !is.null(nseq) && length(nseq) >= rk[k] && !is.null(nseq[[rk[k]]])) nseq[[rk[k]]] else nm0
        as.numeric(.rr_or(v, numeric(0)))
      })
      n_boot <- if (!is.null(nm0)) length(nm0) else NA_integer_
      stat_name <- if (statistic == "program") "program strength (loading norm after removing the all-gene mean)" else "factor strength (loading norm)"
      tab$top_genes <- if ("top_genes" %in% names(t)) as.character(t$top_genes) else NA
    } else stop("Test ", i, " is not an output of rff_factor_test() or rff_program_test().")
    label <- if (nzchar(nm[i])) nm[i] else if (kind == "factor") "rff_factor_test" else sprintf("rff_program_test (null = \"%s\")", null)
    tab$sig <- !is.na(tab$p_use) & tab$p_use <= alpha
    out[[i]] <- list(label = label, kind = kind, null_type = null, sequential = isTRUE(sequential),
                     n_boot = n_boot, highpass = attr(t, "highpass", exact = TRUE), bandwidth = attr(t, "bandwidth", exact = TRUE),
                     max_iter = attr(t, "max_iter", exact = TRUE), elapsed = attr(t, "elapsed", exact = TRUE), stat_name = stat_name,
                     table = tab, obs = I(as.numeric(obs)), nulls = nulls,
                     scores = attr(t, "scores", exact = TRUE), loadings = attr(t, "loadings", exact = TRUE),
                     caveat = .rr_caveat(kind, null, attr(t, "highpass", exact = TRUE)))
  }
  out
}

# Pick one window of a spatial_rff_joint fit for rff_report(); convert a joint
# test (rff_program_test_joint(), discovery mode) to that window's view.
.rr_joint_window <- function(jf, binned, test, window) {
  wn <- names(jf$fits); if (is.null(wn)) wn <- paste0("window", seq_along(jf$fits))
  if (is.null(window)) window <- wn[which.max(rowSums(jf$amplitude))]
  if (is.numeric(window)) window <- wn[window]
  if (length(window) != 1 || !window %in% wn) stop("`window` must name (or index) one window of the joint fit.")
  w <- match(window, wn)
  if (!inherits(binned, "binned_transcripts")) {
    if (!is.list(binned)) stop("`binned` must be a binned_transcripts object or a list of them (one per window).")
    binned <- if (!is.null(names(binned)) && window %in% names(binned)) binned[[window]] else binned[[w]]
  }
  f <- jf$fits[[w]]
  if (is.data.frame(test) && identical(attr(test, "null", exact = TRUE), "shift_joint")) {
    sc <- attr(test, "scores", exact = TRUE); ws <- attr(test, "window_share", exact = TRUE)
    U <- if (!is.null(names(sc)) && window %in% names(sc)) sc[[window]] else NULL
    Fm <- f$field_grid[f$in_tissue, colnames(f$L), drop = FALSE]
    cc <- if (!is.null(U)) abs(suppressWarnings(stats::cor(U, Fm))) else matrix(NA_real_, nrow(test), ncol(Fm))
    cc[!is.finite(cc)] <- 0
    t2 <- data.frame(component = test$component,
                     share = if (!is.null(ws) && window %in% rownames(ws)) unname(ws[window, test$component]) else test$share,
                     p = test$p, factor = colnames(Fm)[apply(cc, 1, which.max)], r_factor = apply(cc, 1, max),
                     top_genes = test$top_genes, stringsAsFactors = FALSE)
    attr(t2, "scores") <- U; attr(t2, "loadings") <- attr(test, "loadings", exact = TRUE)
    attr(t2, "null_shares") <- attr(test, "null_shares", exact = TRUE); attr(t2, "null") <- "shift_joint"
    attr(t2, "sequential") <- TRUE
    test <- t2
  }
  amp <- jf$amplitude
  meta <- list(joint_fit = sprintf("shared loadings over %d windows (fit_spatial_rff_joint); this report shows window %s. Variance shares in the test table are this window's; p-values are the joint test's.",
                                   length(wn), window),
               amplitude_this_window = paste(sprintf("%s %.2f", colnames(amp), amp[w, ]), collapse = "; "),
               amplitude_median_all_windows = paste(sprintf("%s %.2f", colnames(amp), apply(amp, 2, stats::median)), collapse = "; "))
  list(fit = f, binned = binned, test = test, meta = meta)
}

.rr_caveat <- function(kind, null, highpass = NULL) {
  if (identical(null, "shift_joint")) {
    return(paste(
      "Joint gene-shift test (rff_program_test_joint): processed residuals of all windows are pooled into one gene x gene",
      "covariance (shared loadings); surrogates shift every gene independently within every window. Significance means a",
      "program shared by several genes and recurring across windows. p-values refer to the pooled data, not to this window;",
      "shown shares are this window's. Validate programs in held-out data with the confirmatory mode (`loadings =`)."))
  }
  if (kind == "factor") {
    return(paste(
      "Parametric bootstrap: counts are simulated from the fitted null model (bin area, offset, gene intercepts,",
      "cellularity field and dispersions, no factors), the model is refitted to each simulated data set and every",
      "observed factor is compared with the largest statistic of the refits (family-wise over factors; with",
      "sequential testing the k-th factor is compared with refits to data that contain the k - 1 stronger factors).",
      "Significance means: stronger than any factor fitted to data without factor structure. The null assumes that",
      "what the offset and the cellularity field do not explain is independent noise; gene-specific spatially",
      "autocorrelated structure (usual in real tissue) makes it anti-conservative. Refits must use the fit's max_iter."))
  }
  if (identical(null, "shift")) {
    return(paste(
      "Gene-shift surrogates: every gene's counts are moved together with its null mean by an independent random",
      "rigid transformation of the grid, which keeps each gene's own spatial autocorrelation and destroys only the",
      "alignment between genes; component k is compared with the k-th variance share of the surrogates (parallel",
      "analysis). Significance means: a spatial program shared by several genes beyond what each gene's own",
      "structure produces. Gene-specific fields at scales larger than the program can mask it",
      if (is.null(highpass)) "(consider `highpass`)." else sprintf("(here `highpass` = %g was used).", highpass),
      "On real tissue with an irregular mask or a single dominant structure, borderline p-values should be read with care."))
  }
  paste(
    "Parametric null: counts are simulated from the fitted null model (offset, intercepts, cellularity field and",
    "the all-gene part of the factors), independent noise around a smooth mean. When genes have their own spatially",
    "autocorrelated residual structure (usual in real tissue) almost every tissue is called: significance then means",
    "only 'structure the offset misses', not necessarily a shared gene program. Prefer null = \"shift\" on real data.")
}

# ---- offset --------------------------------------------------------------------------

.rr_offset_dev <- function(fit, binned, keep, genes) {
  O <- fit$offset_matrix
  if (is.null(O) || !all(genes %in% colnames(O))) return(NULL)
  res <- stats::setNames(rep(NA_real_, length(genes)), genes)
  dev <- function(y, mu) 2 * colSums(ifelse(y > 0, y * log(y / mu), 0) - (y - mu))
  for (blk in split(seq_along(genes), ceiling(seq_along(genes) / 50))) {
    Y <- as.matrix(binned$counts[keep, genes[blk], drop = FALSE])
    E <- exp(O[keep, genes[blk], drop = FALSE])
    E <- sweep(E, 2, colSums(Y) / pmax(colSums(E), 1e-300), "*")
    m0 <- matrix(colMeans(Y), nrow(Y), ncol(Y), byrow = TRUE)
    d0 <- dev(Y, pmax(m0, 1e-300)); d1 <- dev(Y, pmax(E, 1e-300))
    res[blk] <- ifelse(d0 > 0, 1 - d1 / d0, NA)
  }
  res
}

.rr_offset_desc <- function(fit, offset_info) {
  off <- .rr_or(fit$offset, if (!is.null(fit$offset_matrix)) "matrix" else "area")
  oi <- .rr_or(offset_info, list())
  om <- fit$offset_matrix
  method <- .rr_or(oi$method, .rr_or(attr(om, "method", exact = TRUE), NA))
  k <- .rr_or(oi$k, attr(om, "k", exact = TRUE))
  txt <- switch(off,
    area = "No known structure: only the bin area is used as offset; overall cellularity is modelled by the latent field f0.",
    smoothed_total = sprintf("Smoothed total transcript density (Gaussian bandwidth %s) as offset; no cellularity field.",
                             .rr_or(fit$settings$offset_bandwidth, "?")),
    matrix = paste0("Per-bin, per-gene log offset (known structure, 'residual RFLVM')",
                    if (!is.na(method)) sprintf(", made with rff_offset(method = \"%s\"%s)", method,
                                                if (!is.null(k) && method %in% c("kmeans", "pca")) sprintf(", k = %s", k) else "") else
                      ". How it was made is not recorded in the fit (pass `offset_info = list(method = ...)`)",
                    "."),
    paste("Offset:", off))
  list(type = off, method = method, k = k, text = txt, data_driven = isTRUE(method %in% c("kmeans", "pca")))
}

# fraction of in-tissue bins on the mask boundary (4-neighbourhood)
.rr_mask_boundary <- function(keep, nx, ny) {
  m <- matrix(keep, nx, ny)
  if (!any(m)) return(list(boundary = NA, fill = 0))
  pad <- matrix(FALSE, nx + 2, ny + 2); pad[2:(nx + 1), 2:(ny + 1)] <- m
  inner <- pad[1:nx, 2:(ny + 1)] & pad[3:(nx + 2), 2:(ny + 1)] & pad[2:(nx + 1), 1:ny] & pad[2:(nx + 1), 3:(ny + 2)]
  bnd <- m & !inner
  # boundary excluding the grid edge (a rectangular window is regular)
  ib <- bnd; ib[c(1, nx), ] <- FALSE; ib[, c(1, ny)] <- FALSE
  list(boundary = sum(ib) / sum(m), fill = mean(m))
}

# ---- known programs ------------------------------------------------------------------

.rr_known <- function(programs, known_programs, genes, top_n) {
  kn <- names(known_programs)
  rows <- list(); json_rows <- list()
  for (p in programs) {
    v <- p$loading
    if (all(!is.finite(v))) next
    n_top <- min(top_n, max(5L, ceiling(0.1 * length(genes))))
    top <- .rr_top(v, n_top)
    cells <- list()
    best <- NA; best_r <- 0; best_ph <- 1
    for (k in kn) {
      set <- intersect(known_programs[[k]], genes)
      ind <- as.numeric(genes %in% set)
      r <- if (length(set) >= 1 && length(set) < length(genes)) suppressWarnings(stats::cor(v, ind)) else NA
      ov <- length(intersect(top, set))
      ml <- if (length(set)) mean(v[set]) else NA
      ph <- if (ov > 0) stats::phyper(ov - 1, length(set), length(genes) - length(set), length(top), lower.tail = FALSE) else 1
      cells[[k]] <- list(r = r, overlap = ov, p_overlap = ph, n_set = length(set), n_given = length(known_programs[[k]]), mean_loading = ml,
                         genes = I(intersect(top, set)))
      rows[[length(rows) + 1]] <- data.frame(program = p$id, known_program = k, n_genes_in_panel = length(set),
                                             n_genes_given = length(known_programs[[k]]), cor_loading_membership = r,
                                             top_overlap = ov, top_n = length(top), p_overlap = ph, mean_loading = ml,
                                             overlapping_genes = paste(intersect(top, set), collapse = ";"))
      if (is.finite(r) && abs(r) > abs(best_r)) { best_r <- r; best <- k; best_ph <- ph }
    }
    status <- if (!is.na(best) && abs(best_r) >= 0.3 && best_ph < 0.05) "known" else
      if (!is.na(best) && (abs(best_r) >= 0.3 || best_ph < 0.05)) "partial" else "novel"
    json_rows[[p$id]] <- list(id = p$id, cells = cells, best = best, best_r = best_r, status = status)
  }
  tab <- if (length(rows)) do.call(rbind, rows) else NULL
  list(table = tab, json = list(sets = I(kn), rows = unname(json_rows)), status = vapply(json_rows, `[[`, "", "status"),
       best = lapply(json_rows, function(x) x[c("best", "best_r")]))
}

# ---- maps ------------------------------------------------------------------------------

.rr_maps <- function(fit, binned, keep, genes, programs, comp_scores, max_pixels, map_genes, n_map_genes, gene_tot) {
  g <- fit$grid; nx <- g$nx; ny <- g$ny; n <- nx * ny
  f <- max(1L, as.integer(ceiling(sqrt(n / max_pixels))))
  if (f > 1) warning(sprintf("The %d x %d grid is shown block-averaged by %d x %d bins in the report (max_pixels = %d).",
                             nx, ny, f, f, as.integer(max_pixels)))
  dnx <- as.integer(ceiling(nx / f)); dny <- as.integer(ceiling(ny / f))
  ix <- (seq_len(n) - 1L) %% nx; iy <- (seq_len(n) - 1L) %/% nx
  blk <- (ix %/% f) + dnx * (iy %/% f) + 1L
  kb <- blk[keep]; nb <- tabulate(kb, dnx * dny)
  agg <- function(v, fun = c("mean", "sum")) {
    fun <- match.arg(fun)
    vv <- v[keep]; ok <- is.finite(vv)
    s <- numeric(dnx * dny)
    if (any(ok)) { r <- rowsum(vv[ok], kb[ok]); s[as.integer(rownames(r))] <- r[, 1] }
    cntok <- tabulate(kb[ok], dnx * dny)
    out <- if (fun == "mean") s / pmax(cntok, 1) else s
    out[nb == 0 | cntok == 0] <- NA
    out
  }
  layers <- list()
  add <- function(id, name, group, scale, v, unit = "", desc = "") {
    layers[[length(layers) + 1]] <<- list(id = id, name = name, group = group, scale = scale, unit = unit, desc = desc,
                                          q = .rr_quant(agg(v), 16L))
  }
  fg <- fit$field_grid
  for (p in programs) {
    if (p$kind == "program" && !is.null(comp_scores[[p$id]]))
      add(p$id, p$label, "program", "div", comp_scores[[p$id]], if (!is.na(.rr_or(p$map, NA))) "field (sd units)" else "score",
          sprintf("Program %s (%s): %s", p$label, p$pkg, if (!is.na(.rr_or(p$map, NA))) "smooth field of its program map." else "score map of its detection axis (no matching program map)."))
    if (p$kind == "factor" && p$id %in% colnames(fg))
      add(p$id, p$label, "factor", "div", fg[, p$id], "field (sd units)",
          sprintf("Program map %s (%s): smooth field with unit variance over the tissue; log-intensity effect on gene j is L[j] x field.", p$label, p$pkg))
    if (p$kind == "component" && !is.null(comp_scores[[p$id]]))
      add(p$id, p$label, "component", "div", comp_scores[[p$id]], "score",
          sprintf("Detection axis %s (%s): score of the smoothed-residual PCA component used by rff_program_test().", p$label, p$pkg))
  }
  if (isTRUE(fit$has_density) && "density" %in% colnames(fg))
    add("density", "cellularity f0", "context", "div", fg[, "density"], "field (sd units)",
        sprintf("Cellularity field f0 shared by all genes (log scale x sigma0 = %.2f).", .rr_or(fit$sigma0, NA)))
  tot <- Matrix::rowSums(binned$counts)
  layers[[length(layers) + 1]] <- list(id = "total", name = "total counts", group = "context", scale = "seq",
                                       unit = "transcripts / display pixel", desc = "Observed transcripts (all genes) per bin, summed per display pixel.",
                                       q = .rr_quant(log1p(agg(tot, "sum")), 16L), tf = "expm1")
  if (!is.null(fit$offset_matrix)) {
    om <- rowMeans(fit$offset_matrix); om[!keep] <- NA
    add("offset", "offset (mean)", "context", "seq", om - mean(om[keep]), "log",
        "Mean over genes of the known-structure log offset (centred).")
  } else if (identical(fit$offset, "smoothed_total") && !is.null(fit$offset_grid)) {
    add("offset", "offset (smoothed total)", "context", "seq", fit$offset_grid, "log", "Log smoothed total density used as offset.")
  }
  layers <- layers[order(match(vapply(layers, `[[`, "", "group"), c("program", "factor", "component", "context")))]
  # gene maps: observed counts and null expectation (offset + intercept + cellularity)
  sel <- unique(c(intersect(map_genes, genes),
                  unlist(lapply(programs, function(p) .rr_top(p$loading, 4))),
                  names(sort(gene_tot, decreasing = TRUE))))
  sel <- sel[seq_len(min(length(sel), n_map_genes))]
  la <- 2 * log(g$bin_size)
  dens <- if (isTRUE(fit$has_density)) fit$sigma0 * fg[, "density"] else 0
  gm <- list()
  for (gn in sel) {
    y <- numeric(n); y[keep] <- as.numeric(binned$counts[keep, gn])
    eta <- la + (if (!is.null(fit$offset_matrix)) fit$offset_matrix[, gn] else if (identical(fit$offset, "smoothed_total")) fit$offset_grid - la else 0) +
      fit$alpha[[gn]] + dens
    gm[[length(gm) + 1]] <- list(gene = gn, obs = .rr_quant(log1p(agg(y, "sum")), 8L),
                                 exp = .rr_quant(log1p(agg(exp(eta), "sum")), 8L))
  }
  list(layers = layers, genes = gm, f = f, dnx = dnx, dny = dny, n_map_genes = length(sel))
}

# ---- cells -------------------------------------------------------------------------------

.rr_cells <- function(fit, cells, label_col, xy, comp_scores, max_cells, f, alias = list()) {
  cells <- as.data.frame(cells)
  miss <- setdiff(c(xy, label_col), names(cells))
  if (length(miss)) stop("Columns not found in `cells`: ", paste(miss, collapse = ", "))
  n_out <- 0L
  fv <- withCallingHandlers(rff_fields(fit, cells, x_col = xy[1], y_col = xy[2]), warning = function(w) {
    m <- conditionMessage(w)
    if (grepl("outside the fitted tissue", m)) { n_out <<- as.integer(sub(" of.*", "", m)); invokeRestart("muffleWarning") }
  })
  g <- fit$grid
  x <- as.numeric(cells[[xy[1]]]); y <- as.numeric(cells[[xy[2]]])
  ix <- floor((x - g$xmin) / g$bin_size); iy <- floor((y - g$ymin) / g$bin_size)
  inside <- ix >= 0 & ix < g$nx & iy >= 0 & iy < g$ny
  bidx <- ifelse(inside, ix + 1 + g$nx * iy, NA)
  for (id in names(comp_scores)) fv[[id]] <- comp_scores[[id]][bidx]
  for (id in names(alias)) if (!is.null(fv[[alias[[id]]]])) fv[[id]] <- fv[[alias[[id]]]]
  lab <- if (is.null(label_col)) rep("all cells", nrow(cells)) else as.character(cells[[label_col]])
  lab[is.na(lab)] <- "NA"
  tabl <- sort(table(lab), decreasing = TRUE)
  levs <- names(tabl)
  fields <- names(fv)
  summ <- lapply(fields, function(fl) {
    v <- fv[[fl]]
    lapply(levs, function(l) {
      vv <- v[lab == l]; vv <- vv[is.finite(vv)]
      if (!length(vv)) return(list(n = 0L))
      q <- stats::quantile(vv, c(.05, .25, .5, .75, .95), names = FALSE)
      list(n = length(vv), mean = mean(vv), q05 = q[1], q25 = q[2], q50 = q[3], q75 = q[4], q95 = q[5])
    })
  })
  names(summ) <- fields
  nc <- nrow(cells)
  pick <- if (nc > max_cells) unique(round(seq(1, nc, length.out = max_cells))) else seq_len(nc)
  pts <- list(x = .rr_quant(x[pick], 16L), y = .rr_quant(y[pick], 16L),
              label = I(match(lab[pick], levs) - 1L),
              values = lapply(stats::setNames(fields, fields), function(fl) .rr_quant(fv[[fl]][pick], 8L)))
  out <- data.frame(cells[, c(xy, label_col), drop = FALSE], fv, check.names = FALSE)
  list(json = list(labels = I(levs), counts = I(as.integer(tabl)), fields = I(fields), summary = summ,
                   points = pts, n = nc, n_drawn = length(pick), n_outside = n_out, label_col = label_col),
       table = out, labels = levs, counts = tabl, summ = summ, n_outside = n_out, fields = fields)
}

# ---- findings ------------------------------------------------------------------------------

.rr_findings <- function(fit, tests, programs, sig, alpha, N, J, uniform, strength, st, gene_mean,
                         cellres, known, off_desc, mask, maps, offset_info, top_genes, pc_r = list(), prog_str = NULL,
                         progtab = NULL, min_r = 0.5) {
  F <- list()
  add <- function(level, text, head = NULL, pin = FALSE) F[[length(F) + 1]] <<- list(level = level, head = head, text = text, pin = pin)
  fnames <- names(strength)
  lb <- function(id) vapply(id, function(i) .rr_or(programs[[i]]$label, i), "")
  lbp <- function(id) vapply(id, function(i) if (is.null(programs[[i]])) i else sprintf("%s (%s)", programs[[i]]$label, programs[[i]]$pkg), "")
  # significance
  if (!length(tests)) {
    add("warning", "No significance test was supplied, so no program map can be called significant. Program strengths alone are not evidence: fits to data without any program also produce non-zero program maps. Run rff_program_test(fit, binned, null = \"shift\") (recommended on real data) or rff_factor_test().",
        "No test")
  }
  for (t in tests) {
    ns <- sum(t$table$sig); nt <- sum(!is.na(t$table$p_use))
    what <- if (t$kind == "factor") "program maps (fit factors)" else "detection axes (program test PCs)"
    lvl <- "info"
    add(lvl, sprintf("%s: %d of %d %s significant at alpha = %g (%s null, %s, %s surrogates / bootstrap data sets%s).",
                     t$label, ns, nrow(t$table), what, alpha, t$null_type,
                     if (t$sequential) "sequential p-values" else "max-statistic p-values",
                     ifelse(is.na(t$n_boot), "?", t$n_boot),
                     if (!is.null(t$highpass)) sprintf(", highpass %g", t$highpass) else ""),
        sprintf("Details: %d significant", ns))
    if (is.finite(t$n_boot) && 1 / (t$n_boot + 1) > alpha)
      add("warning", sprintf("%s used only %d null data sets: the smallest attainable p-value (%.3f) is above alpha = %g, so nothing can be significant.",
                             t$label, t$n_boot, 1 / (t$n_boot + 1), alpha), "Too few null draws")
    if (t$kind == "program" && identical(t$null_type, "parametric"))
      add("caution", "rff_program_test() with the parametric null: significance means 'structure the offset misses' (independent-noise null); with gene-specific spatially autocorrelated residuals it calls detection axes in almost every tissue. Confirm with null = \"shift\".",
          "Calibration")
    if (t$kind == "program" && identical(t$null_type, "shift"))
      add("info", "Gene-shift null: a significant detection axis is a spatial pattern shared by several genes beyond each gene's own autocorrelation. It can still be structure that the offset should have contained (a missing cell type or domain): judge it by its genes and map.",
          "Reading the gene-shift null")
    if (t$kind == "factor")
      add("info", "rff_factor_test() compares program maps (fit factors) with refits to data simulated without factors (independent noise around the null mean); in real tissue, gene-specific autocorrelated structure makes it anti-conservative. In simulations it was also less powerful than rff_program_test(null = \"shift\").",
          "Calibration")
    if (!is.null(t$max_iter) && !is.null(st$max_iter) && t$max_iter != st$max_iter)
      add("warning", sprintf("%s refits used max_iter = %s but the fit used %s: null statistics are not comparable.", t$label, t$max_iter, st$max_iter),
          "max_iter mismatch")
  }
  # headline: the Programs table
  main_ids <- names(programs)[vapply(programs, function(p) p$kind == "program", TRUE)]
  if (!is.null(progtab)) {
    nst <- table(factor(progtab$status, c("Confirmed", "Candidate", "Exploratory", "Not supported", "Cellularity/technical")))
    kinds0 <- vapply(tests, `[[`, "", "kind")
    how <- if ("program" %in% kinds0) {
      tp <- tests[[which(kinds0 == "program")[1]]]
      sprintf("detected by rff_program_test (%s null, %s draws) at alpha = %g", tp$null_type, tp$n_boot, alpha)
    } else if ("factor" %in% kinds0) sprintf("detected by rff_factor_test at alpha = %g (no program test: the gene-shift program test is recommended)", alpha) else "no test"
    if (!length(tests)) {
      add("warning", sprintf("Exploratory maps, no significance: %d program map(s) shown as exploratory programs. Run rff_programs() or rff_program_test() to decide which programs exist.",
                             nst[["Exploratory"]]), "Exploratory only", pin = TRUE)
    } else {
      add(if (nst[["Confirmed"]] + nst[["Candidate"]] > 0) "ok" else "info",
          sprintf("%d confirmed and %d candidate program(s), %s.%s%s", nst[["Confirmed"]], nst[["Candidate"]], how,
                  if (nst[["Not supported"]]) sprintf(" %d strong program map(s) not supported by the test.", nst[["Not supported"]]) else "",
                  if (nst[["Cellularity/technical"]]) sprintf(" %d cellularity-like map(s) set aside.", nst[["Cellularity/technical"]]) else ""),
          sprintf("%d program%s", nst[["Confirmed"]] + nst[["Candidate"]], if (nst[["Confirmed"]] + nst[["Candidate"]] == 1) "" else "s"), pin = TRUE)
    }
  }
  for (id in main_ids) {
    p <- programs[[id]]
    tg <- .rr_top(p$loading, 6)
    txt <- sprintf("Strongest genes (sign of loading): %s; the sign of a field is arbitrary, genes with the same sign go up together where the field is high.",
                   paste(sprintf("%s (%s)", tg, ifelse(p$loading[tg] > 0, "+", "-")), collapse = ", "))
    txt <- paste(sprintf("%s%s.", if (is.finite(.rr_or(p$p, NA))) sprintf("p = %s; ", .rr_p(p$p)) else "", p$pkg), txt)
    if (!is.na(.rr_or(p$map, NA)))
      txt <- paste(txt, sprintf("Length scale %.3g%s.", p$lengthscale, if (is.finite(.rr_or(p$r, NA))) sprintf("; detection axis and program map agree with r = %.2f", p$r) else ""))
    if (identical(p$status, "Candidate"))
      txt <- paste(txt, sprintf("No program map matches it with |r| >= %.1f: the RFLVM does not represent this structure well; genes and map are the detection axis's own.", min_r))
    if (!is.null(cellres)) {
      sm <- cellres$summ[[id]]
      if (!is.null(sm) && length(cellres$labels) > 1) {
        mv <- vapply(sm, function(s) .rr_or(s$mean, NA_real_), 0)
        ok <- is.finite(mv) & vapply(sm, function(s) s$n, 0L) >= 10
        if (any(ok)) {
          o <- order(-mv); o <- o[ok[o]]
          txt <- paste(txt, sprintf("Field highest in cells labelled %s; lowest in %s.",
                                    paste(cellres$labels[utils::head(o, 2)], collapse = ", "),
                                    cellres$labels[utils::tail(o, 1)]))
        }
      }
    }
    if (!is.null(known) && !is.null(known$status[[id]])) {
      b <- known$best[[id]]
      txt <- paste(txt, switch(known$status[[id]],
                               known = sprintf("Matches known program '%s' (r = %.2f): probably not new.", b$best, b$best_r),
                               partial = sprintf("Partial overlap with '%s' (r = %.2f).", b$best, b$best_r),
                               novel = "No clear overlap with the known programs given (candidate novel program)."))
    }
    top2 <- sort(p$loading^2, decreasing = TRUE)
    if (length(top2) >= 3 && sum(top2[1:2]) / sum(top2) > 0.7)
      txt <- paste(txt, "Dominated by one or two genes: check them for technical artefacts (probe, segmentation, low counts).")
    add(switch(p$status, Confirmed = "ok", Candidate = "caution", "info"), txt, sprintf("%s: %s", id, p$status))
  }
  # details: detection axes and program maps
  kinds <- vapply(tests, `[[`, "", "kind")
  has_ft <- "factor" %in% kinds; has_pt <- "program" %in% kinds
  if (!is.null(progtab) && any(progtab$status == "Not supported")) {
    ns <- progtab[progtab$status == "Not supported", ]
    add("caution", sprintf("Program map(s) %s look strong (%s) but no significant detection axis matches them (|r| >= %.1f): not supported by the program test; listed under 'Not supported'.",
                           paste(sprintf("%s (fit %s)", ns$map_label, ns$program_map), collapse = ", "),
                           if (has_ft) "significant in rff_factor_test" else "large program strength", min_r), "Not supported")
  }
  if (has_ft && !has_pt)
    add("info", "Only rff_factor_test() was run: it answers 'is this program map (fit factor) stronger than factors fitted to data without factor structure (parametric null)?'. It does not test for a program shared by several genes beyond each gene's own spatial structure; run rff_program_test(null = \"shift\") for that.",
        "Question answered")
  if (has_pt && !has_ft)
    add("info", "Only rff_program_test() was run: it answers 'how many spatial gene programs does the smoothed residual contain beyond the null?' (detection axes). Program maps are not tested individually; read them through the detection axes they match.",
        "Question answered")
  # cellularity-like factors
  cl <- if (!is.null(progtab)) progtab$program_map[progtab$status == "Cellularity/technical"] else names(uniform)[uniform > 0.5 & strength > 0.1 * max(strength)]
  if (length(cl))
    add("caution", sprintf("Cellularity-like program map(s) %s: most of the loading is shared by all genes (uniform share > 0.5), i.e. extra density, not a gene program; listed under 'Cellularity / technical'.",
                           paste(lbp(cl), collapse = ", ")), "Cellularity")
  # ARD
  if (isTRUE(st$ard > 0)) {
    shr <- names(strength)[strength < 0.1 * max(strength)]
    if (length(shr)) add("info", sprintf("Program map(s) %s were shrunk to (near) zero by the ARD penalty (ard = %g): effectively unused.",
                                         paste(lbp(shr), collapse = ", "), st$ard), "ARD")
    else add("caution", sprintf("No program map was shrunk by ARD (ard = %g, %d factors): n_factors may be too small to hold all structure.",
                                st$ard, length(strength)), "ARD")
  }
  # few-gene / low-count programs among the strongest factors
  for (p in programs) {
    if (p$kind != "factor" || p$strength < 0.1 * max(strength)) next
    if (!is.null(progtab) && length(tests) && !(p$id %in% progtab$program_map[!is.na(progtab$program)])) next
    v <- p$loading - mean(p$loading); top2 <- sort(v^2, decreasing = TRUE)
    if (length(top2) >= 3 && sum(top2[1:2]) / sum(top2) > 0.7 && !isTRUE(sig[[p$id]]))
      add("caution", sprintf("%s is dominated by %s: a single-gene axis may be technical.", lbp(p$id),
                             paste(names(top2)[1:2], collapse = " and ")), "Few-gene program map")
    tg <- .rr_top(p$loading, 5)
    if (length(gene_mean) >= 10 && stats::median(gene_mean[tg]) < stats::quantile(gene_mean, 0.2))
      add("caution", sprintf("The top genes of %s are low-count genes (median %.3g transcripts per bin): loadings of sparse genes are noisy.",
                             lbp(p$id), stats::median(gene_mean[tg])), "Low counts")
  }
  # fit convergence
  conv <- fit$convergence
  if (!is.null(conv)) {
    it <- if (!is.null(fit$iterations)) fit$iterations[1] else NA
    if (conv == 0) add("ok", sprintf("Optimiser converged (%s function evaluations).", it), "Fit")
    else if (conv == 1) add("caution", sprintf("The fit stopped at max_iter = %s (normal for this model): loadings and strengths depend on max_iter, so compare fits and bootstrap refits only at the same max_iter.",
                                               .rr_or(st$max_iter, "?")), "Fit stopped at max_iter")
    else add("warning", sprintf("Optimiser returned code %s (%s).", conv, .rr_or(fit$message, "")), "Fit")
  }
  if (N < 500) add("caution", sprintf("Only %d in-tissue bins: program maps and tests have little data.", N), "Few bins")
  if (J < 10) add("caution", sprintf("Only %d genes: programs cannot be told apart well from cellularity.", J), "Few genes")
  if (is.finite(mask$boundary) && mask$boundary > 0.25)
    add("caution", sprintf("Irregular or fragmented tissue mask (%.0f%% of in-tissue bins lie on an internal edge): the gene-shift null moves genes across the edge and borderline p-values are less reliable.",
                           100 * mask$boundary), "Tissue mask")
  # offset
  if (identical(off_desc$type, "area"))
    add("info", "No known structure was removed (offset = \"area\"): the program maps describe all spatial structure, including cell-type composition and domains. Use rff_offset() for a residual analysis.",
        "Offset")
  if (off_desc$data_driven)
    add("caution", sprintf("The offset was estimated from the same data (method = %s%s): the richer it is, the more of any program it absorbs, and structure it misses is called like any other. Check neighbouring values of k.",
                           off_desc$method, if (!is.null(off_desc$k)) paste0(", k = ", off_desc$k) else ""), "Data-driven offset")
  if (identical(off_desc$type, "matrix") && is.na(off_desc$method))
    add("info", "A matrix offset was used, but how it was made is not recorded; pass offset_info = list(method = ..., covariates = ...) to document it.", "Offset")
  if (maps$f > 1) add("info", sprintf("Maps are block-averaged by %d x %d bins for display.", maps$f, maps$f), "Display")
  if (!is.null(cellres) && cellres$n_outside > 0)
    add("caution", sprintf("%d cells lie outside the fitted tissue: their field values are extrapolated.", cellres$n_outside), "Cells")
  F
}

# ---- reproduce -------------------------------------------------------------------------------

.rr_is_call_to <- function(x, fn) {
  is.call(x) && (as.character(x[[1]])[length(as.character(x[[1]]))] %in% fn)
}

.rr_dep <- function(x) paste(deparse(x, width.cutoff = 80), collapse = "\n")

.rr_fit_call <- function(fit, st, off_desc) {
  off_arg <- switch(off_desc$type, area = NULL, smoothed_total = "\"smoothed_total\"", matrix = "off")
  a <- c(sprintf("n_factors = %s", .rr_or(st$n_factors, ncol(fit$L))),
         if (!is.null(st$lengthscales)) sprintf("lengthscales = %s", .rr_dep(st$lengthscales)),
         if (!is.null(off_arg)) sprintf("offset = %s", off_arg),
         if (isTRUE(fit$has_density)) sprintf("density_lengthscale = %s", .rr_dep(st$density_lengthscale)) else if (identical(off_desc$type, "area")) "density_lengthscale = NULL",
         if (!is.null(st$ard)) sprintf("ard = %s", st$ard),
         if (!is.null(st$basis)) sprintf("basis = \"%s\"", st$basis),
         if (!identical(st$basis, "grid") && !is.null(st$n_features)) sprintf("n_features = %s", st$n_features),
         if (!is.null(st$learn_lengthscales)) sprintf("learn_lengthscales = %s", st$learn_lengthscales),
         if (!is.null(st$l1) && length(st$l1) == 1) sprintf("l1 = %s", st$l1),
         if (!is.null(st$family)) sprintf("family = \"%s\"", st$family),
         if (!is.null(st$factor_init)) sprintf("factor_init = \"%s\"", st$factor_init),
         if (!is.null(st$max_iter)) sprintf("max_iter = %s", st$max_iter),
         if (!is.null(st$seed)) sprintf("seed = %s", st$seed))
  paste0("fit <- fit_spatial_rff(binned,\n  ", paste(a, collapse = ",\n  "), ")")
}

.rr_reproduce <- function(cl, fit, tests, off_desc, offset_info, out_dir, file, st, g, res_obj = NULL) {
  if (!is.null(res_obj)) {
    rc <- res_obj$call; rc[[1]] <- as.name("rff_programs")
    if (!is.null(rc$binned)) rc$binned <- as.name("binned")
    cc <- cl; cc$fit <- as.name("res"); if (!is.null(cc$binned)) cc$binned <- as.name("binned"); cc[[1]] <- as.name("rff_report")
    return(paste(c("library(cohalu)",
                   "# 1. Bin the transcripts (the report was made from this grid)",
                   sprintf("binned <- bin_transcripts(tx, bin_size = %g)   # %d x %d bins", g$bin_size, g$nx, g$ny),
                   "# 2. Offset, fit, program test and Programs table in one call",
                   paste0("res <- ", .rr_dep(rc)),
                   if (!is.null(res_obj$fit_call)) c("#    (the fit it made is equivalent to:)", paste0("#    ", strsplit(.rr_dep(res_obj$fit_call), "\n")[[1]])),
                   if (!is.null(offset_info$call)) paste0("#    offset: ", if (is.character(offset_info$call)) offset_info$call else .rr_dep(offset_info$call)),
                   "# 3. Report", .rr_dep(cc),
                   "# The exact fit used here is saved with the report:",
                   sprintf("# fit <- readRDS(\"%s\")", file.path(basename(out_dir), "fit.rds"))), collapse = "\n"))
  }
  lines <- c("library(cohalu)",
             "# 1. Bin the transcripts (the report was made from this grid)",
             sprintf("binned <- bin_transcripts(tx, bin_size = %g)   # %d x %d bins", g$bin_size, g$nx, g$ny))
  if (identical(off_desc$type, "matrix")) {
    oc <- offset_info$call
    lines <- c(lines, "# 2. Known structure as offset",
               if (!is.null(oc)) paste0("off <- ", if (is.character(oc)) oc else .rr_dep(oc)) else
                 if (!is.na(off_desc$method) && off_desc$method %in% c("kmeans", "pca"))
                   sprintf("off <- rff_offset(binned, method = \"%s\", k = %s)", off_desc$method, .rr_or(off_desc$k, 6)) else
                   "off <- rff_offset(binned, covariates)   # covariates: not recorded in the fit")
  }
  fit_expr <- cl$fit
  lines <- c(lines, "# 3. Fit",
             if (.rr_is_call_to(fit_expr, "fit_spatial_rff")) paste0("fit <- ", .rr_dep(fit_expr)) else .rr_fit_call(fit, st, off_desc))
  te <- cl$test
  if (length(tests)) {
    lines <- c(lines, "# 4. Test")
    if (.rr_is_call_to(te, c("rff_program_test", "rff_factor_test"))) lines <- c(lines, paste0("test <- ", .rr_dep(te))) else {
      for (i in seq_along(tests)) {
        t <- tests[[i]]
        lines <- c(lines, if (t$kind == "factor")
          sprintf("test%d <- rff_factor_test(fit, binned, n_boot = %s%s)", i, t$n_boot, if (t$sequential) ", sequential = TRUE" else "") else
            sprintf("test%d <- rff_program_test(fit, binned, n_boot = %s, null = \"%s\"%s%s)", i, t$n_boot, t$null_type,
                    if (!is.null(t$highpass)) sprintf(", highpass = %g", t$highpass) else "",
                    if (!t$sequential) ", sequential = FALSE" else ""))
      }
      if (length(tests) > 1) lines <- c(lines, sprintf("test <- list(%s)", paste0("test", seq_along(tests), collapse = ", ")))
      else lines[length(lines)] <- sub("^test1 <-", "test <-", lines[length(lines)])
    }
  }
  cc <- cl; cc$fit <- as.name("fit"); if (!is.null(cc$test)) cc$test <- as.name("test")
  if (!is.null(cc$binned)) cc$binned <- as.name("binned")
  cc[[1]] <- as.name("rff_report")
  lines <- c(lines, "# 5. Report", .rr_dep(cc),
             "# The exact fit used here is saved with the report:",
             sprintf("# fit <- readRDS(\"%s\")", file.path(basename(out_dir), "fit.rds")))
  paste(lines, collapse = "\n")
}

# ---- export ------------------------------------------------------------------------------------

.rr_export <- function(out_dir, fit, binned, keep, genes, programs, tests, comp_scores, cellres, known,
                       gene_tot, gene_mean, dev_expl, sig, pmin_of, repro, st, cl, alpha, g, progtab = NULL) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  files <- list()
  reg <- function(name, desc, cols) {
    path <- file.path(out_dir, name)
    files[[length(files) + 1]] <<- list(name = name, path = path, desc = desc, cols = cols, size = file.size(path))
  }
  wcsv <- function(df, name) utils::write.csv(df, file.path(out_dir, name), row.names = FALSE)
  if (!is.null(progtab)) {
    pt <- progtab; attr(pt, "loadings") <- NULL; attr(pt, "fields") <- NULL
    wcsv(pt, "programs.csv")
    reg("programs.csv", "The main result: one row per program (Confirmed / Candidate / Exploratory), plus strong program maps not supported by the test and cellularity-like maps.",
        "program (P1, P2, ...; NA for secondary rows), status, detected_by, p (from the detection axis, or the factor test), p_factor_test, detection_axis (PC of rff_program_test), program_map (factor of the fit), r (|r| between axis score and map field), best_map / best_r, lengthscale, program_strength, uniform_share, variance_share, n_genes (|loading| >= 25% of the largest), top_up / top_down, key, map_label (M*), axis_label (D*). Loadings per program: the P columns of loadings.csv; maps: the P columns of fields.csv")
  }
  fac <- vapply(programs, function(p) p$kind == "factor", TRUE)
  pf <- programs[fac]
  fs <- data.frame(factor = names(pf), display_name = vapply(pf, `[[`, "", "label"), strength = vapply(pf, `[[`, 0, "strength"),
                   program_strength = vapply(pf, `[[`, 0, "program_strength"),
                   uniform_share = vapply(pf, `[[`, 0, "uniform_share"),
                   lengthscale = vapply(pf, `[[`, 0, "lengthscale"),
                   p_min = unname(pmin_of[names(pf)]), significant = unname(sig[names(pf)]),
                   top_up = vapply(pf, function(p) paste(.rr_top(p$loading, 8, 1), collapse = ";"), ""),
                   top_down = vapply(pf, function(p) paste(.rr_top(p$loading, 8, -1), collapse = ";"), ""))
  wcsv(fs, "factors.csv")
  reg("factors.csv", "One row per fitted factor (program map).", "factor, display_name (M1, M2, ... as shown in the report), strength (loading norm), program_strength (norm after removing the all-gene mean), uniform_share (share common to all genes; near 1 = cellularity-like), lengthscale, p_min (smallest p over the tests supplied), significant (at alpha), top_up / top_down (genes with the largest positive / negative loadings)")
  ld <- data.frame(gene = genes, total_counts = unname(gene_tot), mean_per_bin = unname(gene_mean),
                   alpha = unname(fit$alpha[genes]),
                   dispersion = if (!is.null(fit$dispersion)) unname(fit$dispersion[genes]) else NA,
                   offset_deviance_explained = if (!is.null(dev_expl)) unname(dev_expl) else NA,
                   do.call(cbind, lapply(programs, function(p) stats::setNames(data.frame(unname(p$loading)), p$id))),
                   check.names = FALSE)
  wcsv(ld, "loadings.csv")
  reg("loadings.csv", "One row per gene.", "gene, total_counts, mean_per_bin (in-tissue), alpha (log intercept), dispersion (NB), offset_deviance_explained (share of Poisson deviance explained by the offset; NA without a matrix offset), then one column per factor (program maps M1.., loadings L on the log scale) and per program-test PC (detection axes D1.., unit-norm PCA loadings); see display_names.csv")
  fg <- fit$field_grid
  fl <- data.frame(bin = which(keep), x = binned$coords$x[keep], y = binned$coords$y[keep],
                   total_counts = Matrix::rowSums(binned$counts[keep, , drop = FALSE]),
                   fg[keep, , drop = FALSE], check.names = FALSE)
  for (id in names(comp_scores)) fl[[id]] <- comp_scores[[id]][keep]
  wcsv(fl, "fields.csv")
  reg("fields.csv", "One row per in-tissue bin.", "bin (index into the full nx x ny grid, x fastest), x, y (bin centre), total_counts, density (cellularity field f0, if fitted), factor1..K (unit-variance fields = program maps M1..MK), and the PC scores of rff_program_test() (detection axes D1..) if given; see display_names.csv")
  if (!is.null(cellres)) {
    wcsv(cellres$table, "cells.csv")
    reg("cells.csv", "One row per supplied cell.", "cell coordinates and label as given, then the field values at the cell (rff_fields(): density, factor1..K = program maps) and detection-axis scores (PC1..) (value of the bin containing the cell)")
  }
  for (i in seq_along(tests)) {
    t <- tests[[i]]
    tb <- t$table; tb$statistic <- as.numeric(t$obs)
    tb <- data.frame(tb[1], display_name = vapply(tb$id, function(i) .rr_or(programs[[i]]$label, i), ""), tb[-1], check.names = FALSE)
    nm <- sprintf("test%d_%s.csv", i, t$kind)
    wcsv(tb, nm)
    reg(nm, sprintf("Results of %s.", t$label), if (t$kind == "factor") "id / factor, display_name (M1, ...), p (max-null), p_sequential, p_use (p used for the flag), top_genes, sig, statistic (observed program strength)" else
      "id, display_name (D1, ...), component, share (variance share), p, p_use, factor / r_factor (best-matching program map, as fit factor name, and |r|), top_genes, sig, statistic (= share)")
    nl <- do.call(rbind, lapply(seq_along(t$nulls), function(k) if (length(t$nulls[[k]]))
      data.frame(id = t$table$id[k], display_name = .rr_or(programs[[t$table$id[k]]]$label, t$table$id[k]), draw = seq_along(t$nulls[[k]]), null_statistic = t$nulls[[k]])))
    if (!is.null(nl)) {
      nm2 <- sprintf("test%d_null.csv", i); wcsv(nl, nm2)
      reg(nm2, sprintf("Null statistics of %s (the distribution each program was compared with).", t$label), "id, display_name, draw, null_statistic")
    }
  }
  if (!is.null(known$table)) {
    kt <- known$table
    kt <- data.frame(kt[1], display_name = vapply(kt$program, function(i) .rr_or(programs[[i]]$label, i), ""), kt[-1], check.names = FALSE)
    wcsv(kt, "known_programs.csv")
    reg("known_programs.csv", "Overlap of every program map / detection axis with every known program.", "program, display_name, known_program, n_genes_in_panel, n_genes_given, cor_loading_membership (correlation of the loadings with set membership), top_overlap (set genes among the top_n genes by |loading|), top_n, p_overlap (hypergeometric), mean_loading, overlapping_genes")
  }
  dn <- data.frame(id = names(programs), display_name = vapply(programs, `[[`, "", "label"),
                   kind = c(factor = "program map", component = "detection axis", program = "program")[vapply(programs, `[[`, "", "kind")],
                   package_name = vapply(programs, `[[`, "", "pkg"))
  wcsv(dn, "display_names.csv")
  reg("display_names.csv", "How the names shown in the report map to the R objects.", "id (column name in loadings.csv / fields.csv / cells.csv and in the R objects: factorK of the fit, PCk of rff_program_test(), Pk = merged program), display_name (P = program, M = program map, D = detection axis), kind, package_name (what the program is made of)")
  saveRDS(fit, file.path(out_dir, "fit.rds"))
  reg("fit.rds", "The spatial_rff_fit object (readRDS()).", "all estimates, fields on the grid, settings")
  sj <- list(fit_settings = st, family = fit$family, offset = .rr_or(fit$offset, "area"), basis = .rr_or(fit$basis, "random"),
             convergence = fit$convergence, message = fit$message, iterations = fit$iterations,
             n_bins = sum(keep), n_genes = length(genes), grid = g, alpha = alpha,
             tests = lapply(tests, function(t) list(label = t$label, kind = t$kind, null = t$null_type, sequential = t$sequential,
                                                     n_boot = t$n_boot, highpass = t$highpass, elapsed_sec = t$elapsed)),
             report_call = .rr_dep(cl))
  writeLines(.rr_json(sj), file.path(out_dir, "settings.json"))
  reg("settings.json", "Fit and test settings, convergence, report call (JSON).", "fit_settings, family, offset, basis, convergence, iterations, grid, tests, report_call")
  writeLines(repro, file.path(out_dir, "reproduce.R"))
  reg("reproduce.R", "R code of the analysis (as in 'How to reproduce').", "-")
  files
}

# ---- static HTML --------------------------------------------------------------------------------

.rr_kv <- function(k, v) paste0("<tr><th>", .rr_h(k), "</th><td>", v, "</td></tr>")
.rr_table <- function(rows, cls = "kv") paste0("<table class=\"", cls, "\">", paste(rows, collapse = ""), "</table>")

.rr_html_overview <- function(N, J, K, g, fit, tests, sig, alpha, maps) {
  nsig <- sum(sig, na.rm = TRUE)
  tile <- function(v, l, s = "") sprintf("<div class=\"tile\"><div class=\"tv\">%s</div><div class=\"tl\">%s</div>%s</div>", v, .rr_h(l),
                                         if (nzchar(s)) sprintf("<div class=\"ts\">%s</div>", .rr_h(s)) else "")
  tt <- if (length(tests)) paste(unique(vapply(tests, function(t) if (t$kind == "factor") "factor test" else paste0(t$null_type, " null"), "")), collapse = ", ") else "none"
  paste0("<div class=\"tiles\">",
         tile(format(N, big.mark = ","), "in-tissue bins", sprintf("%g x %g grid of %g", g$nx, g$ny, g$bin_size)),
         tile(J, "genes"), tile(K, "program maps", sprintf("basis %s, %s", .rr_or(fit$basis, "random"), fit$family)),
         tile(if (length(tests)) nsig else "-", "programs (confirmed + candidate)", sprintf("test: %s, alpha %g", tt, alpha)),
         tile(.rr_num(N * g$bin_size^2 / 1e6, 3), "tissue area", "mm2 if coordinates are in um"),
         "</div>")
}

.rr_html_input <- function(fit, binned, keep, genes, gene_tot, bin_tot, dev_expl, off_desc, offset_info,
                           cells, cellres, known_programs, meta, st, tests, g, cl, mask) {
  N <- sum(keep)
  dropped <- setdiff(binned$genes, genes)
  out <- character(0)
  if (!is.null(meta) && length(meta)) {
    out <- c(out, "<h3>Sample</h3>", .rr_table(mapply(function(k, v) .rr_kv(k, .rr_h(paste(v, collapse = ", "))),
                                                        names(meta), meta)))
  }
  gt <- sort(gene_tot, decreasing = TRUE)
  data_rows <- c(
    .rr_kv("Transcripts in the bins", sprintf("%s (all genes of the binned object); %s in the %d modelled genes, in-tissue",
                                              format(sum(binned$counts), big.mark = ","), format(sum(gene_tot), big.mark = ","), length(genes))),
    .rr_kv("Genes", sprintf("%d modelled%s <details><summary>list</summary><div class=\"genelist\">%s</div></details>", length(genes),
                            if (length(dropped)) sprintf("; %d in the binned object not modelled (%s)", length(dropped),
                                                         .rr_h(paste(utils::head(dropped, 20), collapse = ", "))) else "",
                            .rr_h(paste(genes, collapse = ", ")))),
    .rr_kv("Grid", sprintf("%d x %d bins of %g (coordinate units); origin (%g, %g)", g$nx, g$ny, g$bin_size, g$xmin, g$ymin)),
    .rr_kv("In-tissue bins", sprintf("%s of %s (%.0f%%); tissue area %s coordinate units^2 (%s mm^2 if um); internal mask edge %.0f%% of bins",
                                     format(N, big.mark = ","), format(g$nx * g$ny, big.mark = ","), 100 * N / (g$nx * g$ny),
                                     format(round(N * g$bin_size^2), big.mark = ","), .rr_num(N * g$bin_size^2 / 1e6, 3),
                                     100 * .rr_or(mask$boundary, NA))),
    .rr_kv("Counts per in-tissue bin", sprintf("mean %.2f, median %g, 5-95%% %g-%g, max %g; %.0f%% empty <div id=\"hist_counts\" class=\"minihist\"></div>",
                                              mean(bin_tot), stats::median(bin_tot), stats::quantile(bin_tot, .05), stats::quantile(bin_tot, .95),
                                              max(bin_tot), 100 * mean(bin_tot == 0))),
    .rr_kv("Top genes by counts", .rr_h(paste(sprintf("%s (%s)", names(gt)[1:min(10, length(gt))], format(gt[1:min(10, length(gt))], big.mark = ",")), collapse = ", ")))
  )
  out <- c(out, "<h3>Transcripts and bins</h3>", .rr_table(data_rows))
  # offset
  orow <- c(.rr_kv("Type", .rr_h(off_desc$text)))
  oi <- .rr_or(offset_info, list())
  if (!is.null(oi$covariates)) {
    cv <- as.data.frame(oi$covariates)
    desc <- vapply(names(cv), function(nm) {
      v <- cv[[nm]]
      if (is.numeric(v)) sprintf("%s: numeric, mean %.3g, sd %.3g, range %.3g-%.3g", nm, mean(v, na.rm = TRUE), stats::sd(v, na.rm = TRUE),
                                 min(v, na.rm = TRUE), max(v, na.rm = TRUE)) else {
        tb <- sort(table(v), decreasing = TRUE)
        sprintf("%s: %d levels (%s)", nm, length(tb), paste(sprintf("%s %d", names(tb)[1:min(8, length(tb))], tb[1:min(8, length(tb))]), collapse = ", "))
      }
    }, "")
    orow <- c(orow, .rr_kv("Covariates", paste(.rr_h(desc), collapse = "<br>")))
  }
  for (nm in setdiff(names(oi), c("covariates", "method", "k", "call")))
    orow <- c(orow, .rr_kv(nm, .rr_h(paste(format(oi[[nm]]), collapse = ", "))))
  if (!is.null(oi$call)) orow <- c(orow, .rr_kv("Call", paste0("<code>", .rr_h(if (is.character(oi$call)) oi$call else .rr_dep(oi$call)), "</code>")))
  if (!is.null(dev_expl)) {
    de <- sort(dev_expl, decreasing = TRUE)
    orow <- c(orow, .rr_kv("Deviance explained by the offset",
                           sprintf("per gene (Poisson, vs a constant rate): median %.2f, IQR %.2f-%.2f. Highest: %s. Lowest: %s.",
                                   stats::median(de, na.rm = TRUE), stats::quantile(de, .25, na.rm = TRUE), stats::quantile(de, .75, na.rm = TRUE),
                                   .rr_h(paste(sprintf("%s %.2f", names(de)[1:min(5, length(de))], de[1:min(5, length(de))]), collapse = ", ")),
                                   .rr_h(paste(sprintf("%s %.2f", rev(names(de))[1:min(5, length(de))], rev(de)[1:min(5, length(de))]), collapse = ", ")))))
  }
  out <- c(out, "<h3>Offset / known structure</h3>", .rr_table(orow))
  if (!is.null(cellres)) {
    crow <- c(.rr_kv("Cells", sprintf("%s (%s drawn on the map)%s", format(cellres$json$n, big.mark = ","), format(cellres$json$n_drawn, big.mark = ","),
                                      if (cellres$n_outside) sprintf("; %d outside the fitted tissue", cellres$n_outside) else "")),
              .rr_kv("Labels", if (is.null(cellres$json$label_col)) "none given" else
                sprintf("column <code>%s</code>: %s", .rr_h(cellres$json$label_col),
                        .rr_h(paste(sprintf("%s %s", names(cellres$counts), cellres$counts), collapse = ", ")))))
    out <- c(out, "<h3>Cells</h3>", .rr_table(crow))
  }
  if (!is.null(known_programs)) {
    krow <- vapply(names(known_programs), function(k) {
      s <- known_programs[[k]]; inn <- intersect(s, genes)
      .rr_kv(k, sprintf("%d genes, %d in the panel: %s", length(s), length(inn), .rr_h(paste(utils::head(inn, 25), collapse = ", "))))
    }, "")
    out <- c(out, "<h3>Known programs</h3>", .rr_table(krow))
  }
  # arguments
  fmt <- function(v) if (is.null(v)) "NULL" else if (length(v) > 6) sprintf("%s ... (%d values)", paste(format(utils::head(v, 6)), collapse = ", "), length(v)) else paste(format(v), collapse = ", ")
  arow <- c(vapply(names(st), function(k) .rr_kv(k, .rr_h(fmt(st[[k]]))), ""),
            .rr_kv("offset", .rr_h(.rr_or(fit$offset, "area"))),
            .rr_kv("fitted length scales", .rr_h(paste(signif(fit$lengthscales, 3), collapse = ", "))),
            if (isTRUE(fit$has_density)) .rr_kv("cellularity field", sprintf("sigma0 %.3g, length scale %.3g", fit$sigma0, fit$density_lengthscale)),
            .rr_kv("convergence", .rr_h(sprintf("%s (%s); %s evaluations", .rr_or(fit$convergence, "?"), .rr_or(fit$message, ""),
                                                paste(.rr_or(fit$iterations, "?"), collapse = "/")))))
  out <- c(out, "<h3>Fit arguments</h3>", .rr_table(arow))
  if (length(tests)) {
    trow <- vapply(tests, function(t) .rr_kv(t$label, .rr_h(sprintf("%s; null %s; %s; n_boot %s%s%s%s", if (t$kind == "factor") "rff_factor_test" else "rff_program_test",
                                                                   t$null_type, if (t$sequential) "sequential" else "not sequential", t$n_boot,
                                                                   if (!is.null(t$highpass)) sprintf("; highpass %g", t$highpass) else "",
                                                                   if (!is.null(t$bandwidth)) sprintf("; bandwidth %g", t$bandwidth) else "",
                                                                   if (!is.null(t$elapsed)) sprintf("; %.1f s", t$elapsed) else ""))), "")
    out <- c(out, "<h3>Test arguments</h3>", .rr_table(trow))
  } else out <- c(out, "<h3>Test arguments</h3><p class=\"mut\">No test supplied.</p>")
  args <- as.list(cl)[-1]
  out <- c(out, "<h3>Report arguments</h3>", .rr_table(vapply(names(args), function(k) .rr_kv(k, paste0("<code>", .rr_h(substr(.rr_dep(args[[k]]), 1, 300)), "</code>")), "")))
  paste(out, collapse = "\n")
}

.rr_html_findings <- function(F) {
  lab <- c(ok = "OK", caution = "Caution", warning = "Warning", info = "Info")
  ico <- c(ok = "&#10003;", caution = "!", warning = "&#9888;", info = "i")
  rk <- match(vapply(F, `[[`, "", "level"), c("warning", "ok", "caution", "info"))
  rk[vapply(F, function(f) isTRUE(f$pin), TRUE)] <- 0
  ord <- order(rk)
  items <- vapply(F[ord], function(f) sprintf("<li class=\"fd %s\"><span class=\"badge %s\"><span class=\"bi\">%s</span>%s</span><div>%s%s</div></li>",
                                               f$level, f$level, ico[[f$level]], lab[[f$level]],
                                               if (!is.null(f$head)) sprintf("<b>%s.</b> ", .rr_h(f$head)) else "", .rr_h(f$text)), "")
  paste0("<ul class=\"findings\">", paste(items, collapse = ""), "</ul>")
}

.rr_html_outputs <- function(files, out_dir, file, export) {
  if (!export) return("<p class=\"mut\">No files were exported (<code>export = FALSE</code>).</p>")
  rel <- basename(out_dir)
  sz <- function(b) if (is.na(b)) "?" else if (b > 1e6) sprintf("%.1f MB", b / 1e6) else sprintf("%.1f kB", b / 1e3)
  rows <- vapply(files, function(f) sprintf("<tr><td><a href=\"%s\">%s</a></td><td class=\"num\" style=\"white-space:nowrap\">%s</td><td>%s<div class=\"mut small\">%s</div></td></tr>",
                                             .rr_h(paste0(rel, "/", f$name)), .rr_h(f$name), sz(f$size), .rr_h(f$desc), .rr_h(f$cols)), "")
  paste0("<p>Written next to this report in <code>", .rr_h(out_dir), "</code> (links are relative to the report: keep the folder beside it). The report itself is <code>",
         .rr_h(file), "</code>.</p><table class=\"grid\"><thead><tr><th>File</th><th>Size</th><th>Contents and columns</th></tr></thead><tbody>",
         paste(rows, collapse = ""), "</tbody></table>")
}

.rr_html_methods <- function(fit, tests) {
  paste0(
    "<p>For in-tissue bin <i>b</i> (centre <i>u</i>) and gene <i>j</i>, counts are modelled as ",
    if (identical(fit$family, "nb")) "negative binomial" else "Poisson", " with</p>",
    "<p class=\"eq\">log E[y<sub>bj</sub>] = log(area) + offset<sub>bj</sub> + &alpha;<sub>j</sub> + &sigma;<sub>0</sub> f<sub>0</sub>(u<sub>b</sub>) + &Sigma;<sub>k</sub> L<sub>jk</sub> f<sub>k</sub>(u<sub>b</sub>)</p>",
    "<p>The offset holds structure already known (cell-type composition, domains, an embedding; <code>rff_offset()</code>), f<sub>0</sub> is a cellularity field shared by all genes, and f<sub>1</sub>..f<sub>K</sub> are smooth spatial factors - shown as <b>program maps</b> M1..MK - (Gaussian processes with RBF covariance of length scale &ell;<sub>k</sub>, ",
    if (identical(fit$basis, "grid")) "on the Fourier basis of the bin grid" else "approximated by random Fourier features", "). ",
    "Each field is scaled to unit variance over the tissue, so the loading L<sub>jk</sub> is the change of gene j's log intensity per standard deviation of the field: a loading of 0.7 means the gene is about 2-fold (e<sup>0.7</sup>) higher where the field is +1 sd. ",
    "<b>Program strength</b> is the norm of a factor's loadings after removing their mean over genes; the removed part moves all genes together and cannot be told apart from cellularity (<b>uniform share</b> near 1 = cellularity-like).</p>",
    "<h3>Detection axes vs program maps</h3>", .rr_html_fvp(open = TRUE),
    "<h3>How to read the panels</h3><ul>",
    "<li><b>Program map and detection axis tables</b>: click a row to select a program map (M, fit factor) or detection axis (D, program test PC) everywhere; click a header to sort. Significance comes only from a test.</li>",
    "<li><b>Maps</b>: red = high field (genes with positive loading are up), blue = low. Drag to pan, wheel to zoom, double-click to reset; hover for values. The gene view shows observed counts, the null expectation (offset + intercept + cellularity, without factors), their log2 ratio (what the program maps and noise must explain) and the fitted program-map contribution &Sigma;<sub>k</sub> L<sub>jk</sub> f<sub>k</sub>.</li>",
    "<li><b>Loadings</b>: the heatmap is coloured per column (scaled to the column's largest |loading|); hover for raw values. Factor loadings are on the log scale; component loadings of <code>rff_program_test()</code> are unit-norm PCA loadings of smoothed residuals.</li>",
    "<li><b>Cells</b>: distribution (5, 25, 50, 75, 95% and mean) of each field at the cells of each label.</li>",
    "<li><b>Known programs</b>: correlation of the loadings with membership of each known set and overlap with the top genes; 'known' needs |r| &ge; 0.3 and a significant overlap of the top genes (hypergeometric p &lt; 0.05), 'partial' one of the two, 'novel' neither.</li>",
    "<li><b>Tests</b>: the observed statistic (line) against the null distribution (bars); p = (1 + #null &ge; observed) / (n + 1).</li></ul>",
    "<p class=\"mut\">Reference: Gundersen GW, Zhang MM, Engelhardt BE (2021). Latent variable modeling with random features. AISTATS.</p>")
}

.rr_html_session <- function(pkg_version, t_start) {
  si <- utils::sessionInfo()
  oth <- if (length(si$otherPkgs)) paste(vapply(si$otherPkgs, function(p) paste(p$Package, p$Version), ""), collapse = ", ") else "none"
  sprintf("%s; %s; cohalu %s; attached: %s; %d namespaces loaded; report built in %.1f s.",
          .rr_h(si$R.version$version.string), .rr_h(si$platform), .rr_h(pkg_version), .rr_h(oth),
          length(si$loadedOnly) + length(si$otherPkgs), as.numeric(difftime(Sys.time(), t_start, units = "secs")))
}

.rr_html_fvp <- function(open = FALSE) {
  paste0(
    "<details class=\"help\"", if (open) " open" else "", "><summary>Detection axes vs program maps: what is the difference?</summary><div class=\"helpbody\">",
    "<p><b>Program map</b> (M1, M2, ...; package name: factor) = RFLVM component of <code>fit_spatial_rff()</code>: a smooth field f<sub>k</sub>(u) plus gene loadings L<sub>jk</sub>; its p-value comes from <code>rff_factor_test()</code> (refits to data simulated without factors, parametric bootstrap). ",
    "<b>Detection axis</b> (D1, D2, ...; package name: PC of <code>rff_program_test()</code>) = principal component of the Gaussian-smoothed Pearson residuals under the RFLVM null mean (offset, gene intercepts, cellularity field and the gene-shared part of each factor), after removing the shared multiplicative direction and standardising genes; its p-value comes from comparing its variance share with null data (parametric simulation or gene-shift surrogates) in <code>rff_program_test()</code>.</p>",
    "<table class=\"grid\"><thead><tr><th></th><th>Program map M (fit_spatial_rff factor)</th><th>Detection axis D (rff_program_test PC)</th></tr></thead><tbody>",
    "<tr><th>Role</th><td>Interpretation: smooth map and log-scale gene effects</td><td>Detection: whether and how many programs exist</td></tr>",
    "<tr><th>Strengths</th><td>Denoised field at any point (cells via rff_fields()); loadings are fold changes; length scale</td><td>Better calibrated (gene-shift null keeps each gene's own autocorrelation); more powerful in simulations; no refitting</td></tr>",
    "<tr><th>Weaknesses</th><td>Depends on max_iter, ARD and initialisation; its test is anti-conservative with gene-specific structure</td><td>Unit-norm loadings on smoothed residuals (no fold change); map only on bins; components can mix programs</td></tr>",
    "</tbody></table>",
    "<ol><li>Count the significant detection axes: that is the number of programs supported.</li>",
    "<li>For each significant detection axis, open its best-matching program map (r_factor), e.g. 'D1 is significant &rarr; read it through its best-matching program map M2 (r = 0.81)': if |r| &ge; 0.5, read the program through that map and its loadings.</li>",
    "<li>If |r| &lt; 0.5, the RFLVM did not capture it well: interpret the detection-axis loadings and score map directly. Program maps without a matching significant detection axis are not supported.</li></ol>",
    "</div></details>")
}
