# Nonparametric bootstrap (Section 2.6 of Locher, Richardson, and Stensrud,
# 2026). All thresholds, boundary points, and bounds are smooth functions of
# the identified quantities, so it suffices to bootstrap the identified
# quantities once and re-evaluate each statistic on every replicate.

#' Bootstrap replicates of the identified quantities
#'
#' @param y,a,z Outcome, treatment, and instrument vectors.
#' @param n_boot Number of bootstrap replicates.
#' @param strata `"none"` for iid resampling of observations, `"instrument"`
#'   to resample within instrument arms (keeping the arm sizes fixed, as in a
#'   randomised trial).
#' @return A matrix with `n_boot` rows and one column per identified quantity.
#' @keywords internal
#' @noRd
bootstrap_q <- function(y, a, z, n_boot, strata = "none") {
  g <- iv_cells(a, z)
  n <- length(y)
  groups <- if (strata == "instrument") split(seq_len(n), z) else list(seq_len(n))

  out <- matrix(NA_real_, n_boot, length(q_names), dimnames = list(NULL, q_names))
  for (b in seq_len(n_boot)) {
    idx <- unlist(lapply(groups, function(i) i[sample.int(length(i), length(i), replace = TRUE)]),
                  use.names = FALSE)
    cs <- cell_stats(y[idx], g[idx])
    out[b, ] <- identify_from_stats(cs$n, cs$s)
  }
  out
}

#' Estimate a statistic and its bootstrap percentile confidence interval
#'
#' @param stat_fun Function of the identified quantities `q` returning a named
#'   numeric vector.
#' @param q Identified quantities (point estimates).
#' @param boot Matrix of bootstrap replicates of `q`, or `NULL`.
#' @param level Confidence level.
#' @return A data frame with columns `term`, `estimate` and, if `boot` is
#'   supplied, `ci_lower` and `ci_upper`. The attribute `n_dropped` gives the
#'   number of replicates in which some statistic could not be evaluated
#'   (e.g. because a compliance type was absent from the resample).
#' @keywords internal
#' @noRd
evaluate_stat <- function(stat_fun, q, boot = NULL, level = 0.95) {
  est <- stat_fun(q)
  out <- data.frame(term = names(est), estimate = unname(est), stringsAsFactors = FALSE)
  if (is.null(boot)) return(out)

  k <- length(est)
  reps <- matrix(NA_real_, nrow(boot), k)
  for (b in seq_len(nrow(boot))) {
    reps[b, ] <- tryCatch(
      { val <- stat_fun(boot[b, ]); if (length(val) == k) val else rep(NA_real_, k) },
      error = function(e) rep(NA_real_, k)
    )
  }

  alpha <- (1 - level) / 2
  pct <- function(v, p) if (any(is.finite(v))) stats::quantile(v[is.finite(v)], p, names = FALSE) else NA_real_
  out$ci_lower <- apply(reps, 2, pct, p = alpha)
  out$ci_upper <- apply(reps, 2, pct, p = 1 - alpha)

  # Only count failures for statistics that are defined at the point estimate.
  defined <- is.finite(est)
  attr(out, "n_dropped") <- sum(!apply(is.finite(reps[, defined, drop = FALSE]), 1, all))
  attr(out, "replicates") <- reps
  out
}

#' Warn once if bootstrap replicates had to be dropped
#' @keywords internal
#' @noRd
warn_dropped <- function(..., n_boot) {
  dropped <- max(vapply(list(...), function(x) {
    if (is.null(x)) 0 else attr(x, "n_dropped") %||% 0
  }, numeric(1)))
  if (dropped > 0) {
    warning(sprintf(paste0(
      "%d of %d bootstrap replicates could not be evaluated (e.g. a compliance type ",
      "was absent from the resample) and were dropped."), dropped, n_boot), call. = FALSE)
  }
}

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Validate inference arguments
#' @keywords internal
#' @noRd
check_inference_args <- function(ci, n_boot, level) {
  if (!is.logical(ci) || length(ci) != 1L || is.na(ci)) {
    stop("`ci` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(n_boot) || length(n_boot) != 1L || is.na(n_boot) || n_boot < 2) {
    stop("`n_boot` must be a number >= 2.", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || is.na(level) || level <= 0 || level >= 1) {
    stop("`level` must be a number strictly between 0 and 1.", call. = FALSE)
  }
  invisible(TRUE)
}
