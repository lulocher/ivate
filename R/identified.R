# Identified quantities under the IV assumptions (1)-(5) of the manuscript.
#
# Internally, the identified quantities are stored as a named numeric vector
#   q = c(pi_C, pi_A, pi_N, mu_C0, mu_C1, mu_A1, mu_N0)
# so that bootstrap replicates can be stacked into a matrix. mu_A1 (mu_N0) is
# NA when there are no always-takers (never-takers).

q_names <- c("pi_C", "pi_A", "pi_N", "mu_C0", "mu_C1", "mu_A1", "mu_N0")

q_descriptions <- c(
  pi_C  = "Share of compliers",
  pi_A  = "Share of always-takers",
  pi_N  = "Share of never-takers",
  mu_C0 = "E[Y(0) | complier]",
  mu_C1 = "E[Y(1) | complier]",
  mu_A1 = "E[Y(1) | always-taker]",
  mu_N0 = "E[Y(0) | never-taker]",
  late  = "LATE = E[Y(1) - Y(0) | complier]"
)

#' Cell index for the four (instrument, treatment) cells
#'
#' Cells are ordered (z, a) = (0, 0), (0, 1), (1, 0), (1, 1).
#' @keywords internal
#' @noRd
iv_cells <- function(a, z) 1L + as.integer(a) + 2L * as.integer(z)

#' Sufficient statistics: cell counts and outcome sums
#' @keywords internal
#' @noRd
cell_stats <- function(y, g) {
  list(
    n = tabulate(g, nbins = 4L),
    s = as.vector(rowsum(c(y, 0, 0, 0, 0), c(g, 1:4), reorder = TRUE))
  )
}

#' Identified quantities from cell counts and outcome sums
#'
#' Implements the identification results in Appendix A of Locher, Richardson,
#' and Stensrud (2026):
#' \eqn{\pi_A = P(A = 1 \mid Z = 0)}, \eqn{\pi_N = P(A = 0 \mid Z = 1)},
#' \eqn{\mu_A^{(1)} = E[Y \mid A = 1, Z = 0]},
#' \eqn{\mu_N^{(0)} = E[Y \mid A = 0, Z = 1]}, and the complier means are
#' obtained by subtracting the always-/never-taker contributions from the
#' corresponding mixtures.
#'
#' @param n,s Length-4 vectors of cell counts and outcome sums, ordered as in
#'   `iv_cells()`.
#' @return A named numeric vector (see `q_names`).
#' @keywords internal
#' @noRd
identify_from_stats <- function(n, s) {
  pi_A <- n[2] / (n[1] + n[2])
  pi_N <- n[3] / (n[3] + n[4])
  pi_C <- 1 - pi_A - pi_N
  m    <- s / n

  mu_A1 <- if (isTRUE(pi_A > 0)) m[2] else NA_real_
  mu_N0 <- if (isTRUE(pi_N > 0)) m[3] else NA_real_

  at_part <- if (isTRUE(pi_A > 0)) pi_A * mu_A1 else 0
  nt_part <- if (isTRUE(pi_N > 0)) pi_N * mu_N0 else 0
  mu_C1   <- ((1 - pi_N) * m[4] - at_part) / pi_C
  mu_C0   <- ((1 - pi_A) * m[1] - nt_part) / pi_C

  stats::setNames(c(pi_C, pi_A, pi_N, mu_C0, mu_C1, mu_A1, mu_N0), q_names)
}

#' Identified quantities from data vectors
#' @keywords internal
#' @noRd
identify_quantities <- function(y, a, z) {
  cs <- cell_stats(y, iv_cells(a, z))
  identify_from_stats(cs$n, cs$s)
}

#' Identified quantities as a list, including the LATE
#'
#' This is the environment in which user-supplied restriction bounds (given as
#' one-sided formulas or functions) are evaluated.
#' @keywords internal
#' @noRd
q_as_list <- function(q) {
  c(as.list(q), list(late = unname(q["mu_C1"] - q["mu_C0"])))
}

#' Compliance setting implied by the identified type shares
#' @keywords internal
#' @noRd
iv_setting <- function(q) {
  at <- q[["pi_A"]] > 0
  nt <- q[["pi_N"]] > 0
  if (at && nt) "two-sided" else if (nt) "never-takers only" else if (at) "always-takers only" else "full compliance"
}

#' Sanity checks on the identified quantities
#'
#' Errors if the instrument is not relevant in the assumed direction and warns
#' if the plug-in estimates violate the IV model restrictions (identified
#' means outside the outcome range). This may indicate that the IV
#' assumptions fail, but can also arise from sampling variability.
#' @keywords internal
#' @noRd
check_identified <- function(q, y_range) {
  if (!(q[["pi_C"]] > 0)) {
    stop(sprintf(paste0(
      "The estimated share of compliers is %s <= 0. The instrument must increase ",
      "treatment uptake (P(A = 1 | Z = 1) > P(A = 1 | Z = 0)); consider recoding ",
      "the instrument as 1 - Z."), fmt_num(q[["pi_C"]])), call. = FALSE)
  }
  means <- q[c("mu_C0", "mu_C1", "mu_A1", "mu_N0")]
  means <- means[!is.na(means)]
  outside <- means < y_range[1] | means > y_range[2]
  if (any(outside)) {
    warning(sprintf(paste0(
      "The plug-in estimates violate the IV model restrictions: identified mean(s) ",
      "outside the outcome range [%s, %s]: %s. This may indicate that the IV ",
      "assumptions fail, but can also arise from sampling variability (e.g. a weak ",
      "instrument, a rare outcome, or a small sample)."),
      fmt_num(y_range[1]), fmt_num(y_range[2]),
      paste(names(means)[outside], collapse = ", ")), call. = FALSE)
  }
  invisible(q)
}

#' Violations of the IV model restrictions
#'
#' Under the IV assumptions, the observed-data distribution must imply a
#' positive share of compliers and complier means within the outcome range.
#' For binary outcomes these are the inequalities (11)-(13) of Richardson,
#' Evans, and Robins (2011) characterising the IV model; for other outcomes
#' they are necessary conditions. Plug-in estimates can violate them in
#' finite samples (e.g. with weak instruments, rare outcomes, or small
#' samples), yielding values incompatible with the model.
#'
#' @param qmat Matrix (one row per sample or draw) or vector of identified
#'   quantities.
#' @param y_range Outcome range.
#' @return A logical matrix with one column per restriction; `TRUE` means
#'   violated.
#' @keywords internal
#' @noRd
iv_restriction_violations <- function(qmat, y_range) {
  if (is.null(dim(qmat))) qmat <- matrix(qmat, 1, dimnames = list(NULL, names(qmat)))
  tol <- 1e-12
  outside <- function(m) is.na(m) | m < y_range[1] - tol | m > y_range[2] + tol
  cbind(
    pi_C  = !(is.finite(qmat[, "pi_C"]) & qmat[, "pi_C"] > 0),
    mu_C0 = outside(qmat[, "mu_C0"]),
    mu_C1 = outside(qmat[, "mu_C1"])
  )
}

#' Check the IV model restrictions in the sample and across bootstrap samples
#'
#' @return A list with the violations in the original sample (`sample`,
#'   logical) and, if bootstrap replicates are supplied, the share of
#'   replicates violating each restriction and any of them (`bootstrap`).
#' @keywords internal
#' @noRd
model_check <- function(q, boot = NULL, y_range, binary = FALSE, warn_share = 0.05) {
  out <- list(sample = iv_restriction_violations(q, y_range)[1, ], y_range = y_range)
  if (!is.null(boot)) {
    v <- iv_restriction_violations(boot, y_range)
    out$bootstrap <- c(colMeans(v), any = mean(apply(v, 1, any)))
    out$n_boot <- nrow(boot)
    if (out$bootstrap[["any"]] > warn_share) {
      warning(sprintf(paste0(
        "%s%% of bootstrap samples violate the IV model restrictions (pi_C > 0, complier ",
        "means within the outcome range); the bootstrap intervals may be unreliable.%s"),
        fmt_num(100 * out$bootstrap[["any"]]),
        if (binary) " Consider inference = \"bayes\"." else ""), call. = FALSE)
    }
  }
  out
}

