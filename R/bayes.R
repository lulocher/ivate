# Bayesian inference for the identified quantities, following Section 3 of
# Richardson, Evans, and Robins (2011). The identified part of the model is
# parametrised directly by the observed distributions p(y, a | z = 0) and
# p(y, a | z = 1), each given a Dirichlet prior. Under the IV assumptions
# (randomisation, exclusion, monotonicity) these distributions must satisfy
# the inequalities (11)-(13) of RER, which are equivalent to
#   pi_C >= 0, 0 <= mu_C0 <= 1, 0 <= mu_C1 <= 1.
# The posterior (and prior) restricted to this set is sampled by conjugacy and
# rejection. Only binary outcomes are supported.

# Cell order within an instrument arm: (a, y) = (0, 0), (0, 1), (1, 0), (1, 1).
arm_cells <- c("a0y0", "a0y1", "a1y0", "a1y1")

#' Stop unless the outcome is binary
#' @keywords internal
#' @noRd
check_binary_outcome <- function(y, y_range = NULL) {
  if (!is_binary(y)) {
    stop(paste0("Bayesian inference (Richardson, Evans, and Robins, 2011) requires a binary ",
                "(0/1) outcome. Use inference = \"bootstrap\" for other outcomes."), call. = FALSE)
  }
  if (!is.null(y_range) && !isTRUE(all.equal(as.numeric(y_range), c(0, 1)))) {
    stop("With Bayesian inference the outcome range is fixed at c(0, 1).", call. = FALSE)
  }
  invisible(TRUE)
}

#' Cell counts per instrument arm
#' @return A list with length-4 count vectors `z0` and `z1` (order `arm_cells`).
#' @keywords internal
#' @noRd
arm_counts <- function(y, a, z) {
  idx <- 1L + as.integer(y) + 2L * as.integer(a)
  list(z0 = tabulate(idx[z == 0], nbins = 4L), z1 = tabulate(idx[z == 1], nbins = 4L))
}

#' Structural zeros implied by the compliance setting
#'
#' Without always-takers, P(A = 1 | Z = 0) = 0; without never-takers,
#' P(A = 0 | Z = 1) = 0. These cells receive no prior mass.
#' @keywords internal
#' @noRd
structural_cells <- function(setting) {
  list(
    z0 = if (setting %in% c("never-takers only", "full compliance")) c(FALSE, FALSE, TRUE, TRUE) else rep(FALSE, 4),
    z1 = if (setting %in% c("always-takers only", "full compliance")) c(TRUE, TRUE, FALSE, FALSE) else rep(FALSE, 4)
  )
}

#' Draw from Dirichlet distributions with structural zeros
#' @return An n x 4 matrix of cell probabilities.
#' @keywords internal
#' @noRd
rdirichlet_cells <- function(n, alpha, structural) {
  out <- matrix(0, n, 4L)
  free <- which(!structural)
  g <- vapply(free, function(k) stats::rgamma(n, shape = alpha[k]), numeric(n))
  g <- matrix(g, nrow = n)
  out[, free] <- g / rowSums(g)
  out
}

#' Identified quantities from cell probabilities (vectorised over draws)
#' @param P0,P1 n x 4 matrices of p(a, y | z = 0) and p(a, y | z = 1).
#' @return An n x 7 matrix with columns `q_names`.
#' @keywords internal
#' @noRd
q_from_probs <- function(P0, P1) {
  pi_A <- P0[, 3] + P0[, 4]
  pi_N <- P1[, 1] + P1[, 2]
  pi_C <- 1 - pi_A - pi_N
  mu_A1 <- ifelse(pi_A > 0, P0[, 4] / pi_A, NA_real_)
  mu_N0 <- ifelse(pi_N > 0, P1[, 2] / pi_N, NA_real_)
  pimu_C1 <- P1[, 4] - P0[, 4]
  pimu_C0 <- P0[, 2] - P1[, 2]
  mu_C1 <- ifelse(pi_C > 0, pimu_C1 / pi_C, NA_real_)
  mu_C0 <- ifelse(pi_C > 0, pimu_C0 / pi_C, NA_real_)
  out <- cbind(pi_C, pi_A, pi_N, mu_C0, mu_C1, mu_A1, mu_N0, pimu_C0, pimu_C1)
  colnames(out) <- q_names
  out
}

#' Which draws violate the IV inequalities (11)-(13) of RER (2011)
#' @return A logical matrix with columns `(11)`, `(12)`, `(13)`.
#' @keywords internal
#' @noRd
iv_violations <- function(P0, P1) {
  tol <- 1e-12
  cbind(
    "(11)" = (P1[, 3] + P1[, 4]) < (P0[, 3] + P0[, 4]) - tol,
    "(12)" = P1[, 2] > P0[, 2] + tol | P1[, 1] > P0[, 1] + tol,
    "(13)" = P1[, 4] < P0[, 4] - tol | P1[, 3] < P0[, 3] - tol
  )
}

#' Sample identified quantities from the constrained Dirichlet posterior
#'
#' @param counts Output of `arm_counts()`, or `NULL` to sample from the prior.
#' @param setting Compliance setting (determines structural zeros).
#' @param n_draws Number of accepted draws.
#' @param prior_obs Dirichlet concentration parameter (a single positive
#'   number, applied to every non-structural cell).
#' @return A list with the accepted draws `q` (matrix), the probabilities of
#'   violating each inequality under the unconstrained Dirichlet (`violation`),
#'   and the number of proposals (`n_proposed`).
#' @keywords internal
#' @noRd
sample_dirichlet_q <- function(counts, setting, n_draws, prior_obs, max_proposals = 1000 * n_draws) {
  structural <- structural_cells(setting)
  if (is.null(counts)) counts <- list(z0 = rep(0, 4), z1 = rep(0, 4))
  alpha0 <- counts$z0 + prior_obs
  alpha1 <- counts$z1 + prior_obs
  if (any(counts$z0[structural$z0] > 0) || any(counts$z1[structural$z1] > 0)) {
    stop("Observations found in structurally empty cells.", call. = FALSE)
  }

  accepted <- matrix(NA_real_, 0, length(q_names))
  n_viol <- c("(11)" = 0, "(12)" = 0, "(13)" = 0, any = 0)
  n_proposed <- 0
  batch <- max(1000L, 2L * n_draws)
  while (nrow(accepted) < n_draws && n_proposed < max_proposals) {
    P0 <- rdirichlet_cells(batch, alpha0, structural$z0)
    P1 <- rdirichlet_cells(batch, alpha1, structural$z1)
    viol <- iv_violations(P0, P1)
    bad <- apply(viol, 1, any)
    n_viol <- n_viol + c(colSums(viol), any = sum(bad))
    n_proposed <- n_proposed + batch
    if (any(!bad)) accepted <- rbind(accepted, q_from_probs(P0[!bad, , drop = FALSE], P1[!bad, , drop = FALSE]))
  }
  if (nrow(accepted) < n_draws) {
    stop(sprintf(paste0(
      "Only %d of %d proposals satisfied the IV inequalities: the data are (nearly) ",
      "incompatible with the IV assumptions (posterior probability of a violation %s)."),
      nrow(accepted), n_proposed, fmt_num(n_viol[["any"]] / n_proposed)), call. = FALSE)
  }
  list(q = accepted[seq_len(n_draws), , drop = FALSE],
       violation = n_viol / n_proposed,
       n_proposed = n_proposed)
}

#' Posterior draws of the identified quantities with diagnostics
#' @keywords internal
#' @noRd
posterior_q <- function(y, a, z, setting, n_draws, prior_obs) {
  check_binary_outcome(y)
  post <- sample_dirichlet_q(arm_counts(y, a, z), setting, n_draws, prior_obs)
  if (post$violation[["any"]] > 0.5) {
    warning(sprintf(paste0(
      "Posterior probability that the IV inequalities are violated is %s (under the ",
      "unconstrained Dirichlet posterior): the data cast doubt on the IV assumptions."),
      fmt_num(post$violation[["any"]])), call. = FALSE)
  }
  post
}

#' Validate Bayesian inference arguments
#' @keywords internal
#' @noRd
check_bayes_args <- function(n_draws, prior_obs) {
  if (!is.numeric(n_draws) || length(n_draws) != 1L || is.na(n_draws) || n_draws < 2) {
    stop("`n_draws` must be a number >= 2.", call. = FALSE)
  }
  if (!is.numeric(prior_obs) || length(prior_obs) != 1L || is.na(prior_obs) || prior_obs <= 0) {
    stop("`prior_obs` must be a single positive number (the Dirichlet concentration).",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Draws of the identified quantities for the chosen inference method
#'
#' @return `NULL` if `ci = FALSE`, otherwise a list with the draw matrix
#'   `draws` and, for Bayesian inference, the `bayes` diagnostics.
#' @keywords internal
#' @noRd
uncertainty_draws <- function(d, setting, ci, inference, n_boot, boot_strata, n_draws, prior_obs) {
  if (inference == "bayes") {
    post <- posterior_q(d$y, d$a, d$z, setting, n_draws, prior_obs)
    return(list(draws = post$q, bayes = list(
      prior_obs = prior_obs, n_draws = n_draws, violation = post$violation,
      n_proposed = post$n_proposed)))
  }
  if (!ci) return(NULL)
  list(draws = bootstrap_q(d$y, d$a, d$z, n_boot, boot_strata), bayes = NULL)
}
