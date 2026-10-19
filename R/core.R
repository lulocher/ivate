# Closed-form expressions from Locher, Richardson, and Stensrud (2026). All
# functions take the identified quantities `q` (a named numeric
# vector, see identified.R) so that they can be re-evaluated on every
# bootstrap replicate.

# Contribution of the compliers to the ATE, pi_C * LATE, written in terms of
# the complier contributions so that it is defined when pi_C = 0 (it then
# equals the intention-to-treat effect on the outcome).
complier_term <- function(q) q[["pimu_C1"]] - q[["pimu_C0"]]

# Contribution of a type with share `p` and effect `d`; zero when the type is
# absent (its potential outcome means are then undefined).
type_term <- function(p, d) if (isTRUE(p > 0)) p * d else rep(0, length(d))

#' ATE as a function of the sensitivity parameters (the ATE decomposition)
#' @keywords internal
#' @noRd
ate_at <- function(q, theta_A0, theta_N1) {
  unname(
    complier_term(q) +
      type_term(q[["pi_A"]], q[["mu_A1"]] - theta_A0) +
      type_term(q[["pi_N"]], theta_N1 - q[["mu_N0"]])
  )
}

#' Threshold on theta_N1 at which ATE = cutoff (two- and one-sided case)
#'
#' `theta_A0` is ignored when there are no always-takers.
#' @keywords internal
#' @noRd
threshold_N1 <- function(q, cutoff, theta_A0 = NA_real_) {
  unname(
    (cutoff - complier_term(q) -
       type_term(q[["pi_A"]], q[["mu_A1"]] - theta_A0) +
       q[["pi_N"]] * q[["mu_N0"]]) / q[["pi_N"]]
  )
}

#' Threshold on theta_A0 at which ATE = cutoff
#'
#' `theta_N1` is ignored when there are no never-takers.
#' @keywords internal
#' @noRd
threshold_A0 <- function(q, cutoff, theta_N1 = NA_real_) {
  unname(
    (complier_term(q) +
       q[["pi_A"]] * q[["mu_A1"]] +
       type_term(q[["pi_N"]], theta_N1 - q[["mu_N0"]]) - cutoff) / q[["pi_A"]]
  )
}

#' Coefficients of the boundary B (two-sided non-compliance)
#'
#' The boundary \eqn{\Psi_{ATE}(\theta_A^{(0)}, \theta_N^{(1)}) = x} is a line,
#' expressed here in four equivalent ways: solved for either sensitivity
#' parameter, on the level scale and on the scale of the implied type-specific
#' effects \eqn{\delta_A = \mu_A^{(1)} - \theta_A^{(0)}} and
#' \eqn{\delta_N = \theta_N^{(1)} - \mu_N^{(0)}}.
#' @keywords internal
#' @noRd
boundary_coefficients <- function(q, cutoff) {
  cl <- complier_term(q)
  pA <- q[["pi_A"]]; pN <- q[["pi_N"]]
  c(
    level_N1_intercept  = (cutoff - cl - pA * q[["mu_A1"]] + pN * q[["mu_N0"]]) / pN,
    level_N1_slope      = pA / pN,
    level_A0_intercept  = (cl + pA * q[["mu_A1"]] - pN * q[["mu_N0"]] - cutoff) / pA,
    level_A0_slope      = pN / pA,
    effect_N_intercept  = (cutoff - cl) / pN,
    effect_N_slope      = -pA / pN,
    effect_A_intercept  = (cutoff - cl) / pA,
    effect_A_slope      = -pN / pA
  )
}

#' Least favourable value of a sensitivity parameter within an interval
#'
#' The ATE is decreasing in theta_A0 and increasing in theta_N1. For the
#' criterion ATE > x the least favourable value is therefore the upper end for
#' theta_A0 and the lower end for theta_N1; for ATE < x it is the reverse.
#' @keywords internal
#' @noRd
least_favourable <- function(interval, parameter, direction) {
  upper_end <- xor(parameter == "theta_A0", direction == "less")
  if (upper_end) interval[2] else interval[1]
}

#' Implied comparisons at a threshold value (as in the comparison figure of the manuscript)
#'
#' For theta_N1 = tau: the implied effect among never-takers, the implied
#' ratio, and the implied heterogeneity between never-takers and compliers
#' under treatment, each with an identified benchmark. For theta_A0 = tau the
#' always-taker analogues are returned.
#' @keywords internal
#' @noRd
implied_comparisons <- function(q, tau, parameter) {
  late <- q[["mu_C1"]] - q[["mu_C0"]]
  ratio_C <- q[["mu_C1"]] / q[["mu_C0"]]
  if (parameter == "theta_N1") {
    c(effect_implied          = tau - q[["mu_N0"]],
      effect_benchmark        = late,
      ratio_implied           = tau / q[["mu_N0"]],
      ratio_benchmark         = ratio_C,
      heterogeneity_implied   = tau - q[["mu_C1"]],
      heterogeneity_benchmark = q[["mu_N0"]] - q[["mu_C0"]])
  } else {
    c(effect_implied          = q[["mu_A1"]] - tau,
      effect_benchmark        = late,
      ratio_implied           = q[["mu_A1"]] / tau,
      ratio_benchmark         = ratio_C,
      heterogeneity_implied   = tau - q[["mu_C0"]],
      heterogeneity_benchmark = q[["mu_A1"]] - q[["mu_C1"]])
  }
}

# Labels for the implied-comparison table.
implied_labels <- list(
  theta_N1 = data.frame(
    comparison = c("effect", "ratio", "heterogeneity"),
    implied    = c("theta_N1 - mu_N0", "theta_N1 / mu_N0", "theta_N1 - mu_C1"),
    benchmark  = c("mu_C1 - mu_C0", "mu_C1 / mu_C0", "mu_N0 - mu_C0"),
    stringsAsFactors = FALSE
  ),
  theta_A0 = data.frame(
    comparison = c("effect", "ratio", "heterogeneity"),
    implied    = c("mu_A1 - theta_A0", "mu_A1 / theta_A0", "theta_A0 - mu_C0"),
    benchmark  = c("mu_C1 - mu_C0", "mu_C1 / mu_C0", "mu_A1 - mu_C1"),
    stringsAsFactors = FALSE
  )
)
