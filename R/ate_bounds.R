#' Bounds on the average treatment effect under IV assumptions
#'
#' Computes bounds on the average treatment effect (ATE) by plugging in the
#' most extreme feasible values for the unidentified sensitivity parameters
#' \eqn{\theta_A = \mu_A^{(0)}} (mean potential outcome under control for
#' always-takers) and \eqn{\theta_N = \mu_N^{(1)}} (mean potential outcome
#' under treatment for never-takers). When the outcome is binary, these
#' coincide with the sharp Balke-Pearl bounds (see Appendix B of Locher and
#' Stensrud (2025)).
#'
#' The ATE decomposes as (equation (8) in Locher and Stensrud (2025)):
#' \deqn{
#'   \Psi_{ATE} = \pi_C(\mu_C^{(1)} - \mu_C^{(0)})
#'              + \pi_A(\mu_A^{(1)} - \theta_A)
#'              + \pi_N(\theta_N   - \mu_N^{(0)}).
#' }
#' Since \eqn{\Psi_{ATE}} is decreasing in \eqn{\theta_A} and increasing in
#' \eqn{\theta_N}, the bounds are obtained by evaluating at the corners of the
#' assumed feasible set (equation (12)):
#' \deqn{
#'   \Psi_{ATE} \in
#'     [\Psi_{ATE}(u_A, \ell_N),\; \Psi_{ATE}(\ell_A, u_N)].
#' }
#' Restricting \eqn{\theta_A} and \eqn{\theta_N} to a credible sub-interval
#' of their feasible range narrows the bounds. If the resulting interval
#' excludes the decision threshold \eqn{x}, the policy decision can be
#' resolved under the imposed assumptions.
#'
#' @param data A data frame containing the outcome, treatment, and instrument
#'   variables.
#' @param formula A two-part formula of the form
#'   \code{outcome ~ treatment | instrument}.
#' @param bounds_y0_always_takers Numeric vector of length 2 giving the
#'   assumed lower and upper bounds \eqn{[\ell_A, u_A]} for
#'   \eqn{\theta_A = \mu_A^{(0)}}. Defaults to \code{c(0, 1)}.
#' @param bounds_y1_never_takers Numeric vector of length 2 giving the
#'   assumed lower and upper bounds \eqn{[\ell_N, u_N]} for
#'   \eqn{\theta_N = \mu_N^{(1)}}. Defaults to \code{c(0, 1)}.
#'
#' @return A named numeric vector with elements \code{"Lower bound for ATE"}
#'   and \code{"Upper bound for ATE"}.
#'
#' @references
#' Angrist, J. D., Imbens, G. W., and Rubin, D. B. (1996). Identification of
#' causal effects using instrumental variables. \emph{Journal of the American
#' Statistical Association}, 91(434), 444--455.
#' \doi{10.1080/01621459.1996.10476902}
#'
#' Balke, A. and Pearl, J. (1997). Bounds on treatment effects from studies
#' with imperfect compliance. \emph{Journal of the American Statistical
#' Association}, 92(439), 1171--1176.
#' \doi{10.2307/2965583}
#'
#' Locher, L. and Stensrud, M. J. (2025). Informative inference on the average
#' treatment effect with instrumental variables.
#'
#' @examples
#' # Simulate a dataset with two-sided non-compliance
#' set.seed(42)
#' n   <- 2000
#' z   <- rbinom(n, 1, 0.5)
#' a   <- ifelse(z == 1, rbinom(n, 1, 0.7), rbinom(n, 1, 0.2))
#' y   <- rbinom(n, 1, 0.3 + 0.2 * a)
#' dat <- data.frame(y = y, a = a, z = z)
#'
#' ate_bounds(dat, y ~ a | z)
#'
#' # Tighten the bounds by restricting the sensitivity parameters
#' ate_bounds(dat, y ~ a | z,
#'            bounds_y0_always_takers = c(0.1, 0.5),
#'            bounds_y1_never_takers  = c(0.2, 0.6))
#'
#' @export
ate_bounds <- function(data, formula,
                       bounds_y0_always_takers = c(0, 1),
                       bounds_y1_never_takers  = c(0, 1)) {
  d  <- validate_iv_data(data, formula)
  kc <- get_known_quantities(d$outcome, d$treatment, d$instrument)

  ate_lower <- kc$p_c * (kc$y1_c - kc$y0_c) +
    kc$p_at * (kc$y1_at - bounds_y0_always_takers[2]) +
    kc$p_nt * (bounds_y1_never_takers[1] - kc$y0_nt)

  ate_upper <- kc$p_c * (kc$y1_c - kc$y0_c) +
    kc$p_at * (kc$y1_at - bounds_y0_always_takers[1]) +
    kc$p_nt * (bounds_y1_never_takers[2] - kc$y0_nt)

  c("Lower bound for ATE" = ate_lower, "Upper bound for ATE" = ate_upper)
}
