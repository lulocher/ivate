#' Bounds on the ATE under restrictions on the sensitivity parameters
#'
#' Propagates restrictions on the unidentified quantities
#' \eqn{\theta_A^{(0)} = E[Y(0) \mid A]} and
#' \eqn{\theta_N^{(1)} = E[Y(1) \mid N]} into bounds on the average treatment
#' effect (Section 2.5 of Locher, Richardson, and Stensrud, 2026).
#'
#' Each restriction translates into an interval for the corresponding
#' sensitivity parameter; restrictions can be stated directly on the
#' parameter or on derived quantities such as the treatment effect, the
#' ratio, or the outcome heterogeneity relative to compliers (see
#' [restrictions]). The feasible set is the intersection of these intervals
#' with the outcome range. Since the ATE is decreasing in
#' \eqn{\theta_A^{(0)}} and increasing in \eqn{\theta_N^{(1)}}, the bounds are
#' obtained at the corners of the feasible set (equation (11)):
#' \deqn{\Psi_{ATE} \in [\Psi_{ATE}(u_A, \ell_N),\ \Psi_{ATE}(\ell_A, u_N)].}
#'
#' Without restrictions and with a binary outcome, the bounds coincide with
#' the Balke-Pearl bounds under monotonicity (Appendix B of the manuscript);
#' these unrestricted bounds are always reported for comparison.
#'
#' Confidence intervals for each bound are percentile intervals from the
#' nonparametric bootstrap (Section 2.6); restrictions that depend on
#' identified quantities (e.g. `restrict_effect(lower = ~ late)`) are
#' re-evaluated in every replicate. With `inference = "bayes"` (binary
#' outcomes only), bounds are summarised by posterior medians and credible
#' intervals under a Dirichlet prior on the observed distribution restricted
#' to the IV model (Richardson, Evans, and Robins, 2011), who recommend such
#' posteriors on bounds for partially identified parameters. Call
#' [set.seed()] beforehand for reproducible intervals.
#'
#' @inheritParams ate_sensitivity
#' @param theta_A0,theta_N1 Restrictions on the sensitivity parameters:
#'   `NULL` (unrestricted, the default), a single number (fixes the
#'   parameter), a length-2 interval, a restriction created by
#'   [restrict_level()], [restrict_effect()], [restrict_ratio()], or
#'   [restrict_heterogeneity()], or a list of these (intersected). Ignored
#'   (with a warning) if the corresponding compliance type is absent.
#'
#' @return An object of class `ivate_bounds` with elements
#'   \describe{
#'     \item{`bounds`}{Data frame with the lower and upper bound on the ATE
#'       and their confidence intervals.}
#'     \item{`unrestricted`}{The same for the bounds without restrictions.}
#'     \item{`feasible`}{Data frame with the feasible interval for each
#'       sensitivity parameter and a description of the restrictions.}
#'     \item{`identified`}{Data frame of identified quantities.}
#'     \item{`model_check`}{Violations of the IV model restrictions in the
#'       original sample and the share of bootstrap samples violating them;
#'       see the section "Model check" in [ate_sensitivity()].}
#'   }
#'   together with the inputs.
#'
#' @references
#' Locher, L., Richardson, T., and Stensrud, M. (2026). Informative inference
#' on the average treatment effect with instrumental variables.
#'
#' Balke, A. and Pearl, J. (1997). Bounds on treatment effects from studies
#' with imperfect compliance. *Journal of the American Statistical
#' Association*, 92(439), 1171--1176.
#'
#' Richardson, T. S., Evans, R. J., and Robins, J. M. (2011). Transparent
#' parametrizations of models for potential outcomes. In *Bayesian Statistics
#' 9*, 569--610. Oxford University Press.
#'
#' @seealso [ate_sensitivity()], [restrictions].
#'
#' @examples
#' set.seed(1)
#' n <- 2000
#' z <- rbinom(n, 1, 0.5)
#' a <- ifelse(z == 1, rbinom(n, 1, 0.7), rbinom(n, 1, 0.2))
#' y <- rbinom(n, 1, 0.3 + 0.2 * a)
#' dat <- data.frame(y, a, z)
#'
#' # Unrestricted (Balke-Pearl) bounds
#' ate_bounds(dat, y ~ a | z, n_boot = 200)
#'
#' # Effect among always-takers >= LATE >= effect among never-takers
#' ate_bounds(dat, y ~ a | z,
#'            theta_A0 = restrict_effect(lower = ~ late),
#'            theta_N1 = restrict_effect(upper = ~ late),
#'            ci = FALSE)
#'
#' # Treatment at most doubles the outcome among never-takers
#' ate_bounds(dat, y ~ a | z, theta_N1 = restrict_ratio(upper = 2), ci = FALSE)
#'
#' @export
ate_bounds <- function(data, formula,
                       theta_A0    = NULL,
                       theta_N1    = NULL,
                       y_range     = NULL,
                       ci          = TRUE,
                       n_boot      = 1000,
                       level       = 0.95,
                       boot_strata = c("none", "instrument"),
                       inference   = c("bootstrap", "bayes"),
                       n_draws     = 4000,
                       prior_obs   = 1) {
  boot_strata <- match.arg(boot_strata)
  inference   <- match.arg(inference)
  check_inference_args(ci, n_boot, level)
  if (inference == "bayes") check_bayes_args(n_draws, prior_obs)

  d       <- prepare_iv_data(data, formula)
  if (inference == "bayes") check_binary_outcome(d$y, y_range)
  y_range <- resolve_y_range(d$y, y_range)
  q       <- identify_quantities(d$y, d$a, d$z)
  check_identified(q, y_range)
  setting <- iv_setting(q)

  restr_A <- as_restriction_list(theta_A0, "theta_A0")
  restr_N <- as_restriction_list(theta_N1, "theta_N1")
  if (length(restr_A) && q[["pi_A"]] == 0) {
    warning("There are no always-takers; `theta_A0` is ignored.", call. = FALSE)
    restr_A <- list()
  }
  if (length(restr_N) && q[["pi_N"]] == 0) {
    warning("There are no never-takers; `theta_N1` is ignored.", call. = FALSE)
    restr_N <- list()
  }

  unc    <- uncertainty_draws(d, setting, ci, inference, n_boot, boot_strata, n_draws, prior_obs)
  draws  <- unc$draws
  mcheck <- model_check(q, if (inference == "bootstrap") draws else NULL, y_range, is_binary(d$y))
  center <- if (inference == "bayes") "median" else "plugin"

  restricted   <- evaluate_stat(bounds_stat_fun(restr_A, restr_N, y_range), q, draws, level, center)
  unrestricted <- evaluate_stat(bounds_stat_fun(list(), list(), y_range), q, draws, level, center)
  if (!is.null(draws)) warn_dropped(restricted, unrestricted, n_draws = nrow(draws))

  feasible <- data.frame(
    parameter   = c("theta_A0", "theta_N1"),
    lower       = restricted$estimate[restricted$term %in% c("theta_A0_lower", "theta_N1_lower")],
    upper       = restricted$estimate[restricted$term %in% c("theta_A0_upper", "theta_N1_upper")],
    restriction = c(describe_restrictions(restr_A, "theta_A0"),
                    describe_restrictions(restr_N, "theta_N1")),
    stringsAsFactors = FALSE
  )
  feasible <- feasible[c(q[["pi_A"]] > 0, q[["pi_N"]] > 0), , drop = FALSE]
  rownames(feasible) <- NULL

  structure(list(
    call         = match.call(),
    vars         = d$vars,
    n            = length(d$y),
    setting      = setting,
    y_range      = y_range,
    q            = q,
    identified   = identified_table(q, draws, level, center),
    bounds       = bounds_rows(restricted),
    unrestricted = bounds_rows(unrestricted),
    feasible     = feasible,
    draws        = draws,
    inference    = inference,
    bayes        = unc$bayes,
    model_check  = mcheck,
    level        = level
  ), class = "ivate_bounds")
}

#' Statistic function for the bounds and the feasible intervals
#' @keywords internal
#' @noRd
bounds_stat_fun <- function(restr_A, restr_N, y_range) {
  function(q) {
    iA <- if (isTRUE(q[["pi_A"]] > 0)) feasible_interval(restr_A, "theta_A0", q, y_range) else c(NA_real_, NA_real_)
    iN <- if (isTRUE(q[["pi_N"]] > 0)) feasible_interval(restr_N, "theta_N1", q, y_range) else c(NA_real_, NA_real_)
    c(lower          = ate_at(q, iA[2], iN[1]),
      upper          = ate_at(q, iA[1], iN[2]),
      theta_A0_lower = iA[1], theta_A0_upper = iA[2],
      theta_N1_lower = iN[1], theta_N1_upper = iN[2])
  }
}

#' Extract the bound rows from an evaluate_stat() table
#' @keywords internal
#' @noRd
bounds_rows <- function(tab) {
  out <- tab[tab$term %in% c("lower", "upper"), , drop = FALSE]
  names(out)[1] <- "bound"
  rownames(out) <- NULL
  out
}

#' Describe a list of restrictions
#' @keywords internal
#' @noRd
describe_restrictions <- function(restr, parameter) {
  if (length(restr) == 0) return("none (outcome range)")
  paste(vapply(restr, describe_restriction, character(1), parameter = parameter), collapse = "; ")
}
