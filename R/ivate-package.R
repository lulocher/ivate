#' ivate: Informative inference on the ATE with instrumental variables
#'
#' Tools for sensitivity analysis and partial identification of the average
#' treatment effect (ATE) in instrumental variable (IV) settings with a binary
#' instrument and a binary treatment.
#'
#' Under the standard IV assumptions of Angrist, Imbens, and Rubin (1996),
#' point identification of the ATE is precluded by two unobserved quantities:
#' \eqn{\theta_A^{(0)}}, the mean potential outcome under control among
#' always-takers, and \eqn{\theta_N^{(1)}}, the mean potential outcome under
#' treatment among never-takers. The package implements the two
#' complementary strategies of Locher, Richardson, and Stensrud (2026):
#'
#' * [ate_sensitivity()] derives the thresholds on the sensitivity parameters
#'   at which the ATE crosses a policy-relevant cutoff and expresses them as
#'   implied comparisons with identified quantities.
#' * [ate_bounds()] computes bounds on the ATE under restrictions on the
#'   sensitivity parameters (see [restrictions]), nesting the Balke-Pearl
#'   bounds as a special case for binary outcomes.
#'
#' Both functions quantify uncertainty with the nonparametric bootstrap or,
#' for binary outcomes, by Bayesian inference (`inference = "bayes"`) based on
#' the transparent parametrization of Richardson, Evans, and Robins (2011).
#' [ate_bayes()] additionally gives the posterior of the ATE, either as a
#' function of the sensitivity parameters or under a prior on them (see
#' [priors]). Results can be visualised with
#' [plot()][plot.ivate_sensitivity()]; [expand_counts()] converts published
#' \eqn{2 \times 2 \times 2} tables into individual-level data.
#'
#' @references
#' Locher, L., Richardson, T., and Stensrud, M. (2026). Informative inference
#' on the average treatment effect with instrumental variables.
#'
#' Angrist, J. D., Imbens, G. W., and Rubin, D. B. (1996). Identification of
#' causal effects using instrumental variables. *Journal of the American
#' Statistical Association*, 91(434), 444--455.
#' \doi{10.1080/01621459.1996.10476902}
#'
#' Balke, A. and Pearl, J. (1997). Bounds on treatment effects from studies
#' with imperfect compliance. *Journal of the American Statistical
#' Association*, 92(439), 1171--1176. \doi{10.2307/2965583}
#'
#' Richardson, T. S., Evans, R. J., and Robins, J. M. (2011). Transparent
#' parametrizations of models for potential outcomes. In *Bayesian Statistics
#' 9*, 569--610. Oxford University Press.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom ggplot2 .data
## usethis namespace: end
NULL
