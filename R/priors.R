#' Priors on the sensitivity parameters
#'
#' Construct prior distributions for the unidentified sensitivity parameters
#' \eqn{\theta_A^{(0)}} and \eqn{\theta_N^{(1)}}, for use in [ate_bayes()]. The
#' prior is specified conditionally on the identified quantities: following
#' Richardson, Evans, and Robins (2011, eq. 27), the joint posterior factorises
#' as \eqn{p(\theta \mid q)\, p(q \mid \mathrm{data})}, so every posterior draw
#' of the identified quantities \eqn{q} is paired with a draw of \eqn{\theta}
#' from its prior given that \eqn{q}. Since the data carry no information
#' about \eqn{\theta}, the posterior of the ATE depends on this prior even
#' asymptotically; see Details of [ate_bayes()].
#'
#' * `prior_beta()`: a Beta distribution on the parameter itself.
#' * `prior_uniform()`: a uniform distribution on the feasible set implied by
#'   restrictions (see [restrictions]); without restrictions, uniform on
#'   \eqn{[0, 1]}.
#' * `prior_normal()`: a normal distribution on the chosen scale (as in
#'   [restrictions]: `"level"`, `"effect"`, `"ratio"`, or `"heterogeneity"`),
#'   truncated to parameter values in \eqn{[0, 1]}. For example,
#'   `prior_normal(mean = ~ late, sd = 0.05, on = "effect")` centres the
#'   treatment effect of the unidentified type at the complier effect.
#'
#' A single number passed instead of a prior is a point mass at that value.
#' The two sensitivity parameters are independent given the identified
#' quantities.
#'
#' @param shape1,shape2 Shape parameters of the Beta distribution.
#' @param ... Restrictions created by [restrict_level()], [restrict_effect()],
#'   [restrict_ratio()], or [restrict_heterogeneity()], or numeric intervals.
#' @param mean,sd Mean and standard deviation on the scale `on`. Numbers,
#'   one-sided formulas, or functions of the identified quantities (see
#'   [restrictions]).
#' @param on Scale on which the normal prior is placed.
#'
#' @return An object of class `ivate_prior`.
#'
#' @references
#' Richardson, T. S., Evans, R. J., and Robins, J. M. (2011). Transparent
#' parametrizations of models for potential outcomes. In *Bayesian Statistics
#' 9*, 569--610. Oxford University Press.
#'
#' @examples
#' prior_beta(2, 2)
#' prior_uniform(restrict_effect(lower = 0))
#' prior_normal(mean = ~ late, sd = 0.05, on = "effect")
#'
#' @name priors
NULL

#' @rdname priors
#' @export
prior_beta <- function(shape1 = 1, shape2 = 1) {
  for (s in list(shape1, shape2)) {
    if (!is.numeric(s) || length(s) != 1L || is.na(s) || s <= 0) {
      stop("Beta shape parameters must be single positive numbers.", call. = FALSE)
    }
  }
  structure(list(type = "beta", shape1 = shape1, shape2 = shape2), class = "ivate_prior")
}

#' @rdname priors
#' @export
prior_uniform <- function(...) {
  restr <- as_restriction_list(list(...), "prior_uniform")
  structure(list(type = "uniform", restrictions = restr), class = "ivate_prior")
}

#' @rdname priors
#' @export
prior_normal <- function(mean, sd, on = c("level", "effect", "ratio", "heterogeneity")) {
  on <- match.arg(on)
  for (arg in c("mean", "sd")) {
    v <- get(arg)
    ok <- (is.numeric(v) && length(v) == 1L && !is.na(v)) || inherits(v, "formula") || is.function(v)
    if (!ok) stop(sprintf("`%s` must be a single number, a one-sided formula, or a function.", arg),
                  call. = FALSE)
  }
  if (is.numeric(sd) && sd <= 0) stop("`sd` must be positive.", call. = FALSE)
  structure(list(type = "normal", mean = mean, sd = sd, on = on), class = "ivate_prior")
}

#' Point-mass prior (internal; created from numbers)
#' @keywords internal
#' @noRd
prior_point <- function(value) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || value < 0 || value > 1) {
    stop("A fixed sensitivity parameter must be a single number in [0, 1].", call. = FALSE)
  }
  structure(list(type = "point", value = value), class = "ivate_prior")
}

#' Map a value on a scale to the sensitivity parameter
#' @keywords internal
#' @noRd
scale_to_theta <- function(v, on, parameter, q) {
  is_N <- parameter == "theta_N1"
  switch(on,
    level         = v,
    effect        = if (is_N) q[["mu_N0"]] + v else q[["mu_A1"]] - v,
    ratio         = if (is_N) v * q[["mu_N0"]] else q[["mu_A1"]] / v,
    heterogeneity = if (is_N) q[["mu_C1"]] + v else q[["mu_C0"]] + v
  )
}

#' Range of values on a scale corresponding to a parameter in the unit interval
#' @keywords internal
#' @noRd
scale_range <- function(on, parameter, q) {
  is_N <- parameter == "theta_N1"
  switch(on,
    level         = c(0, 1),
    effect        = if (is_N) c(0, 1) - q[["mu_N0"]] else q[["mu_A1"]] - c(1, 0),
    heterogeneity = if (is_N) c(0, 1) - q[["mu_C1"]] else c(0, 1) - q[["mu_C0"]],
    ratio = {
      ref <- if (is_N) q[["mu_N0"]] else q[["mu_A1"]]
      if (!(ref > 0)) stop("Ratio priors require a positive reference mean.", call. = FALSE)
      if (is_N) c(0, 1 / ref) else c(ref, Inf)
    }
  )
}

#' Draw from a truncated normal distribution by inversion
#' @keywords internal
#' @noRd
rtruncnorm1 <- function(mean, sd, lower, upper) {
  p_lo <- stats::pnorm(lower, mean, sd)
  p_hi <- stats::pnorm(upper, mean, sd)
  if (p_hi - p_lo < 1e-12) return(if (mean < lower) lower else upper)
  min(max(stats::qnorm(stats::runif(1, p_lo, p_hi), mean, sd), lower), upper)
}

#' Draw the sensitivity parameter from its prior given each draw of q
#'
#' @param prior An `ivate_prior`.
#' @param parameter `"theta_A0"` or `"theta_N1"`.
#' @param qmat Matrix of draws of the identified quantities.
#' @return A numeric vector with one draw per row of `qmat` (`NA` where the
#'   prior could not be evaluated, e.g. an empty feasible set).
#' @keywords internal
#' @noRd
sample_prior <- function(prior, parameter, qmat) {
  n <- nrow(qmat)
  switch(prior$type,
    point = rep(prior$value, n),
    beta  = stats::rbeta(n, prior$shape1, prior$shape2),
    uniform = {
      if (length(prior$restrictions) == 0) return(stats::runif(n))
      vapply(seq_len(n), function(s) {
        tryCatch({
          iv <- feasible_interval(prior$restrictions, parameter, qmat[s, ], c(0, 1))
          stats::runif(1, iv[1], iv[2])
        }, error = function(e) NA_real_)
      }, numeric(1))
    },
    normal = vapply(seq_len(n), function(s) {
      tryCatch({
        q  <- qmat[s, ]
        m  <- eval_limit(prior$mean, q)
        sd <- eval_limit(prior$sd, q)
        rg <- scale_range(prior$on, parameter, q)
        v  <- rtruncnorm1(m, sd, rg[1], rg[2])
        min(max(scale_to_theta(v, prior$on, parameter, q), 0), 1)
      }, error = function(e) NA_real_)
    }, numeric(1))
  )
}

#' Normalise a user specification into an ivate_prior
#' @keywords internal
#' @noRd
as_prior <- function(spec, arg) {
  if (inherits(spec, "ivate_prior")) return(spec)
  if (is.numeric(spec) && length(spec) == 1L) return(prior_point(spec))
  if (inherits(spec, "ivate_restriction")) {
    stop(sprintf("Wrap restrictions in prior_uniform() to use them as a prior for `%s`.", arg),
         call. = FALSE)
  }
  stop(sprintf("`%s` must be a prior (see ?priors) or a single number.", arg), call. = FALSE)
}

#' Describe a prior in words
#' @keywords internal
#' @noRd
describe_prior <- function(prior, parameter) {
  switch(prior$type,
    point   = paste(parameter, "=", fmt_num(prior$value)),
    beta    = sprintf("%s ~ Beta(%s, %s)", parameter, fmt_num(prior$shape1), fmt_num(prior$shape2)),
    uniform = if (length(prior$restrictions) == 0) {
      sprintf("%s ~ Uniform(0, 1)", parameter)
    } else {
      sprintf("%s ~ Uniform on {%s}", parameter, describe_restrictions(prior$restrictions, parameter))
    },
    normal  = sprintf("%s ~ Normal(%s, %s), truncated to %s in [0, 1]",
                      scale_expression(prior$on, parameter), fmt_limit(prior$mean),
                      fmt_limit(prior$sd), parameter)
  )
}

#' @export
print.ivate_prior <- function(x, ...) {
  cat("<ivate prior>\n")
  cat("  N:", describe_prior(x, "theta_N1"), "\n")
  cat("  A:", describe_prior(x, "theta_A0"), "\n")
  invisible(x)
}
