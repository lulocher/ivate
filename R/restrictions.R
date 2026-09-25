#' Restrictions on the sensitivity parameters
#'
#' Construct restrictions on a sensitivity parameter, stated on whichever
#' scale is most natural (Section 2.5 of Locher, Richardson, and Stensrud,
#' 2026). A restriction is passed to the `theta_A0` or `theta_N1` argument of
#' [ate_bounds()] or [ate_sensitivity()], which determines the parameter it
#' applies to. Each restriction translates into an interval for the
#' sensitivity parameter; several restrictions (in a list) are intersected,
#' and the result is always intersected with the outcome range.
#'
#' The scales are defined within compliance type (effect, ratio) or relative
#' to compliers (heterogeneity), mirroring the implied comparisons in
#' Figure 1 of the manuscript:
#'
#' | Constructor | `theta_A0` \eqn{(\theta_A^{(0)})} | `theta_N1` \eqn{(\theta_N^{(1)})} |
#' |---|---|---|
#' | `restrict_level()` | \eqn{\theta_A^{(0)}} | \eqn{\theta_N^{(1)}} |
#' | `restrict_effect()` | \eqn{\mu_A^{(1)} - \theta_A^{(0)}} | \eqn{\theta_N^{(1)} - \mu_N^{(0)}} |
#' | `restrict_ratio()` | \eqn{\mu_A^{(1)} / \theta_A^{(0)}} | \eqn{\theta_N^{(1)} / \mu_N^{(0)}} |
#' | `restrict_heterogeneity()` | \eqn{\theta_A^{(0)} - \mu_C^{(0)}} | \eqn{\theta_N^{(1)} - \mu_C^{(1)}} |
#'
#' Effects and ratios are always oriented as treated versus untreated, so
#' that e.g. `restrict_effect(lower = 0)` states that the treatment effect is
#' non-negative for the type in question. Ratio restrictions require a
#' non-negative outcome (e.g. a risk ratio for a binary outcome).
#'
#' `lower` and `upper` may be numbers, one-sided formulas, or functions. A
#' formula such as `~ late` or a function such as `function(q) q$late` is
#' evaluated on the identified quantities, and is re-evaluated in every
#' bootstrap replicate so that its uncertainty is propagated. The available
#' quantities are `pi_C`, `pi_A`, `pi_N`, `mu_C0`, `mu_C1`, `mu_A1`, `mu_N0`,
#' and `late` (\eqn{\mu_C^{(1)} - \mu_C^{(0)}}).
#'
#' Plain numbers can also be passed directly to `theta_A0` / `theta_N1`:
#' a single number fixes the parameter, and a length-2 vector is equivalent
#' to `restrict_level(lower, upper)`.
#'
#' @param lower,upper Lower and upper limit on the chosen scale. Numbers,
#'   one-sided formulas, or functions of the identified quantities.
#'   Defaults to an unrestricted limit.
#'
#' @return An object of class `ivate_restriction`.
#'
#' @examples
#' # Treatment increases the outcome among never-takers at most 1.2-fold
#' restrict_ratio(upper = 1.2)
#'
#' # Treatment effect among always-takers is non-negative
#' restrict_effect(lower = 0)
#'
#' # Effect among always-takers at least as large as among compliers
#' restrict_effect(lower = ~ late)
#'
#' @name restrictions
NULL

new_restriction <- function(scale, lower, upper) {
  check_limit <- function(x, arg) {
    ok <- (is.numeric(x) && length(x) == 1L && !is.na(x)) ||
      inherits(x, "formula") || is.function(x)
    if (!ok) {
      stop(sprintf("`%s` must be a single number, a one-sided formula, or a function.", arg),
           call. = FALSE)
    }
  }
  check_limit(lower, "lower")
  check_limit(upper, "upper")
  if (is.numeric(lower) && is.numeric(upper) && lower > upper) {
    stop("`lower` must not exceed `upper`.", call. = FALSE)
  }
  if (scale == "ratio" && is.numeric(lower) && lower < 0) {
    stop("Ratio restrictions require `lower` >= 0.", call. = FALSE)
  }
  structure(list(scale = scale, lower = lower, upper = upper), class = "ivate_restriction")
}

#' @rdname restrictions
#' @export
restrict_level <- function(lower = -Inf, upper = Inf) new_restriction("level", lower, upper)

#' @rdname restrictions
#' @export
restrict_effect <- function(lower = -Inf, upper = Inf) new_restriction("effect", lower, upper)

#' @rdname restrictions
#' @export
restrict_ratio <- function(lower = 0, upper = Inf) new_restriction("ratio", lower, upper)

#' @rdname restrictions
#' @export
restrict_heterogeneity <- function(lower = -Inf, upper = Inf) {
  new_restriction("heterogeneity", lower, upper)
}

#' Evaluate a restriction limit on the identified quantities
#' @keywords internal
#' @noRd
eval_limit <- function(limit, q) {
  value <- if (is.numeric(limit)) {
    limit
  } else if (inherits(limit, "formula")) {
    eval(limit[[length(limit)]], q_as_list(q), environment(limit))
  } else {
    limit(q_as_list(q))
  }
  if (!is.numeric(value) || length(value) != 1L) {
    stop("Restriction limits must evaluate to a single number.", call. = FALSE)
  }
  value
}

#' Translate one restriction into an interval for the sensitivity parameter
#' @keywords internal
#' @noRd
restriction_interval <- function(r, parameter, q) {
  lo <- eval_limit(r$lower, q)
  up <- eval_limit(r$upper, q)
  is_N <- parameter == "theta_N1"

  switch(r$scale,
    level = c(lo, up),
    effect = if (is_N) q[["mu_N0"]] + c(lo, up) else q[["mu_A1"]] - c(up, lo),
    heterogeneity = if (is_N) q[["mu_C1"]] + c(lo, up) else q[["mu_C0"]] + c(lo, up),
    ratio = {
      ref <- if (is_N) q[["mu_N0"]] else q[["mu_A1"]]
      if (!(ref > 0)) {
        stop("Ratio restrictions require a positive reference mean (", if (is_N) "mu_N0" else "mu_A1",
             ").", call. = FALSE)
      }
      if (lo < 0) stop("Ratio restrictions require `lower` >= 0.", call. = FALSE)
      if (is_N) {
        c(ref * lo, ref * up)
      } else {
        c(if (is.infinite(up)) 0 else ref / up, if (lo == 0) Inf else ref / lo)
      }
    }
  )
}

#' Normalise the user specification of a sensitivity parameter restriction
#'
#' Accepts `NULL`, a number, a length-2 numeric vector, an
#' `ivate_restriction`, or a list of numbers/restrictions.
#' @keywords internal
#' @noRd
as_restriction_list <- function(spec, arg) {
  if (is.null(spec)) return(list())
  if (inherits(spec, "ivate_restriction")) return(list(spec))
  if (is.numeric(spec)) {
    if (anyNA(spec) || !(length(spec) %in% 1:2)) {
      stop(sprintf("Numeric `%s` must be a single value or a length-2 interval.", arg),
           call. = FALSE)
    }
    rng <- range(spec)
    return(list(restrict_level(rng[1], rng[2])))
  }
  if (is.list(spec)) return(unlist(lapply(spec, as_restriction_list, arg = arg), recursive = FALSE))
  stop(sprintf("`%s` must be NULL, numeric, a restriction (see ?restrictions), or a list of these.",
               arg), call. = FALSE)
}

#' Feasible interval for a sensitivity parameter
#'
#' Intersects the outcome range with the intervals implied by all
#' restrictions. Errors if the intersection is empty.
#' @keywords internal
#' @noRd
feasible_interval <- function(restrictions, parameter, q, y_range) {
  interval <- y_range
  for (r in restrictions) {
    ri <- restriction_interval(r, parameter, q)
    interval <- c(max(interval[1], ri[1]), min(interval[2], ri[2]))
  }
  if (interval[1] > interval[2]) {
    stop(sprintf("The restrictions on %s are incompatible (empty feasible set).", parameter),
         call. = FALSE)
  }
  interval
}

#' Describe a restriction in words, e.g. "theta_N1 / mu_N0 <= 1.2"
#' @keywords internal
#' @noRd
describe_restriction <- function(r, parameter) {
  is_N <- parameter == "theta_N1"
  expr <- switch(r$scale,
    level         = parameter,
    effect        = if (is_N) "theta_N1 - mu_N0" else "mu_A1 - theta_A0",
    ratio         = if (is_N) "theta_N1 / mu_N0" else "mu_A1 / theta_A0",
    heterogeneity = if (is_N) "theta_N1 - mu_C1" else "theta_A0 - mu_C0"
  )
  fmt_limit <- function(x) {
    if (is.numeric(x)) fmt_num(x)
    else if (inherits(x, "formula")) paste(deparse(x[[length(x)]]), collapse = "")
    else "f(q)"
  }
  unbounded_lo <- is.numeric(r$lower) && (is.infinite(r$lower) || (r$scale == "ratio" && r$lower == 0))
  unbounded_up <- is.numeric(r$upper) && is.infinite(r$upper)
  if (!unbounded_lo && !unbounded_up && identical(r$lower, r$upper)) {
    return(paste(expr, "=", fmt_limit(r$lower)))
  }
  parts <- c(
    if (!unbounded_lo) paste(fmt_limit(r$lower), "<="),
    expr,
    if (!unbounded_up) paste("<=", fmt_limit(r$upper))
  )
  paste(parts, collapse = " ")
}

#' @export
print.ivate_restriction <- function(x, ...) {
  cat("<ivate restriction on the", x$scale, "scale>\n")
  cat("  N:", describe_restriction(x, "theta_N1"), "\n")
  cat("  A:", describe_restriction(x, "theta_A0"), "\n")
  invisible(x)
}
