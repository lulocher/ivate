#' Sensitivity analysis for the ATE under IV assumptions
#'
#' Under the standard IV assumptions (consistency, relevance, independence,
#' exclusion restriction, and monotonicity), point identification of the
#' average treatment effect (ATE) is precluded by two unidentified quantities:
#' \eqn{\theta_A^{(0)} = E[Y(0) \mid A]}, the mean potential outcome under
#' control among always-takers, and \eqn{\theta_N^{(1)} = E[Y(1) \mid N]}, the
#' mean potential outcome under treatment among never-takers. This function
#' derives the values of these sensitivity parameters at which the decision
#' criterion \eqn{\Psi_{ATE} > x} (or \eqn{\Psi_{ATE} < x}) switches, and
#' expresses them as implied comparisons with identified quantities (Section
#' 2.4 of Locher, Richardson, and Stensrud, 2026).
#'
#' The ATE decomposes as (equation (7))
#' \deqn{\Psi_{ATE} = \pi_C(\mu_C^{(1)} - \mu_C^{(0)})
#'   + \pi_A(\mu_A^{(1)} - \theta_A^{(0)})
#'   + \pi_N(\theta_N^{(1)} - \mu_N^{(0)}),}
#' where all quantities except \eqn{\theta_A^{(0)}} and \eqn{\theta_N^{(1)}}
#' are identified.
#'
#' **One-sided non-compliance.** If only never-takers (or only always-takers)
#' are present, a single scalar threshold \eqn{\tau} fully characterises the
#' decision problem (equation (10)).
#'
#' **Two-sided non-compliance.** The criterion defines a half-plane in
#' \eqn{(\theta_A^{(0)}, \theta_N^{(1)})}; its boundary \eqn{\mathcal{B}} is
#' returned both on the level scale and on the scale of the implied
#' type-specific treatment effects (equations (8) and (9)). Following the
#' two-step strategy of Section 3.2, one of the two parameters can be
#' restricted via `theta_A0` or `theta_N1`: the threshold for the other
#' parameter is then evaluated at the least favourable value within the
#' restriction. For instance, `theta_A0 = 1` fixes
#' \eqn{\theta_A^{(0)} = 1} (the worst case for \eqn{\Psi_{ATE} > x} with a
#' binary outcome), and `theta_A0 = restrict_effect(lower = 0)` assumes a
#' non-negative effect among always-takers, whose least favourable value is
#' \eqn{\theta_A^{(0)} = \mu_A^{(1)}}.
#'
#' **Implied comparisons.** At the threshold, the implied treatment effect
#' and ratio within the unidentified type and the implied outcome
#' heterogeneity relative to compliers are reported, each next to an
#' identified benchmark (Figure 1 of the manuscript). Ratios are only
#' reported for non-negative outcomes.
#'
#' **Inference.** With `inference = "bootstrap"` (default), estimates are
#' plug-in estimates and intervals are percentile intervals from the
#' nonparametric bootstrap (Section 2.6). With `inference = "bayes"` (binary
#' outcomes only), the identified part of the model is given a Dirichlet prior
#' on \eqn{p(y, a \mid z)} restricted to the distributions compatible with the
#' IV assumptions (Richardson, Evans, and Robins, 2011); estimates are
#' posterior medians and intervals are equal-tailed credible intervals. The
#' posterior probability that the IV inequalities are violated is reported as
#' a model check. Call [set.seed()] beforehand for reproducible intervals.
#'
#' **Relevance.** Thresholds, boundaries, bounds, and the ATE depend on the
#' compliers only through \eqn{\pi_C \mu_C^{(1)}} and \eqn{\pi_C \mu_C^{(0)}},
#' which remain identified without compliers. The procedures therefore also
#' apply when \eqn{\pi_C = 0}; the complier means, the LATE, and comparisons
#' involving them are then undefined (`NA`), and a warning is issued.
#'
#' **Model check.** In finite samples, the plug-in estimates can violate the
#' testable restrictions of the IV model (\eqn{\pi_C \ge 0} and complier means
#' within the outcome range; for binary outcomes these are the inequalities
#' (11)-(13) of Richardson, Evans, and Robins, 2011), yielding values
#' incompatible with the model. This is checked in the original sample (with
#' a warning) and in every bootstrap sample; the share of violating bootstrap
#' samples is reported and a warning is issued if it exceeds 5 percent, in which
#' case the bootstrap intervals may be unreliable. A violation in the sample
#' may indicate that the IV assumptions fail, but can also arise from sampling
#' variability, e.g. with weak instruments, rare outcomes, or small samples.
#'
#' @param data A data frame.
#' @param formula A formula `outcome ~ treatment | instrument` (or a string).
#'   Treatment and instrument must be binary (0/1), with the instrument
#'   coded such that it increases treatment uptake.
#' @param cutoff The policy-relevant threshold \eqn{x} in the decision
#'   criterion. Defaults to `0`.
#' @param direction `"greater"` for the criterion \eqn{\Psi_{ATE} > x}
#'   (default) or `"less"` for \eqn{\Psi_{ATE} < x}. This affects the
#'   reported inequalities and which value of a restricted parameter is least
#'   favourable, not the location of the threshold.
#' @param theta_A0,theta_N1 Optional restriction on one sensitivity parameter
#'   (two-sided non-compliance only): a single number (fixes the parameter), a
#'   length-2 interval, a restriction created by [restrict_level()],
#'   [restrict_effect()], [restrict_ratio()], or [restrict_heterogeneity()],
#'   or a list of these (intersected).
#' @param y_range Numeric vector `c(y_min, y_max)` giving the range of the
#'   outcome. Defaults to `c(0, 1)` for binary outcomes and to the observed
#'   range otherwise.
#' @param ci Logical; compute bootstrap confidence intervals? Defaults to
#'   `TRUE`. Ignored for Bayesian inference, which always samples the
#'   posterior.
#' @param n_boot Number of bootstrap replicates. Defaults to `1000`.
#' @param level Confidence level. Defaults to `0.95`.
#' @param boot_strata `"none"` (default) resamples observations iid;
#'   `"instrument"` resamples within instrument arms, keeping the arm sizes
#'   fixed as in a randomised trial.
#' @param inference `"bootstrap"` (default) or `"bayes"` (binary outcomes
#'   only); see Details.
#' @param n_draws Number of posterior draws for Bayesian inference. Defaults
#'   to `4000`.
#' @param prior_obs Concentration parameter of the Dirichlet prior on each
#'   cell of \eqn{p(y, a \mid z)} for Bayesian inference. Defaults to `1`
#'   (uniform).
#'
#' @return An object of class `ivate_sensitivity` with elements
#'   \describe{
#'     \item{`setting`}{`"two-sided"`, `"never-takers only"`, or
#'       `"always-takers only"`.}
#'     \item{`identified`}{Data frame of identified quantities (type shares,
#'       potential outcome means, LATE) with confidence intervals.}
#'     \item{`threshold`}{Data frame with the threshold \eqn{\tau} on the
#'       remaining sensitivity parameter and its confidence interval, or
#'       `NULL` in the unrestricted two-sided case.}
#'     \item{`fixed`}{Data frame describing the restricted parameter and its
#'       least favourable value, or `NULL`.}
#'     \item{`implied`}{Data frame of implied comparisons at the threshold
#'       with identified benchmarks, or `NULL`.}
#'     \item{`boundary`}{Data frame with intercepts and slopes of the boundary
#'       \eqn{\mathcal{B}} (two-sided case), or `NULL`.}
#'     \item{`draws`}{Matrix of bootstrap replicates or posterior draws of
#'       the identified quantities, or `NULL`.}
#'     \item{`bayes`}{For Bayesian inference, the prior settings and the
#'       posterior probabilities that the IV inequalities (11)-(13) of
#'       Richardson, Evans, and Robins (2011) are violated.}
#'     \item{`model_check`}{Violations of the IV model restrictions in the
#'       original sample (`sample`, logical) and the share of bootstrap
#'       samples violating each restriction (`bootstrap`); see Details.}
#'   }
#'   together with the inputs. Use [plot()][plot.ivate_sensitivity()] to
#'   visualise the result.
#'
#' @references
#' Locher, L., Richardson, T., and Stensrud, M. (2026). Informative inference
#' on the average treatment effect with instrumental variables.
#'
#' Angrist, J. D., Imbens, G. W., and Rubin, D. B. (1996). Identification of
#' causal effects using instrumental variables. *Journal of the American
#' Statistical Association*, 91(434), 444--455.
#'
#' Richardson, T. S., Evans, R. J., and Robins, J. M. (2011). Transparent
#' parametrizations of models for potential outcomes. In *Bayesian Statistics
#' 9*, 569--610. Oxford University Press.
#'
#' @seealso [ate_bounds()] for bounds under restrictions, [restrictions],
#'   [ate_bayes()] for the posterior of the ATE.
#'
#' @examples
#' set.seed(1)
#' n <- 2000
#' z <- rbinom(n, 1, 0.5)
#'
#' # One-sided non-compliance (never-takers only)
#' a <- ifelse(z == 1, rbinom(n, 1, 0.6), 0)
#' y <- rbinom(n, 1, 0.3 + 0.1 * a)
#' res <- ate_sensitivity(data.frame(y, a, z), y ~ a | z, n_boot = 200)
#' res
#'
#' # Two-sided non-compliance: boundary, then fix theta_A0 at its worst case
#' a2 <- ifelse(z == 1, rbinom(n, 1, 0.7), rbinom(n, 1, 0.2))
#' y2 <- rbinom(n, 1, 0.3 + 0.2 * a2)
#' dat2 <- data.frame(y = y2, a = a2, z = z)
#' ate_sensitivity(dat2, y ~ a | z, ci = FALSE)
#' ate_sensitivity(dat2, y ~ a | z, theta_A0 = 1, ci = FALSE)
#' ate_sensitivity(dat2, y ~ a | z, theta_A0 = restrict_effect(lower = 0), ci = FALSE)
#'
#' # Bayesian inference (binary outcomes)
#' ate_sensitivity(dat2, y ~ a | z, theta_A0 = 1, inference = "bayes", n_draws = 1000)
#'
#' @export
ate_sensitivity <- function(data, formula,
                            cutoff      = 0,
                            direction   = c("greater", "less"),
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
  direction   <- match.arg(direction)
  boot_strata <- match.arg(boot_strata)
  inference   <- match.arg(inference)
  if (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff)) {
    stop("`cutoff` must be a single number.", call. = FALSE)
  }
  check_inference_args(ci, n_boot, level)
  if (inference == "bayes") check_bayes_args(n_draws, prior_obs)

  d       <- prepare_iv_data(data, formula)
  if (inference == "bayes") check_binary_outcome(d$y, y_range)
  y_range <- resolve_y_range(d$y, y_range)
  q       <- identify_quantities(d$y, d$a, d$z)
  check_identified(q, y_range)
  setting <- iv_setting(q)
  if (setting == "full compliance") {
    stop(sprintf("No non-compliance in the data: the ATE is point identified (ATE = %s).",
                 fmt_num(q[["mu_C1"]] - q[["mu_C0"]])), call. = FALSE)
  }

  restr_A <- as_restriction_list(theta_A0, "theta_A0")
  restr_N <- as_restriction_list(theta_N1, "theta_N1")
  roles   <- sensitivity_roles(setting, restr_A, restr_N)
  restr   <- if (identical(roles$fixed, "theta_A0")) restr_A else restr_N

  unc    <- uncertainty_draws(d, setting, ci, inference, n_boot, boot_strata, n_draws, prior_obs)
  draws  <- unc$draws
  mcheck <- model_check(q, if (inference == "bootstrap") draws else NULL, y_range, is_binary(d$y))
  center <- if (inference == "bayes") "median" else "plugin"

  identified <- identified_table(q, draws, level, center)

  boundary <- NULL
  if (setting == "two-sided") {
    boundary <- boundary_table(
      evaluate_stat(function(q) boundary_coefficients(q, cutoff), q, draws, level, center),
      direction
    )
  }

  threshold <- fixed <- implied <- NULL
  stats_1d <- NULL
  if (!is.null(roles$target)) {
    stat_fun <- threshold_stat_fun(roles$target, roles$fixed, restr, cutoff, direction, y_range)
    stats_1d <- evaluate_stat(stat_fun, q, draws, level, center)
    tab      <- split_stat_table(stats_1d)

    threshold <- data.frame(
      parameter = roles$target,
      relation  = criterion_relation(roles$target, direction),
      tab["threshold", -1, drop = FALSE],
      row.names = NULL
    )
    if (!is.null(roles$fixed)) {
      fixed <- data.frame(
        parameter   = roles$fixed,
        restriction = paste(vapply(restr, describe_restriction, character(1),
                                   parameter = roles$fixed), collapse = "; "),
        tab["fixed_value", -1, drop = FALSE],
        row.names = NULL
      )
      names(fixed)[names(fixed) == "estimate"] <- "value"
    }
    implied <- implied_table(tab, roles$target, nonnegative = y_range[1] >= 0)
  }

  if (!is.null(draws)) warn_dropped(stats_1d, n_draws = nrow(draws))

  structure(list(
    call         = match.call(),
    vars         = d$vars,
    n            = length(d$y),
    setting      = setting,
    cutoff       = cutoff,
    direction    = direction,
    y_range      = y_range,
    q            = q,
    identified   = identified,
    threshold    = threshold,
    fixed        = fixed,
    implied      = implied,
    boundary     = boundary,
    target       = roles$target,
    restrictions = if (is.null(roles$fixed)) NULL else stats::setNames(list(restr), roles$fixed),
    draws        = draws,
    inference    = inference,
    bayes        = unc$bayes,
    model_check  = mcheck,
    level        = level
  ), class = "ivate_sensitivity")
}

#' Determine which parameter receives a threshold and which is restricted
#' @keywords internal
#' @noRd
sensitivity_roles <- function(setting, restr_A, restr_N) {
  has_A <- length(restr_A) > 0
  has_N <- length(restr_N) > 0
  switch(setting,
    "two-sided" = {
      if (has_A && has_N) {
        stop(paste0("Restrict at most one of `theta_A0` and `theta_N1`: the threshold is ",
                    "derived for the other. Use ate_bounds() to restrict both."), call. = FALSE)
      }
      if (has_A) list(target = "theta_N1", fixed = "theta_A0")
      else if (has_N) list(target = "theta_A0", fixed = "theta_N1")
      else list(target = NULL, fixed = NULL)
    },
    "never-takers only" = {
      if (has_N) {
        stop(paste0("With never-takers only, theta_N1 is the sole sensitivity parameter. ",
                    "Use ate_bounds() to propagate restrictions on it."), call. = FALSE)
      }
      if (has_A) warning("There are no always-takers; `theta_A0` is ignored.", call. = FALSE)
      list(target = "theta_N1", fixed = NULL)
    },
    "always-takers only" = {
      if (has_A) {
        stop(paste0("With always-takers only, theta_A0 is the sole sensitivity parameter. ",
                    "Use ate_bounds() to propagate restrictions on it."), call. = FALSE)
      }
      if (has_N) warning("There are no never-takers; `theta_N1` is ignored.", call. = FALSE)
      list(target = "theta_A0", fixed = NULL)
    }
  )
}

#' Statistic function for a one-dimensional threshold and implied comparisons
#' @keywords internal
#' @noRd
threshold_stat_fun <- function(target, fixed, restr, cutoff, direction, y_range) {
  function(q) {
    other <- NA_real_
    if (!is.null(fixed)) {
      other <- least_favourable(feasible_interval(restr, fixed, q, y_range), fixed, direction)
    }
    tau <- if (target == "theta_N1") threshold_N1(q, cutoff, other) else threshold_A0(q, cutoff, other)
    c(threshold = tau, fixed_value = other, implied_comparisons(q, tau, target))
  }
}

#' Relation under which the decision criterion holds
#'
#' The ATE increases in theta_N1 and in the type-specific effects, and
#' decreases in theta_A0.
#' @keywords internal
#' @noRd
criterion_relation <- function(lhs, direction) {
  increasing <- lhs != "theta_A0"
  if (xor(increasing, direction == "less")) ">" else "<"
}

#' Convert an evaluate_stat() result to a data frame indexed by term
#' @keywords internal
#' @noRd
split_stat_table <- function(tab) {
  rownames(tab) <- tab$term
  tab
}

#' Table of identified quantities
#' @keywords internal
#' @noRd
identified_table <- function(q, draws, level, center = "plugin") {
  tab <- evaluate_stat(function(q) c(q, late = q[["mu_C1"]] - q[["mu_C0"]]), q, draws, level, center)
  tab <- tab[!is.na(tab$estimate) & !tab$term %in% c("pimu_C0", "pimu_C1"), , drop = FALSE]
  data.frame(quantity = tab$term, description = unname(q_descriptions[tab$term]),
             tab[, -1, drop = FALSE], row.names = NULL, stringsAsFactors = FALSE)
}

#' Table of implied comparisons
#' @keywords internal
#' @noRd
implied_table <- function(tab, target, nonnegative) {
  labels <- implied_labels[[target]]
  rows <- lapply(labels$comparison, function(cmp) {
    imp <- tab[paste0(cmp, "_implied"), , drop = FALSE]
    ben <- tab[paste0(cmp, "_benchmark"), , drop = FALSE]
    out <- data.frame(implied_value = imp$estimate, benchmark_value = ben$estimate)
    if ("ci_lower" %in% names(tab)) {
      out$implied_ci_lower   <- imp$ci_lower
      out$implied_ci_upper   <- imp$ci_upper
      out$benchmark_ci_lower <- ben$ci_lower
      out$benchmark_ci_upper <- ben$ci_upper
    }
    out
  })
  out <- cbind(labels, do.call(rbind, rows))
  out$relative <- ifelse(out$comparison == "ratio", NA_real_, out$implied_value / out$benchmark_value)
  if (!nonnegative) out <- out[out$comparison != "ratio", , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Table of boundary coefficients
#' @keywords internal
#' @noRd
boundary_table <- function(tab, direction) {
  specs <- data.frame(
    scale  = c("level", "level", "effect", "effect"),
    lhs    = c("theta_N1", "theta_A0", "theta_N1 - mu_N0", "mu_A1 - theta_A0"),
    rhs    = c("theta_A0", "theta_N1", "mu_A1 - theta_A0", "theta_N1 - mu_N0"),
    prefix = c("level_N1", "level_A0", "effect_N", "effect_A"),
    stringsAsFactors = FALSE
  )
  rownames(tab) <- tab$term
  rows <- lapply(seq_len(nrow(specs)), function(i) {
    it <- tab[paste0(specs$prefix[i], "_intercept"), , drop = FALSE]
    sl <- tab[paste0(specs$prefix[i], "_slope"), , drop = FALSE]
    out <- data.frame(
      scale     = specs$scale[i],
      lhs       = specs$lhs[i],
      relation  = criterion_relation(if (specs$lhs[i] == "theta_A0") "theta_A0" else "other", direction),
      intercept = it$estimate,
      slope     = sl$estimate,
      rhs       = specs$rhs[i],
      stringsAsFactors = FALSE
    )
    if ("ci_lower" %in% names(tab)) {
      out$intercept_ci_lower <- it$ci_lower
      out$intercept_ci_upper <- it$ci_upper
      out$slope_ci_lower     <- sl$ci_lower
      out$slope_ci_upper     <- sl$ci_upper
    }
    out
  })
  do.call(rbind, rows)
}
