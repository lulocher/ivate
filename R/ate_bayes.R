#' Bayesian inference for the ATE with an instrumental variable
#'
#' Bayesian analysis of the average treatment effect for a binary outcome,
#' based on the transparent parametrization of Richardson, Evans, and Robins
#' (2011). The identified part of the model -- the observed distributions
#' \eqn{p(y, a \mid z)} -- receives a Dirichlet prior restricted to the
#' distributions compatible with the IV assumptions (their inequalities
#' (11)-(13), equivalent to \eqn{\pi_C \ge 0} and complier means in
#' \eqn{[0, 1]}), and its posterior is sampled by conjugacy and rejection.
#' The sensitivity parameters \eqn{\theta_A^{(0)}} and \eqn{\theta_N^{(1)}}
#' do not enter the likelihood, and are handled in one of two ways:
#'
#' * **Sensitivity analysis** (numeric `theta_A0` / `theta_N1`): the posterior
#'   of \eqn{\Psi_{ATE}(\theta_A^{(0)}, \theta_N^{(1)})} is computed at each
#'   supplied value (or on the grid of all combinations), together with the
#'   posterior probability that the decision criterion holds. This is
#'   the Bayesian analogue of [ate_sensitivity()] and of Figure 5 in
#'   Richardson, Evans, and Robins (2011).
#' * **Full Bayesian analysis** (priors, see [priors]): each posterior draw of
#'   the identified quantities is paired with a draw of the sensitivity
#'   parameters from their prior given those quantities (eq. 27 of
#'   Richardson, Evans, and Robins, 2011), yielding the posterior of the ATE.
#'
#' In a full Bayesian analysis the posterior of the ATE depends on the prior
#' for the sensitivity parameters even in infinitely large samples, since the
#' data carry no information about them. There is therefore no default
#' prior: it must be chosen deliberately. For transparency, the posterior of
#' the identified bounds on the ATE and (by default) the prior distribution
#' of the ATE are reported alongside.
#'
#' Only binary outcomes are supported. Estimates are posterior medians and
#' intervals are equal-tailed credible intervals; in a sensitivity analysis
#' over several values, simultaneous credible bands are also returned. Call
#' [set.seed()] beforehand for reproducible results.
#'
#' @inheritParams ate_sensitivity
#' @param theta_A0,theta_N1 For each sensitivity parameter that is relevant
#'   in the data (i.e. always-takers / never-takers are present), either
#'   numeric value(s) in \eqn{[0, 1]} for a sensitivity analysis, or a prior
#'   created by [prior_beta()], [prior_uniform()], or [prior_normal()] for a
#'   full Bayesian analysis. If one parameter is given a prior, a single
#'   number for the other fixes it (point-mass prior).
#' @param level Credible level. Defaults to `0.95`.
#' @param prior_draws Logical; in a full Bayesian analysis, also sample the
#'   prior distribution of the ATE (for comparison with the posterior)?
#'   Defaults to `TRUE`.
#'
#' @return An object of class `ivate_bayes` with elements
#'   \describe{
#'     \item{`mode`}{`"sensitivity"` or `"full"`.}
#'     \item{`identified`}{Posterior medians and credible intervals of the
#'       identified quantities.}
#'     \item{`bounds`}{Posterior of the (unrestricted) bounds on the ATE.}
#'     \item{`grid`}{Sensitivity analysis: data frame with one row per value
#'       of the sensitivity parameters, giving the posterior median, pointwise
#'       (`ci_*`) and simultaneous (`band_*`) credible limits of the ATE, and
#'       the posterior probability `prob` that the decision criterion holds.}
#'     \item{`posterior`, `prior_summary`}{Full Bayesian analysis: summaries
#'       of the posterior and prior distribution of the ATE.}
#'     \item{`priors`}{Full Bayesian analysis: the priors on the sensitivity
#'       parameters.}
#'     \item{`prior_ate`}{Full Bayesian analysis: draws of the ATE from its
#'       prior (if `prior_draws = TRUE`).}
#'     \item{`draws`}{Full Bayesian analysis (or a single value in a
#'       sensitivity analysis): draws of the sensitivity parameters and the
#'       ATE.}
#'     \item{`bayes`}{Prior settings and posterior probabilities that the IV
#'       inequalities are violated (a model check).}
#'   }
#'
#' @references
#' Richardson, T. S., Evans, R. J., and Robins, J. M. (2011). Transparent
#' parametrizations of models for potential outcomes. In *Bayesian Statistics
#' 9*, 569--610. Oxford University Press.
#'
#' Locher, L., Richardson, T., and Stensrud, M. (2026). Informative inference
#' on the average treatment effect with instrumental variables.
#'
#' @seealso [priors], [ate_sensitivity()] and [ate_bounds()] with
#'   `inference = "bayes"`.
#'
#' @examples
#' # Lipid data (one-sided non-compliance), Richardson, Evans, and Robins (2011)
#' lipid <- expand_counts(data.frame(
#'   z = c(0, 0, 0, 0, 1, 1, 1, 1),
#'   a = c(0, 0, 1, 1, 0, 0, 1, 1),
#'   y = c(0, 1, 0, 1, 0, 1, 0, 1),
#'   n = c(158, 14, 0, 0, 52, 12, 23, 78)
#' ))
#' set.seed(1)
#'
#' # Sensitivity analysis: posterior of the ATE as a function of theta_N1
#' sens <- ate_bayes(lipid, y ~ a | z, theta_N1 = seq(0, 1, by = 0.05), n_draws = 1000)
#' sens
#' plot(sens)
#'
#' # Full Bayesian analysis with a prior centred at the complier effect
#' full <- ate_bayes(lipid, y ~ a | z,
#'                   theta_N1 = prior_normal(mean = ~ late, sd = 0.1, on = "effect"),
#'                   n_draws = 1000)
#' full
#' plot(full)
#'
#' @export
ate_bayes <- function(data, formula,
                      theta_A0    = NULL,
                      theta_N1    = NULL,
                      cutoff      = 0,
                      direction   = c("greater", "less"),
                      n_draws     = 4000,
                      prior_obs   = 1,
                      level       = 0.95,
                      prior_draws = TRUE) {
  direction <- match.arg(direction)
  if (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff)) {
    stop("`cutoff` must be a single number.", call. = FALSE)
  }
  check_bayes_args(n_draws, prior_obs)
  check_inference_args(TRUE, 2, level)

  d <- prepare_iv_data(data, formula)
  check_binary_outcome(d$y)
  q <- identify_quantities(d$y, d$a, d$z)
  check_identified(q, c(0, 1))
  setting <- iv_setting(q)
  if (setting == "full compliance") {
    stop(sprintf("No non-compliance in the data: the ATE is point identified (ATE = %s).",
                 fmt_num(q[["mu_C1"]] - q[["mu_C0"]])), call. = FALSE)
  }

  specs <- bayes_specs(setting, theta_A0, theta_N1)
  post  <- posterior_q(d$y, d$a, d$z, setting, n_draws, prior_obs)
  qd    <- post$q

  bounds <- bounds_rows(evaluate_stat(bounds_stat_fun(list(), list(), c(0, 1)), q, qd, level, "median"))

  out <- list(
    call       = match.call(),
    vars       = d$vars,
    n          = length(d$y),
    setting    = setting,
    mode       = specs$mode,
    cutoff     = cutoff,
    direction  = direction,
    level      = level,
    q          = q,
    identified = identified_table(q, qd, level, "median"),
    bounds     = bounds,
    inference  = "bayes",
    bayes      = list(prior_obs = prior_obs, n_draws = n_draws,
                      violation = post$violation, n_proposed = post$n_proposed),
    q_draws    = qd
  )

  if (specs$mode == "sensitivity") {
    grid <- expand.grid(theta_A0 = specs$values$theta_A0 %||% NA_real_,
                        theta_N1 = specs$values$theta_N1 %||% NA_real_)
    summ <- ate_grid_summary(qd, grid, cutoff, direction, level)
    out$grid  <- cbind(grid, summ$table)
    out$draws <- summ$draws
  } else {
    out$priors <- specs$priors
    post_ate   <- ate_prior_draws(qd, specs$priors)
    out$draws  <- post_ate
    out$posterior <- ate_summary(post_ate$ate, cutoff, direction, level)
    if (prior_draws) {
      prior_q   <- sample_dirichlet_q(NULL, setting, n_draws, prior_obs)$q
      prior_ate <- ate_prior_draws(prior_q, specs$priors)
      out$prior_summary <- ate_summary(prior_ate$ate, cutoff, direction, level)
      out$prior_ate     <- prior_ate$ate
    }
  }
  structure(out, class = "ivate_bayes")
}

#' Validate the sensitivity-parameter specifications and determine the mode
#' @keywords internal
#' @noRd
bayes_specs <- function(setting, theta_A0, theta_N1) {
  present <- c(theta_A0 = setting %in% c("two-sided", "always-takers only"),
               theta_N1 = setting %in% c("two-sided", "never-takers only"))
  specs <- list(theta_A0 = theta_A0, theta_N1 = theta_N1)
  for (p in names(specs)) {
    if (!present[[p]] && !is.null(specs[[p]])) {
      warning(sprintf("There are no %s; `%s` is ignored.",
                      if (p == "theta_A0") "always-takers" else "never-takers", p), call. = FALSE)
      specs[p] <- list(NULL)
    }
    if (present[[p]] && is.null(specs[[p]])) {
      stop(sprintf(paste0(
        "Supply `%s`: numeric value(s) in [0, 1] for a sensitivity analysis, or a prior ",
        "(see ?priors) for a full Bayesian analysis."), p), call. = FALSE)
    }
  }
  specs <- specs[present]

  is_prior <- vapply(specs, function(s) inherits(s, c("ivate_prior", "ivate_restriction")), logical(1))
  if (any(is_prior)) {
    priors <- lapply(names(specs), function(p) {
      s <- specs[[p]]
      if (is.numeric(s) && length(s) != 1L) {
        stop(sprintf(paste0("In a full Bayesian analysis a number fixes `%s`; supply a single ",
                            "value or a prior."), p), call. = FALSE)
      }
      as_prior(s, p)
    })
    names(priors) <- names(specs)
    return(list(mode = "full", priors = priors))
  }

  values <- lapply(names(specs), function(p) {
    v <- specs[[p]]
    if (!is.numeric(v) || length(v) == 0 || anyNA(v) || any(v < 0 | v > 1)) {
      stop(sprintf("`%s` must be numeric value(s) in [0, 1] or a prior (see ?priors).", p),
           call. = FALSE)
    }
    sort(unique(v))
  })
  names(values) <- names(specs)
  list(mode = "sensitivity", values = values)
}

#' Linear representation of the ATE on each draw
#'
#' ATE = intercept + slope_A * theta_A0 + slope_N * theta_N1, with absent
#' types contributing zero.
#' @keywords internal
#' @noRd
ate_linear <- function(qd) {
  pA <- qd[, "pi_A"]; pN <- qd[, "pi_N"]
  list(
    intercept = qd[, "pi_C"] * (qd[, "mu_C1"] - qd[, "mu_C0"]) +
      ifelse(pA > 0, pA * qd[, "mu_A1"], 0) - ifelse(pN > 0, pN * qd[, "mu_N0"], 0),
    slope_A = -pA,
    slope_N = pN
  )
}

#' Posterior summaries of the ATE on a grid of sensitivity parameters
#'
#' Processes the grid in chunks to bound memory. Simultaneous bands use the
#' maximum standardised deviation from the posterior median across the grid.
#' @keywords internal
#' @noRd
ate_grid_summary <- function(qd, grid, cutoff, direction, level, chunk = 500L) {
  lin <- ate_linear(qd)
  tA <- ifelse(is.na(grid$theta_A0), 0, grid$theta_A0)
  tN <- ifelse(is.na(grid$theta_N1), 0, grid$theta_N1)
  G <- nrow(grid)
  alpha <- (1 - level) / 2
  ate_chunk <- function(idx) {
    lin$intercept + outer(lin$slope_A, tA[idx]) + outer(lin$slope_N, tN[idx])
  }
  holds <- function(A) if (direction == "greater") A > cutoff else A < cutoff

  chunks <- split(seq_len(G), ceiling(seq_len(G) / chunk))
  tab <- do.call(rbind, lapply(chunks, function(idx) {
    A <- ate_chunk(idx)
    qs <- apply(A, 2, stats::quantile, probs = c(0.5, alpha, 1 - alpha), names = FALSE)
    data.frame(median = qs[1, ], ci_lower = qs[2, ], ci_upper = qs[3, ],
               sd = apply(A, 2, stats::sd), prob = colMeans(holds(A)))
  }))

  if (G > 1) {
    max_dev <- rep(0, nrow(qd))
    for (idx in chunks) {
      A <- ate_chunk(idx)
      dev <- abs(sweep(A, 2, tab$median[idx])) / rep(pmax(tab$sd[idx], 1e-12), each = nrow(A))
      max_dev <- pmax(max_dev, apply(dev, 1, max))
    }
    crit <- stats::quantile(max_dev, level, names = FALSE)
    tab$band_lower <- tab$median - crit * tab$sd
    tab$band_upper <- tab$median + crit * tab$sd
  } else {
    tab$band_lower <- tab$ci_lower
    tab$band_upper <- tab$ci_upper
  }
  tab$sd <- NULL
  tab <- tab[, c("median", "ci_lower", "ci_upper", "band_lower", "band_upper", "prob")]
  rownames(tab) <- NULL

  draws <- if (G == 1) {
    data.frame(theta_A0 = grid$theta_A0, theta_N1 = grid$theta_N1, ate = as.vector(ate_chunk(1L)))
  } else NULL
  list(table = tab, draws = draws)
}

#' Draw the sensitivity parameters from their priors and compute the ATE
#' @keywords internal
#' @noRd
ate_prior_draws <- function(qd, priors) {
  n  <- nrow(qd)
  tA <- if (!is.null(priors$theta_A0)) sample_prior(priors$theta_A0, "theta_A0", qd) else rep(NA_real_, n)
  tN <- if (!is.null(priors$theta_N1)) sample_prior(priors$theta_N1, "theta_N1", qd) else rep(NA_real_, n)
  lin <- ate_linear(qd)
  ate <- lin$intercept + lin$slope_A * ifelse(is.na(tA), 0, tA) + lin$slope_N * ifelse(is.na(tN), 0, tN)

  failed <- (!is.null(priors$theta_A0) & is.na(tA)) | (!is.null(priors$theta_N1) & is.na(tN))
  if (all(failed)) {
    stop("The prior(s) could not be evaluated on any draw (e.g. empty feasible set).", call. = FALSE)
  }
  if (any(failed)) {
    warning(sprintf("The prior(s) could not be evaluated on %d of %d draws, which were dropped.",
                    sum(failed), n), call. = FALSE)
  }
  data.frame(theta_A0 = tA, theta_N1 = tN, ate = ate)[!failed, , drop = FALSE]
}

#' Summary of draws of the ATE
#' @keywords internal
#' @noRd
ate_summary <- function(ate, cutoff, direction, level) {
  alpha <- (1 - level) / 2
  qs <- stats::quantile(ate, c(0.5, alpha, 1 - alpha), names = FALSE)
  data.frame(median = qs[1], mean = mean(ate), ci_lower = qs[2], ci_upper = qs[3],
             prob = mean(if (direction == "greater") ate > cutoff else ate < cutoff))
}

#' @export
print.ivate_bayes <- function(x, ...) {
  print_header(x, "Bayesian analysis of the ATE with an instrumental variable")
  rel <- if (x$direction == "greater") ">" else "<"
  cat(sprintf("  Criterion: ATE %s %s\n", rel, fmt_num(x$cutoff)))
  print_bayes_info(x)

  lo <- x$bounds[x$bounds$bound == "lower", ]
  up <- x$bounds[x$bounds$bound == "upper", ]
  cat(sprintf("\nBounds on the ATE (identified): [%s, %s]\n  %s lower %s, upper %s\n",
              fmt_num(lo$estimate), fmt_num(up$estimate), ci_header(x$level, "bayes"),
              fmt_ci(lo$ci_lower, lo$ci_upper), fmt_ci(up$ci_lower, up$ci_upper)))

  prob_label <- sprintf("P(ATE %s %s)", rel, fmt_num(x$cutoff))
  if (x$mode == "sensitivity") {
    cat("\nPosterior of the ATE given the sensitivity parameters:\n")
    g <- x$grid
    show <- data.frame(check.names = FALSE,
      g[, c("theta_A0", "theta_N1")[c(any(!is.na(g$theta_A0)), any(!is.na(g$theta_N1)))], drop = FALSE],
      median = fmt_num(g$median),
      interval = fmt_ci(g$ci_lower, g$ci_upper),
      prob = fmt_num(g$prob)
    )
    names(show)[names(show) == "interval"] <- ci_header(x$level, "bayes")
    names(show)[names(show) == "prob"] <- prob_label
    n_show <- min(nrow(show), 15L)
    print(show[seq_len(n_show), , drop = FALSE], row.names = FALSE, right = FALSE)
    if (nrow(show) > n_show) cat(sprintf("  ... %d more rows in `$grid`\n", nrow(show) - n_show))
  } else {
    cat("\nPriors on the sensitivity parameters:\n")
    for (p in names(x$priors)) cat("  ", describe_prior(x$priors[[p]], p), "\n", sep = "")
    show_row <- function(s, label) {
      cat(sprintf("  %-10s median %s, %s %s, %s = %s\n", label, fmt_num(s$median),
                  ci_header(x$level, "bayes"), fmt_ci(s$ci_lower, s$ci_upper), prob_label,
                  fmt_num(s$prob)))
    }
    cat("\nATE:\n")
    show_row(x$posterior, "Posterior:")
    if (!is.null(x$prior_summary)) show_row(x$prior_summary, "Prior:")
    cat("\nNote: the data carry no information about the sensitivity parameters, so this\n",
        "posterior depends on their prior even in large samples.\n", sep = "")
  }
  invisible(x)
}
