#' Plot a sensitivity analysis
#'
#' Visualises the result of [ate_sensitivity()] in the style of the
#' threshold and boundary plots of the examples in Locher, Richardson, and
#' Stensrud (2026).
#'
#' * **Threshold plot** (one-sided non-compliance, or two-sided with one
#'   parameter restricted): the ATE as a function of the remaining
#'   sensitivity parameter, with a pointwise confidence (credible) band, the
#'   cutoff as a horizontal line, and the threshold \eqn{\tau} (dashed) with
#'   its confidence interval (shaded).
#' * **Boundary plot** (two-sided non-compliance): a heatmap of the ATE over
#'   \eqn{(\theta_A^{(0)}, \theta_N^{(1)})} with the boundary \eqn{\mathcal{B}}
#'   and its pointwise confidence (credible) band.
#'
#' For results obtained with `inference = "bayes"`, curves show posterior
#' medians and bands are pointwise credible bands.
#'
#' In both plots, the identified potential outcome means are marked as
#' reference ticks on the secondary axes; nearby labels are merged.
#'
#' Plots are `ggplot` objects and can be modified with the usual `ggplot2`
#' syntax, or combined into multi-panel figures with e.g. the \pkg{patchwork}
#' package.
#'
#' @param x An `ivate_sensitivity` object.
#' @param type `"auto"` (default) draws the threshold plot when a scalar
#'   threshold is available and the boundary plot otherwise; `"threshold"` or
#'   `"boundary"` force a type (the boundary plot requires two-sided
#'   non-compliance). `"comparisons"` draws the roster of estimated, assumed,
#'   and implied quantities (the comparison figure of the manuscript); see
#'   [plot_comparisons()], which also offers intervals.
#' @param range Numeric vector of length 2 giving the plotted range of the
#'   sensitivity parameter(s). Defaults to the outcome range.
#' @param ci Logical; draw confidence (credible) bands if available?
#'   Defaults to `TRUE`.
#' @param n_grid Number of grid points per axis.
#' @param label_tol Reference labels closer than `label_tol` times the width
#'   of `range` are merged.
#' @param ... Passed to [plot_comparisons()] for `type = "comparisons"`
#'   (e.g. `digits`); otherwise unused.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' set.seed(1)
#' n <- 2000
#' z <- rbinom(n, 1, 0.5)
#' a <- ifelse(z == 1, rbinom(n, 1, 0.6), 0)
#' y <- rbinom(n, 1, 0.3 + 0.1 * a)
#' res <- ate_sensitivity(data.frame(y, a, z), y ~ a | z, n_boot = 200)
#' plot(res)
#' plot(res, range = c(0.2, 0.6))
#'
#' a2 <- ifelse(z == 1, rbinom(n, 1, 0.7), rbinom(n, 1, 0.2))
#' y2 <- rbinom(n, 1, 0.3 + 0.2 * a2)
#' res2 <- ate_sensitivity(data.frame(y = y2, a = a2, z = z), y ~ a | z, n_boot = 200)
#' plot(res2)
#' plot(res2, type = "comparisons")
#'
#' @export
plot.ivate_sensitivity <- function(x, type = c("auto", "threshold", "boundary", "comparisons"),
                                   range = NULL, ci = TRUE, n_grid = 101,
                                   label_tol = 0.03, ...) {
  autoplot.ivate_sensitivity(x, type = type, range = range, ci = ci,
                             n_grid = n_grid, label_tol = label_tol, ...)
}

#' @rdname plot.ivate_sensitivity
#' @param object An `ivate_sensitivity` object.
#' @importFrom ggplot2 autoplot
#' @export
autoplot.ivate_sensitivity <- function(object, type = c("auto", "threshold", "boundary", "comparisons"),
                                       range = NULL, ci = TRUE, n_grid = 101,
                                       label_tol = 0.03, ...) {
  type <- match.arg(type)
  if (type == "comparisons") return(plot_comparisons(object, ...))
  if (type == "auto") type <- if (is.null(object$target)) "boundary" else "threshold"
  if (type == "boundary" && object$setting != "two-sided") {
    stop("The boundary plot requires two-sided non-compliance.", call. = FALSE)
  }
  if (type == "threshold" && is.null(object$target)) {
    stop(paste0("No scalar threshold available. Restrict `theta_A0` or `theta_N1` in ",
                "ate_sensitivity() to obtain one."), call. = FALSE)
  }
  if (is.null(range)) range <- object$y_range
  if (!is.numeric(range) || length(range) != 2L || range[1] >= range[2]) {
    stop("`range` must be a numeric vector c(lower, upper) with lower < upper.", call. = FALSE)
  }
  draws <- if (ci) object$draws else NULL

  if (type == "threshold") {
    plot_threshold(object, range, draws, n_grid, label_tol)
  } else {
    plot_boundary(object, range, draws, n_grid, label_tol)
  }
}

# -- Internal ----------------------------------------------------------------

param_label <- list(
  theta_A0 = quote(theta[A]^"(0)"),
  theta_N1 = quote(theta[N]^"(1)")
)

mean_labels <- c(mu_C0 = "E[Y(0)|C]", mu_C1 = "E[Y(1)|C]", mu_A1 = "E[Y(1)|A]", mu_N0 = "E[Y(0)|N]")

#' Identified potential outcome means used as reference ticks
#' @keywords internal
#' @noRd
reference_means <- function(q, range) {
  vals <- q[names(mean_labels)]
  keep <- !is.na(vals) & vals >= range[1] & vals <= range[2]
  list(vals = unname(vals[keep]), labels = unname(mean_labels[keep]))
}

#' Group nearby secondary-axis labels to reduce overlap
#'
#' Ticks stay at their original positions; within each group of positions
#' closer than `tol * diff(range)`, the combined label is shown at the first
#' position and the others are left blank.
#' @keywords internal
#' @noRd
collapse_nearby_labels <- function(vals, labels, range, tol = 0.03) {
  if (length(vals) == 0) return(list(vals = numeric(0), labels = character(0)))
  ord    <- order(vals)
  vals   <- vals[ord]
  labels <- labels[ord]
  groups <- cumsum(c(TRUE, diff(vals) > tol * diff(range)))
  new_labels <- character(length(vals))
  for (g in unique(groups)) {
    idx <- which(groups == g)
    new_labels[idx[1]] <- paste(labels[idx], collapse = "\n")
  }
  list(vals = vals, labels = new_labels)
}

#' Function returning the least favourable value of the restricted parameter
#' @keywords internal
#' @noRd
fixed_value_fun <- function(x) {
  if (is.null(x$restrictions)) return(function(q) NA_real_)
  fixed <- names(x$restrictions)
  restr <- x$restrictions[[1]]
  function(q) least_favourable(feasible_interval(restr, fixed, q, x$y_range), fixed, x$direction)
}

#' Pointwise percentile band (and median) over draws
#' @keywords internal
#' @noRd
pointwise_band <- function(draws, curve_fun, n_points, level) {
  reps <- matrix(NA_real_, nrow(draws), n_points)
  for (b in seq_len(nrow(draws))) {
    reps[b, ] <- tryCatch(curve_fun(draws[b, ]), error = function(e) rep(NA_real_, n_points))
  }
  alpha <- (1 - level) / 2
  pct <- function(v, p) if (any(is.finite(v))) stats::quantile(v[is.finite(v)], p, names = FALSE) else NA_real_
  list(lower  = apply(reps, 2, pct, p = alpha),
       upper  = apply(reps, 2, pct, p = 1 - alpha),
       median = apply(reps, 2, pct, p = 0.5))
}

#' Central values of the identified quantities used for plotting
#'
#' Plug-in estimates for the bootstrap, posterior medians for Bayesian
#' inference.
#' @keywords internal
#' @noRd
center_q <- function(x) {
  if (!identical(x$inference, "bayes")) return(x$q)
  est <- stats::setNames(x$identified$estimate, x$identified$quantity)
  q <- x$q
  q[intersect(names(q), names(est))] <- est[intersect(names(q), names(est))]
  q
}

secondary_axis <- function(ref) {
  ggplot2::dup_axis(name = NULL, breaks = ref$vals, labels = ref$labels)
}

sensitivity_theme <- function() {
  ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      axis.text.x.top           = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 0, size = 7),
      axis.ticks.x.top          = ggplot2::element_line(colour = "gray30"),
      axis.ticks.length.x.top   = ggplot2::unit(6, "pt"),
      axis.text.y.right         = ggplot2::element_text(size = 7),
      axis.ticks.y.right        = ggplot2::element_line(colour = "gray30"),
      axis.ticks.length.y.right = ggplot2::unit(6, "pt")
    )
}

#' ATE as a function of the target parameter (threshold plot)
#' @keywords internal
#' @noRd
plot_threshold <- function(x, range, draws, n_grid, label_tol) {
  target   <- x$target
  fixed_of <- fixed_value_fun(x)
  grid     <- seq(range[1], range[2], length.out = max(n_grid, 2L) * 3L)

  ate_curve <- function(q) {
    other <- fixed_of(q)
    if (target == "theta_N1") ate_at(q, other, grid) else ate_at(q, grid, other)
  }
  df <- data.frame(param = grid, ate = ate_curve(x$q))
  band <- if (!is.null(draws)) pointwise_band(draws, ate_curve, length(grid), x$level) else NULL
  if (!is.null(band) && identical(x$inference, "bayes")) df$ate <- band$median

  tau <- x$threshold$estimate
  ref <- reference_means(center_q(x), range)
  tau_in_range <- tau >= range[1] && tau <= range[2]
  sec <- collapse_nearby_labels(
    c(ref$vals, if (tau_in_range) tau),
    c(ref$labels, if (tau_in_range) "tau*"),
    range, tol = label_tol
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$param, y = .data$ate))

  ci_visible <- has_ci(x$threshold) &&
    x$threshold$ci_upper >= range[1] && x$threshold$ci_lower <= range[2]
  if (!is.null(band) && ci_visible) {
    p <- p + ggplot2::annotate("rect",
      xmin = max(x$threshold$ci_lower, range[1]), xmax = min(x$threshold$ci_upper, range[2]),
      ymin = -Inf, ymax = Inf, fill = "steelblue", alpha = 0.15)
  }
  if (!is.null(band)) {
    df_band <- data.frame(param = grid, lower = band$lower, upper = band$upper)
    p <- p + ggplot2::geom_ribbon(
      data = df_band, ggplot2::aes(x = .data$param, ymin = .data$lower, ymax = .data$upper),
      inherit.aes = FALSE, fill = "gray30", alpha = 0.3)
  }

  p <- p +
    ggplot2::geom_line(colour = "black", linewidth = 0.7) +
    ggplot2::geom_hline(yintercept = x$cutoff, colour = "black", linewidth = 0.8, alpha = 0.6)
  if (length(ref$vals) > 0) {
    p <- p + ggplot2::geom_vline(xintercept = ref$vals, linetype = "dotted", colour = "gray30")
  }
  if (tau_in_range) {
    p <- p + ggplot2::geom_vline(xintercept = tau, linetype = "dashed", colour = "gray30")
  }

  p +
    ggplot2::scale_x_continuous(name = param_label[[target]], sec.axis = secondary_axis(sec)) +
    ggplot2::scale_y_continuous(name = "ATE") +
    ggplot2::coord_cartesian(xlim = range) +
    sensitivity_theme()
}

#' ATE heatmap with the boundary B (boundary plot)
#' @keywords internal
#' @noRd
plot_boundary <- function(x, range, draws, n_grid, label_tol) {
  seq_r <- seq(range[1], range[2], length.out = n_grid)
  grid <- expand.grid(theta_A0 = seq_r, theta_N1 = seq_r)
  grid$ate <- ate_at(center_q(x), grid$theta_A0, grid$theta_N1)

  line_fun <- function(q) {
    bc <- boundary_coefficients(q, x$cutoff)
    bc[["level_N1_intercept"]] + bc[["level_N1_slope"]] * seq_r
  }
  df_line <- data.frame(theta_A0 = seq_r, theta_N1 = line_fun(x$q))
  band <- if (!is.null(draws)) pointwise_band(draws, line_fun, n_grid, x$level) else NULL
  if (!is.null(band) && identical(x$inference, "bayes")) df_line$theta_N1 <- band$median

  ref <- reference_means(center_q(x), range)
  sec <- collapse_nearby_labels(ref$vals, ref$labels, range, tol = label_tol)

  p <- ggplot2::ggplot(grid, ggplot2::aes(x = .data$theta_A0, y = .data$theta_N1)) +
    ggplot2::geom_raster(ggplot2::aes(fill = .data$ate), interpolate = TRUE, alpha = 0.6)

  if (!is.null(band)) {
    df_band <- data.frame(theta_A0 = seq_r, lower = band$lower, upper = band$upper)
    p <- p + ggplot2::geom_ribbon(
      data = df_band, ggplot2::aes(x = .data$theta_A0, ymin = .data$lower, ymax = .data$upper),
      inherit.aes = FALSE, fill = "gray30", alpha = 0.3)
  }

  p +
    ggplot2::geom_line(data = df_line, colour = "black", linewidth = 0.8) +
    ggplot2::scale_fill_viridis_c(name = "ATE") +
    ggplot2::scale_x_continuous(name = param_label$theta_A0, sec.axis = secondary_axis(sec)) +
    ggplot2::scale_y_continuous(name = param_label$theta_N1, sec.axis = secondary_axis(sec)) +
    ggplot2::coord_cartesian(xlim = range, ylim = range, expand = FALSE) +
    sensitivity_theme()
}
