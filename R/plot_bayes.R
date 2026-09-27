#' Plot a Bayesian analysis of the ATE
#'
#' Visualises the result of [ate_bayes()].
#'
#' * **Sensitivity analysis over one parameter**: the posterior median of the
#'   ATE as a function of the sensitivity parameter, with pointwise credible
#'   band (shaded) and simultaneous credible band (dashed), the cutoff, and
#'   the posterior medians of the identified bounds (dotted); cf. Figure 5 of
#'   Richardson, Evans, and Robins (2011). With `type = "probability"`, the
#'   posterior probability that the decision criterion holds is shown
#'   instead.
#' * **Sensitivity analysis over both parameters**: a heatmap of the
#'   posterior probability that the criterion holds (default) or of the
#'   posterior median ATE (`type = "ate"`), with contours at probability 0.5
#'   or at the cutoff, respectively.
#' * **Full Bayesian analysis** (or a single value): the posterior density of
#'   the ATE, optionally with its prior density, the cutoff, and the posterior
#'   medians of the identified bounds; cf. Figure 3 of Richardson, Evans, and
#'   Robins (2011).
#'
#' @param x,object An `ivate_bayes` object.
#' @param type `"auto"` (default), `"ate"`, or `"probability"`; see Details.
#'   Ignored for density plots.
#' @param show_prior Logical; overlay the prior density of the ATE in a full
#'   Bayesian analysis? Defaults to `TRUE`.
#' @param ... Unused.
#'
#' @return A `ggplot` object.
#' @export
plot.ivate_bayes <- function(x, type = c("auto", "ate", "probability"), show_prior = TRUE, ...) {
  autoplot.ivate_bayes(x, type = type, show_prior = show_prior, ...)
}

#' @rdname plot.ivate_bayes
#' @export
autoplot.ivate_bayes <- function(object, type = c("auto", "ate", "probability"),
                                 show_prior = TRUE, ...) {
  type <- match.arg(type)
  if (object$mode == "full") return(plot_bayes_density(object, show_prior))

  g <- object$grid
  varying <- c("theta_A0", "theta_N1")[c(length(unique(stats::na.omit(g$theta_A0))) > 1,
                                        length(unique(stats::na.omit(g$theta_N1))) > 1)]
  if (length(varying) == 0) return(plot_bayes_density(object, FALSE))
  if (length(varying) == 1) {
    plot_bayes_curve(object, varying, if (type == "auto") "ate" else type)
  } else {
    plot_bayes_surface(object, if (type == "auto") "probability" else type)
  }
}

prob_axis_label <- function(x) {
  sprintf("P(ATE %s %s | data)", if (x$direction == "greater") ">" else "<", fmt_num(x$cutoff))
}

bayes_theme <- function() ggplot2::theme_minimal(base_size = 13)

#' Fixed values of the non-varying parameter, as a subtitle
#' @keywords internal
#' @noRd
fixed_subtitle <- function(x, varying) {
  other <- setdiff(c("theta_A0", "theta_N1"), varying)
  v <- unique(stats::na.omit(x$grid[[other]]))
  if (length(v) == 1) paste(other, "=", fmt_num(v)) else NULL
}

#' Posterior of the ATE along one sensitivity parameter
#' @keywords internal
#' @noRd
plot_bayes_curve <- function(x, param, type) {
  g <- x$grid
  g$param <- g[[param]]
  p <- ggplot2::ggplot(g, ggplot2::aes(x = .data$param))
  if (type == "ate") {
    bounds <- x$bounds$estimate
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
                           fill = "gray30", alpha = 0.3) +
      ggplot2::geom_line(ggplot2::aes(y = .data$band_lower), linetype = "dashed", colour = "gray30") +
      ggplot2::geom_line(ggplot2::aes(y = .data$band_upper), linetype = "dashed", colour = "gray30") +
      ggplot2::geom_line(ggplot2::aes(y = .data$median), colour = "black", linewidth = 0.8) +
      ggplot2::geom_hline(yintercept = x$cutoff, colour = "black", linewidth = 0.8, alpha = 0.6) +
      ggplot2::geom_hline(yintercept = bounds, linetype = "dotted", colour = "gray30") +
      # Zoom on the credible band; bound lines far outside it would flatten the curve.
      ggplot2::coord_cartesian(ylim = range(c(g$band_lower, g$band_upper, x$cutoff))) +
      ggplot2::labs(y = "ATE")
  } else {
    p <- p +
      ggplot2::geom_hline(yintercept = 0.5, linetype = "dotted", colour = "gray30") +
      ggplot2::geom_line(ggplot2::aes(y = .data$prob), colour = "black", linewidth = 0.8) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(y = prob_axis_label(x))
  }
  p + ggplot2::labs(x = param_label[[param]], subtitle = fixed_subtitle(x, param)) + bayes_theme()
}

#' Posterior probability (or median ATE) over both sensitivity parameters
#' @keywords internal
#' @noRd
plot_bayes_surface <- function(x, type) {
  g <- x$grid
  g$fill <- if (type == "probability") g$prob else g$median
  p <- ggplot2::ggplot(g, ggplot2::aes(x = .data$theta_A0, y = .data$theta_N1)) +
    ggplot2::geom_raster(ggplot2::aes(fill = .data$fill), interpolate = TRUE, alpha = 0.8)
  if (type == "probability") {
    p <- p +
      ggplot2::geom_contour(ggplot2::aes(z = .data$prob), breaks = 0.5, colour = "black", linewidth = 0.8) +
      ggplot2::scale_fill_viridis_c(name = prob_axis_label(x), limits = c(0, 1))
  } else {
    p <- p +
      ggplot2::geom_contour(ggplot2::aes(z = .data$median), breaks = x$cutoff, colour = "black", linewidth = 0.8) +
      ggplot2::scale_fill_viridis_c(name = "Posterior\nmedian ATE")
  }
  p +
    ggplot2::labs(x = param_label$theta_A0, y = param_label$theta_N1) +
    ggplot2::coord_cartesian(expand = FALSE) +
    bayes_theme()
}

#' Posterior (and prior) density of the ATE
#' @keywords internal
#' @noRd
plot_bayes_density <- function(x, show_prior) {
  df <- data.frame(ate = x$draws$ate, distribution = "Posterior")
  if (show_prior && !is.null(x$prior_ate)) {
    df <- rbind(df, data.frame(ate = x$prior_ate, distribution = "Prior"))
  }
  df$distribution <- factor(df$distribution, levels = c("Posterior", "Prior"))
  ggplot2::ggplot(df, ggplot2::aes(x = .data$ate, linetype = .data$distribution)) +
    ggplot2::geom_vline(xintercept = x$bounds$estimate, linetype = "dotted", colour = "gray30") +
    ggplot2::geom_vline(xintercept = x$cutoff, colour = "black", linewidth = 0.8, alpha = 0.6) +
    ggplot2::geom_density(colour = "black", linewidth = 0.7, key_glyph = "path") +
    ggplot2::scale_linetype_manual(values = c(Posterior = "solid", Prior = "dashed"), name = NULL,
                                   guide = if (nlevels(droplevels(df$distribution)) > 1) "legend" else "none") +
    ggplot2::labs(x = "ATE", y = "Density",
                  subtitle = if (x$mode == "sensitivity") fixed_subtitle_all(x) else NULL) +
    bayes_theme()
}

fixed_subtitle_all <- function(x) {
  g <- x$grid
  parts <- c(if (!is.na(g$theta_A0[1])) paste("theta_A0 =", fmt_num(g$theta_A0[1])),
             if (!is.na(g$theta_N1[1])) paste("theta_N1 =", fmt_num(g$theta_N1[1])))
  paste(parts, collapse = ", ")
}
