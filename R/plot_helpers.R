# Internal plotting helpers — not exported.

#' Sensitivity heatmap for two-sided non-compliance
#'
#' Produces a heatmap of the ATE over the two-dimensional sensitivity parameter
#' space \eqn{(\theta_A, \theta_N)}, with a contour at \code{threshold}
#' (shown in dark red) and reference lines for all identified potential outcome
#' means. Y(1) reference lines are dashed; Y(0) reference lines are dotted.
#'
#' @param grid          A data frame with columns named by \code{x_var},
#'   \code{y_var}, and \code{z_var}.
#' @param x_var         Name of the column to use as the x-axis variable.
#' @param x_label       Label for the x-axis (may be an \code{expression}).
#' @param y_var         Name of the column to use as the y-axis variable.
#' @param y_label       Label for the y-axis (may be an \code{expression}).
#' @param z_var         Name of the column to use as the fill variable
#'   (default \code{"ate"}).
#' @param threshold     The ATE value at which to draw the contour line.
#' @param value_range   Numeric vector of length 4:
#'   \code{c(x_min, x_max, y_min, y_max)}.
#' @param vline_x       Numeric vector of x-positions for vertical reference
#'   lines.
#' @param vline_labels  Character vector of labels for \code{vline_x}.
#' @param vline_lty     Character vector of linetypes for \code{vline_x}.
#' @param hline_y       Numeric vector of y-positions for horizontal reference
#'   lines.
#' @param hline_labels  Character vector of labels for \code{hline_y}.
#' @param hline_lty     Character vector of linetypes for \code{hline_y}.
#' @param show_title    Logical; whether to display the plot title.
#' @param show_subtitle Logical; whether to display the plot subtitle.
#' @param label_tol     Tolerance for collapsing nearby secondary-axis labels
#'   (as a fraction of the axis range).
#' @return A \code{ggplot} object.
#' @keywords internal
plot_sensitivity_two_sided <- function(
    grid,
    x_var, x_label,
    y_var, y_label,
    z_var          = "ate",
    threshold,
    value_range,
    vline_x,  vline_labels,  vline_lty,
    hline_y,  hline_labels,  hline_lty,
    show_title     = TRUE,
    show_subtitle  = TRUE,
    label_tol      = 0.05
) {
  stopifnot(
    length(vline_x) == length(vline_labels),
    length(vline_x) == length(vline_lty),
    length(hline_y) == length(hline_labels),
    length(hline_y) == length(hline_lty)
  )

  x_range   <- value_range[1:2]
  y_range   <- value_range[3:4]
  cx        <- collapse_nearby_labels(vline_x, vline_labels, x_range, tol = label_tol)
  cy        <- collapse_nearby_labels(hline_y, hline_labels, y_range, tol = label_tol)

  df_vlines <- data.frame(x = vline_x, lty = vline_lty)
  df_hlines <- data.frame(y = hline_y, lty = hline_lty)

  p <- ggplot2::ggplot(
    grid,
    ggplot2::aes(
      x = .data[[x_var]],
      y = .data[[y_var]],
      z = .data[[z_var]]
    )
  ) +
    ggplot2::geom_raster(
      ggplot2::aes(fill = .data[[z_var]]),
      interpolate = TRUE,
      alpha       = 0.6
    ) +
    ggplot2::geom_contour(colour = "darkred", breaks = threshold, linewidth = 1.2) +
    ggplot2::geom_vline(
      data        = df_vlines,
      ggplot2::aes(xintercept = x, linetype = lty),
      colour      = "gray30",
      show.legend = FALSE
    ) +
    ggplot2::geom_hline(
      data        = df_hlines,
      ggplot2::aes(yintercept = y, linetype = lty),
      colour      = "gray30",
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_viridis_c(name = "ATE") +
    ggplot2::scale_x_continuous(
      name   = x_label,
      breaks = seq(x_range[1], x_range[2], length.out = 6),
      labels = function(x) format(round(x, 2)),
      sec.axis = ggplot2::dup_axis(
        name   = NULL,
        breaks = cx$vals,
        labels = cx$labels
      )
    ) +
    ggplot2::scale_y_continuous(
      name   = y_label,
      breaks = seq(y_range[1], y_range[2], length.out = 6),
      labels = function(x) format(round(x, 2)),
      sec.axis = ggplot2::dup_axis(
        name   = NULL,
        breaks = cy$vals,
        labels = cy$labels
      )
    ) +
    ggplot2::coord_cartesian(xlim = x_range, ylim = y_range) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      axis.text.x.top   = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 0, size = 7),
      axis.text.y.right = ggplot2::element_text(size = 7)
    )

  if (show_title)    p <- p + ggplot2::labs(title = "Sensitivity analysis")
  if (show_subtitle) p <- p + ggplot2::labs(subtitle = paste("Red contour: ATE =", threshold))

  p
}

#' Sensitivity colour strip for one-sided non-compliance
#'
#' Produces a horizontal colour strip showing the ATE as a function of a
#' single sensitivity parameter, with a vertical line in dark red at the
#' critical threshold value.
#'
#' @param x_label       Label for the x-axis (may be an \code{expression}).
#' @param target_val    The critical value of the sensitivity parameter at
#'   which \eqn{\Psi_{ATE} = } \code{threshold} (drawn as dark red line).
#' @param known_vals    Numeric vector of x-positions for reference lines.
#' @param known_labels  Character vector of labels for \code{known_vals}.
#' @param known_lty     Character vector of linetypes for \code{known_vals}.
#' @param threshold     The ATE value used in the subtitle.
#' @param value_range   Numeric vector of length 2 giving the x-axis range.
#' @param ate_fun       A function taking a numeric vector of sensitivity
#'   parameter values and returning the corresponding ATE values.
#' @param type_label    Short string appended to the plot title (e.g.,
#'   \code{"Never-Takers"}).
#' @param show_title    Logical; whether to display the plot title.
#' @param show_subtitle Logical; whether to display the plot subtitle.
#' @param label_tol     Tolerance for collapsing nearby secondary-axis labels
#'   (as a fraction of the axis range).
#' @return A \code{ggplot} object.
#' @keywords internal
plot_sensitivity_one_sided <- function(
    x_label,
    target_val,
    known_vals,
    known_labels,
    known_lty,
    threshold,
    value_range,
    ate_fun,
    type_label     = "",
    show_title     = TRUE,
    show_subtitle  = TRUE,
    label_tol      = 0.05
) {
  df_grid     <- data.frame(param = seq(value_range[1], value_range[2], length.out = 300))
  df_grid$ate <- ate_fun(df_grid$param)
  df_known    <- data.frame(x = known_vals, lty = known_lty)

  # Include target_val in collapsing so nearby known labels absorb it without
  # adding stray text (the target_val label is intentionally left blank).
  all_vals   <- c(known_vals, target_val)
  all_labels <- c(known_labels, "")
  cx         <- collapse_nearby_labels(all_vals, all_labels, value_range, tol = label_tol)

  p <- ggplot2::ggplot(df_grid, ggplot2::aes(x = param)) +
    ggplot2::geom_tile(ggplot2::aes(y = 0, fill = ate), height = 1, alpha = 0.5) +
    ggplot2::geom_vline(xintercept = target_val, color = "darkred", linewidth = 1.2) +
    ggplot2::geom_vline(
      data        = df_known,
      ggplot2::aes(xintercept = x, linetype = lty),
      colour      = "gray30",
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_viridis_c(name = "ATE") +
    ggplot2::scale_x_continuous(
      name     = x_label,
      limits   = value_range,
      sec.axis = ggplot2::dup_axis(
        breaks = cx$vals,
        labels = cx$labels,
        name   = ""
      )
    ) +
    ggplot2::labs(y = NULL) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      axis.text.y        = ggplot2::element_blank(),
      axis.ticks.y       = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      axis.text.x.top    = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 0, size = 7)
    )

  if (show_title)    p <- p + ggplot2::labs(title = paste("Sensitivity analysis:", type_label))
  if (show_subtitle) p <- p + ggplot2::labs(subtitle = paste("Red line: ATE =", threshold))

  p
}
