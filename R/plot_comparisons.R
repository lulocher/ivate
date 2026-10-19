#' Plot the roster of identified, assumed, and implied quantities
#'
#' Draws the potential outcome means by compliance type and treatment, in the
#' layout of Figure 1 of Locher, Richardson, and Stensrud (2026), with the
#' values of a sensitivity analysis filled in (as in the figures of its
#' examples). Rows are \eqn{E[Y(1)]} and
#' \eqn{E[Y(0)]}; columns are compliers, always-takers, and never-takers
#' (types absent from the data are dropped). Each position shows the value
#' of the corresponding quantity; the unidentified sensitivity parameters sit
#' at \eqn{E[Y(0)]} for always-takers and \eqn{E[Y(1)]} for never-takers.
#' Vertical arrows indicate within-type comparisons and curved arrows
#' cross-type comparisons; arrows are dashed when they involve an
#' unidentified quantity.
#'
#' Quantities are coloured by their status:
#' * **estimated** (green): identified potential outcome means;
#' * **assumed** (orange): a sensitivity parameter fixed by the user or set
#'   to its least favourable value under a restriction (step 1 or 2 of the
#'   two-step strategy);
#' * **implied** (red): the threshold \eqn{\tau} of the remaining sensitivity
#'   parameter, shown with the relation under which the decision criterion
#'   holds (e.g. "> 0.46"). Without a restriction in the two-sided case, both sensitivity
#'   parameters are implied and are shown with the estimated boundary
#'   \eqn{\mathcal{B}}.
#'
#' Arrows take the colour of the "strongest" status among their two ends
#' (implied over assumed over estimated). For Bayesian results, values are
#' posterior medians.
#'
#' @param x An `ivate_sensitivity` object from [ate_sensitivity()].
#' @param ci Logical; show confidence (credible) intervals below the values?
#'   Defaults to `FALSE`.
#' @param digits Number of significant digits shown for all values (e.g.
#'   0.70, 0.011, 0.00084 with the default `2`).
#' @param colours Named character vector with colours for `estimated`,
#'   `assumed`, and `implied` quantities.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' set.seed(1)
#' n <- 2000
#' z <- rbinom(n, 1, 0.5)
#' a <- ifelse(z == 1, rbinom(n, 1, 0.7), rbinom(n, 1, 0.2))
#' y <- rbinom(n, 1, 0.3 + 0.2 * a)
#' dat <- data.frame(y, a, z)
#'
#' # Two-sided, unrestricted: both sensitivity parameters on the boundary
#' plot_comparisons(ate_sensitivity(dat, y ~ a | z, ci = FALSE))
#'
#' # Step 1 of the two-step strategy: theta_A0 fixed at its worst case
#' res <- ate_sensitivity(dat, y ~ a | z, theta_A0 = 1, n_boot = 200)
#' plot_comparisons(res, ci = TRUE)
#'
#' @export
plot_comparisons <- function(x, ci = FALSE, digits = 2,
                             colours = c(estimated = "#2E7D32", assumed = "#E08E0B",
                                         implied = "#C62828")) {
  if (!inherits(x, "ivate_sensitivity")) {
    stop("`x` must be the result of ate_sensitivity().", call. = FALSE)
  }
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits < 1) {
    stop("`digits` must be a single number >= 1 (significant digits).", call. = FALSE)
  }
  if (!all(c("estimated", "assumed", "implied") %in% names(colours))) {
    stop("`colours` must be named `estimated`, `assumed`, and `implied`.", call. = FALSE)
  }

  nodes  <- comparison_nodes(x, ci, digits)
  arrows <- comparison_arrows(nodes)
  status_levels <- c("estimated", "assumed", "implied")
  nodes$status  <- factor(nodes$status, levels = status_levels)
  arrows$status <- factor(arrows$status, levels = status_levels)

  cols <- sort(unique(nodes$x))
  headers <- data.frame(x = cols, y = max(nodes$y) + 1.05,
                        label = type_headers[nodes$type[match(cols, nodes$x)]])
  rows <- data.frame(x = min(cols) - 0.75, y = c(row_y[["1"]], row_y[["0"]]),
                     label = c('E*"["*Y*"(1)]"', 'E*"["*Y*"(0)]"'))

  p <- ggplot2::ggplot() +
    ggplot2::geom_text(data = headers, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
                       fontface = "bold", size = 4.5) +
    ggplot2::geom_text(data = rows, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
                       parse = TRUE, size = 4.5)

  arrow_head <- grid::arrow(length = grid::unit(2.5, "mm"), ends = "both", type = "closed")
  for (i in seq_len(nrow(arrows))) {
    a <- arrows[i, ]
    layer_data <- a
    if (a$curvature == 0) {
      p <- p + ggplot2::geom_segment(data = layer_data,
        ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
                     colour = .data$status),
        linetype = a$linetype, arrow = arrow_head, linewidth = 0.5, show.legend = FALSE)
    } else {
      p <- p + ggplot2::geom_curve(data = layer_data,
        ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
                     colour = .data$status),
        curvature = a$curvature, linetype = a$linetype, arrow = arrow_head,
        linewidth = 0.5, show.legend = FALSE)
    }
  }

  # Values (and intervals below them); the boundary formula is drawn as plotmath
  value_dy    <- if (ci) 0.05 else 0
  value_plain <- nodes[!nodes$value_parse, , drop = FALSE]
  value_math  <- nodes[nodes$value_parse, , drop = FALSE]
  p <- p +
    ggplot2::geom_text(data = value_plain,
      ggplot2::aes(x = .data$x, y = .data$y + value_dy, label = .data$value, colour = .data$status),
      size = 4.3, show.legend = FALSE)
  if (nrow(value_math) > 0) {
    p <- p + ggplot2::geom_text(data = value_math,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$value, colour = .data$status),
      parse = TRUE, size = 3.8, show.legend = FALSE)
  }
  if (ci) {
    ci_nodes <- nodes[nodes$ci_label != "", , drop = FALSE]
    if (nrow(ci_nodes) > 0) {
      p <- p + ggplot2::geom_text(data = ci_nodes,
        ggplot2::aes(x = .data$x, y = .data$y - 0.08, label = .data$ci_label, colour = .data$status),
        size = 3.1, show.legend = FALSE)
    }
  }

  # Invisible points to produce a legend with the three statuses
  present <- status_levels[status_levels %in% c(as.character(nodes$status), as.character(arrows$status))]
  legend_df <- data.frame(x = mean(cols), y = mean(nodes$y),
                          status = factor(present, levels = status_levels))
  p +
    ggplot2::geom_point(data = legend_df,
      ggplot2::aes(x = .data$x, y = .data$y, colour = .data$status), alpha = 0, size = 3) +
    ggplot2::scale_colour_manual(values = colours, limits = present, name = NULL,
      labels = c(estimated = "Estimated", assumed = "Assumed", implied = "Implied")) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1, shape = 15, size = 4))) +
    ggplot2::labs(subtitle = comparison_subtitle(x)) +
    ggplot2::coord_cartesian(xlim = c(min(cols) - 1.05, max(cols) + 0.45),
                             ylim = c(min(nodes$y) - 0.95, max(nodes$y) + 1.2)) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", plot.margin = ggplot2::margin(8, 8, 8, 8))
}

# -- Internal ------------------------------------------------------------------

node_ry <- 0.17
row_y   <- c("1" = 2.2, "0" = 1)
type_headers <- c(C = "Compliers", A = "Always-takers", N = "Never-takers")

#' Nodes of the comparison roster with values and status
#' @keywords internal
#' @noRd
comparison_nodes <- function(x, ci, digits = 2) {
  has_A <- x$q[["pi_A"]] > 0
  has_N <- x$q[["pi_N"]] > 0
  types <- c("C", if (has_A) "A", if (has_N) "N")
  col_x <- stats::setNames(seq(1, by = 1.25, length.out = length(types)), types)

  id <- x$identified
  est <- function(qty) {
    r <- id[id$quantity == qty, ]
    # Complier means are undefined without compliers (pi_C = 0)
    if (nrow(r) == 0) return(list(value = "undefined", ci = ""))
    list(value = fmt_digits(r$estimate, digits),
         ci = if (has_ci(r)) fmt_paren_ci(r$ci_lower, r$ci_upper, digits) else "")
  }
  node <- function(key, type, arm, name, status, value, ci_label = "", circled = FALSE,
                   value_parse = FALSE) {
    data.frame(id = key, type = type, x = col_x[[type]], y = row_y[[arm]], name = name,
               status = status, value = value, ci_label = ci_label, circled = circled,
               value_parse = value_parse, stringsAsFactors = FALSE)
  }

  out <- list(
    node("C1", "C", "1", 'mu[C]^"(1)"', "estimated", est("mu_C1")$value, est("mu_C1")$ci),
    node("C0", "C", "0", 'mu[C]^"(0)"', "estimated", est("mu_C0")$value, est("mu_C0")$ci)
  )
  if (has_A) {
    out <- c(out, list(
      node("A1", "A", "1", 'mu[A]^"(1)"', "estimated", est("mu_A1")$value, est("mu_A1")$ci),
      theta_node(x, "theta_A0", "A0", "A", "0", 'theta[A]^"(0)"', node, digits)
    ))
  }
  if (has_N) {
    out <- c(out, list(
      theta_node(x, "theta_N1", "N1", "N", "1", 'theta[N]^"(1)"', node, digits),
      node("N0", "N", "0", 'mu[N]^"(0)"', "estimated", est("mu_N0")$value, est("mu_N0")$ci)
    ))
  }
  do.call(rbind, out)
}

#' Node for a sensitivity parameter: implied, assumed, or on the boundary
#' @keywords internal
#' @noRd
theta_node <- function(x, parameter, key, type, arm, name, node, digits = 2) {
  if (identical(x$target, parameter)) {
    t <- x$threshold
    ci_label <- if (has_ci(t)) fmt_paren_ci(t$ci_lower, t$ci_upper, digits) else ""
    return(node(key, type, arm, name, "implied", paste(t$relation, fmt_digits(t$estimate, digits)),
                ci_label, circled = TRUE))
  }
  if (!is.null(x$fixed) && identical(x$fixed$parameter, parameter)) {
    f <- x$fixed
    ci_label <- if (has_ci(f) && !isTRUE(f$ci_lower == f$ci_upper)) fmt_paren_ci(f$ci_lower, f$ci_upper, digits) else ""
    return(node(key, type, arm, name, "assumed", fmt_digits(f$value, digits), ci_label, circled = TRUE))
  }
  # Two-sided without restriction: show the estimated boundary.
  b <- x$boundary[x$boundary$scale == "level" & x$boundary$lhs == parameter, ]
  other <- if (parameter == "theta_A0") 'theta[A]^"(0)"' else 'theta[N]^"(1)"'
  rhs <- if (b$rhs == "theta_A0") 'theta[A]^"(0)"' else 'theta[N]^"(1)"'
  sign <- if (b$slope < 0) "-" else "+"
  formula <- sprintf('"%s " * "%s" %s "%s" %%.%% %s', b$relation, fmt_digits(b$intercept, digits), sign,
                     fmt_digits(abs(b$slope), digits), rhs)
  node(key, type, arm, name, "implied", formula, circled = TRUE, value_parse = TRUE)
}

#' Arrows of the comparison roster (as in Figure 1 of the manuscript)
#' @keywords internal
#' @noRd
comparison_arrows <- function(nodes) {
  rank <- c(estimated = 1, assumed = 2, implied = 3)
  get <- function(id) nodes[nodes$id == id, ]
  arrow <- function(from, to, kind) {
    a <- get(from); b <- get(to)
    if (nrow(a) == 0 || nrow(b) == 0) return(NULL)
    status <- names(rank)[max(rank[a$status], rank[b$status])]
    linetype <- if (a$circled || b$circled) "dashed" else "solid"
    top <- a$y == row_y[["1"]]
    if (kind == "vertical") {
      data.frame(x = a$x, y = a$y + node_ry + 0.03, xend = b$x, yend = b$y - node_ry - 0.03,
                 curvature = 0, linetype = linetype, status = status)
    } else {
      long <- abs(b$x - a$x) > 1.3
      dy <- (node_ry + 0.03) * if (top) 1 else -1
      inset <- if (long) 0 else 0.12
      curv <- (if (long) 0.35 else 0.45) * if (top) -1 else 1
      data.frame(x = a$x + inset, y = a$y + dy, xend = b$x - inset, yend = b$y + dy,
                 curvature = curv, linetype = linetype, status = status)
    }
  }
  do.call(rbind, list(
    arrow("C0", "C1", "vertical"), arrow("A0", "A1", "vertical"), arrow("N0", "N1", "vertical"),
    arrow("C1", "A1", "cross"), arrow("A1", "N1", "cross"),
    arrow("C1", "N1", "cross"),
    arrow("C0", "A0", "cross"), arrow("A0", "N0", "cross"),
    arrow("C0", "N0", "cross")
  ))
}

#' Subtitle describing the criterion and the restriction
#' @keywords internal
#' @noRd
comparison_subtitle <- function(x) {
  crit <- sprintf("Criterion: ATE %s %s", if (x$direction == "greater") ">" else "<", fmt_num(x$cutoff))
  if (!is.null(x$fixed)) {
    f <- x$fixed
    fixed_number <- grepl(paste0("^", f$parameter, " = "), f$restriction)
    crit <- paste0(crit, "; ", if (fixed_number) {
      paste("assumed", f$restriction)
    } else {
      paste0(f$parameter, " at its least favourable value under ", f$restriction)
    })
  }
  crit
}

#' Format an interval as "(lower, upper)"
#' @keywords internal
#' @noRd
fmt_paren_ci <- function(lower, upper, digits = 2) {
  paste0("(", fmt_digits(lower, digits), ", ", fmt_digits(upper, digits), ")")
}

#' Format with `digits` significant digits, keeping trailing zeros (0.70)
#' @keywords internal
#' @noRd
fmt_digits <- function(x, digits = 2) {
  vapply(x, function(v) {
    if (is.na(v)) return("NA")
    out <- trimws(formatC(signif(v, digits), format = "fg", digits = digits, flag = "#"))
    sub("\\.$", "", out)
  }, character(1), USE.NAMES = FALSE)
}
