#' Sensitivity analysis for the ATE under IV assumptions
#'
#' Point identification of the average treatment effect (ATE) is precluded
#' under standard IV assumptions by two unidentified quantities:
#' \eqn{\theta_A = \mu_A^{(0)}}, the mean potential outcome under control for
#' always-takers, and \eqn{\theta_N = \mu_N^{(1)}}, the mean potential outcome
#' under treatment for never-takers. This function derives the thresholds on
#' these sensitivity parameters at which the ATE equals a policy-relevant
#' decision threshold, and produces a corresponding visualisation.
#'
#' The ATE decomposes as (equation (7) in Locher, Richardson, and Stensrud (2026)):
#' \deqn{
#'   \Psi_{ATE} = \pi_C(\mu_C^{(1)} - \mu_C^{(0)})
#'              + \pi_A(\mu_A^{(1)} - \theta_A)
#'              + \pi_N(\theta_N   - \mu_N^{(0)}).
#' }
#'
#' Three compliance scenarios are handled automatically:
#' \describe{
#'   \item{Two-sided non-compliance}{Both always-takers and never-takers are
#'     present. The decision criterion \eqn{\Psi_{ATE} \ge x} defines a
#'     half-plane in \eqn{(\theta_A, \theta_N)} space. The boundary of this
#'     half-plane is returned as a character string and visualised as a
#'     heatmap with a dark red contour at \code{decision_threshold}.}
#'   \item{Always-takers only (\eqn{\pi_N = 0})}{Only \eqn{\theta_A} is
#'     unidentified. A scalar threshold \eqn{\tau_A} is returned together with
#'     the implied additive treatment effect (and, for binary outcomes, the
#'     implied risk ratio) among always-takers. Visualised as a colour strip.}
#'   \item{Never-takers only (\eqn{\pi_A = 0})}{Only \eqn{\theta_N} is
#'     unidentified. A scalar threshold \eqn{\tau_N} is returned together with
#'     the implied additive treatment effect (and, for binary outcomes, the
#'     implied risk ratio) among never-takers. Visualised as a colour strip.}
#' }
#'
#' @param data A data frame containing the outcome, treatment, and instrument
#'   variables.
#' @param formula A two-part formula of the form
#'   \code{outcome ~ treatment | instrument}.
#' @param decision_threshold Numeric scalar. The policy-relevant threshold
#'   \eqn{x} such that the decision criterion is \eqn{\Psi_{ATE} > x}.
#'   Defaults to \code{0}.
#' @param plot_outcome_range Numeric vector of length 2 giving the range of
#'   the outcome used for the sensitivity parameter grid and plot axes. If
#'   \code{NULL} (the default), the observed range of the outcome is used.
#' @param plot_n_grid Integer. Number of grid points along each axis of the
#'   sensitivity parameter grid (two-sided case only). Defaults to \code{101}.
#' @param plot_show_title Logical; whether to display the plot title.
#'   Defaults to \code{TRUE}.
#' @param plot_show_subtitle Logical; whether to display the plot subtitle.
#'   Defaults to \code{TRUE}.
#' @param plot_label_tol Numeric. Tolerance for collapsing nearby
#'   secondary-axis labels, expressed as a fraction of the axis range.
#'   Defaults to \code{0.05}.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{\code{identified_quantities}}{A named list of identified type
#'       shares (\eqn{\pi_C}, \eqn{\pi_A}, \eqn{\pi_N}) and conditional
#'       potential outcome means (\eqn{\mu_C^{(0)}}, \eqn{\mu_C^{(1)}},
#'       \eqn{\mu_A^{(1)}}, \eqn{\mu_N^{(0)}}).}
#'     \item{\code{threshold}}{The sensitivity threshold(s) at which
#'       \eqn{\Psi_{ATE} = } \code{decision_threshold}. In the two-sided
#'       case, a character string giving the boundary equations for both
#'       \eqn{\theta_A} and \eqn{\theta_N}. In the one-sided cases, a named
#'       list containing the scalar threshold (\code{theta_A} or
#'       \code{theta_N}), the \code{implied_effect_additive}, and, for binary
#'       outcomes, the \code{implied_effect_risk_ratio}.}
#'     \item{\code{plot}}{A \code{ggplot} object visualising the sensitivity
#'       analysis.}
#'   }
#'   In the fully compliant case (no non-compliance), \code{NULL} is returned
#'   invisibly and the point-identified ATE is printed as a message.
#'
#' @references
#' Angrist, J. D., Imbens, G. W., and Rubin, D. B. (1996). Identification of
#' causal effects using instrumental variables. \emph{Journal of the American
#' Statistical Association}, 91(434), 444--455.
#' \doi{10.1080/01621459.1996.10476902}
#'
#' Balke, A. and Pearl, J. (1997). Bounds on treatment effects from studies
#' with imperfect compliance. \emph{Journal of the American Statistical
#' Association}, 92(439), 1171--1176.
#' \doi{10.2307/2965583}
#'
#' Locher, L., Richardson, T., and Stensrud, M. (2026). Informative inference on the average
#' treatment effect with instrumental variables.
#'
#' @examples
#' # ── Two-sided non-compliance ──────────────────────────────────────────────
#' set.seed(42)
#' n   <- 2000
#' z   <- rbinom(n, 1, 0.5)
#' a   <- ifelse(z == 1, rbinom(n, 1, 0.7), rbinom(n, 1, 0.2))
#' y   <- rbinom(n, 1, 0.3 + 0.2 * a)
#' dat <- data.frame(y = y, a = a, z = z)
#'
#' res <- ate_sensitivity(dat, y ~ a | z)
#' res$threshold
#' res$plot
#'
#' # ── One-sided non-compliance (never-takers only) ───────────────────────────
#' a_nt <- ifelse(z == 1, rbinom(n, 1, 0.6), 0L)
#' dat2 <- data.frame(y = y, a = a_nt, z = z)
#'
#' res2 <- ate_sensitivity(dat2, y ~ a | z, decision_threshold = 0)
#' res2$threshold
#' res2$plot
#'
#' @export
ate_sensitivity <- function(
    data,
    formula,
    decision_threshold  = 0,
    plot_outcome_range  = NULL,
    plot_n_grid         = 101,
    plot_show_title     = TRUE,
    plot_show_subtitle  = TRUE,
    plot_label_tol      = 0.05
) {
  d <- validate_iv_data(data, formula)

  y <- d$outcome
  a <- d$treatment
  z <- d$instrument

  if (is.null(plot_outcome_range)) plot_outcome_range <- c(min(y), max(y))

  kc     <- get_known_quantities(y, a, z)
  binary <- is_binary(y)
  late   <- kc$y1_c - kc$y0_c

  # Reference lines: all four identified potential outcome means.
  # Y(1) quantities -> dashed; Y(0) quantities -> dotted.
  all_obs        <- c(kc$y0_c, kc$y0_nt, kc$y1_c, kc$y1_at)
  all_obs_labels <- c("E[Y(0)|C]", "E[Y(0)|N]", "E[Y(1)|C]", "E[Y(1)|A]")
  all_obs_lty    <- c("dotted",    "dotted",     "dashed",    "dashed")

  # ── Two-sided non-compliance ──────────────────────────────────────────────
  if (kc$at_exist && kc$nt_exist) {

    y0_at_seq <- seq(plot_outcome_range[1], plot_outcome_range[2], length.out = plot_n_grid)
    y1_nt_seq <- seq(plot_outcome_range[1], plot_outcome_range[2], length.out = plot_n_grid)

    grid     <- expand.grid(y0_at = y0_at_seq, y1_nt = y1_nt_seq)
    grid$ate <- kc$p_c * late +
      kc$p_at * (kc$y1_at - grid$y0_at) +
      kc$p_nt * (grid$y1_nt - kc$y0_nt)

    # Boundary expressed separately for each sensitivity parameter
    # (equations (8) and (9) in Locher, Richardson, and Stensrud (2026)).
    a_coef_A <- (kc$p_at * kc$y1_at + kc$p_c * late - kc$p_nt * kc$y0_nt - decision_threshold) / kc$p_at
    b_coef_A <- kc$p_nt / kc$p_at
    str_A    <- paste0("theta_A = ", fmt(a_coef_A), " + ", fmt(b_coef_A), " * theta_N")

    a_coef_N <- (decision_threshold - kc$p_c * late - kc$p_at * kc$y1_at + kc$p_nt * kc$y0_nt) / kc$p_nt
    b_coef_N <- kc$p_at / kc$p_nt
    str_N    <- paste0("theta_N = ", fmt(a_coef_N), " + ", fmt(b_coef_N), " * theta_A")

    plot_out <- plot_sensitivity_two_sided(
      grid           = grid,
      x_var          = "y0_at",
      x_label        = expression(theta[A]),
      y_var          = "y1_nt",
      y_label        = expression(theta[N]),
      threshold      = decision_threshold,
      value_range    = c(plot_outcome_range, plot_outcome_range),
      vline_x        = all_obs,
      vline_labels   = all_obs_labels,
      vline_lty      = all_obs_lty,
      hline_y        = all_obs,
      hline_labels   = all_obs_labels,
      hline_lty      = all_obs_lty,
      show_title     = plot_show_title,
      show_subtitle  = plot_show_subtitle,
      label_tol      = plot_label_tol
    )

    return(list(
      identified_quantities = kc,
      threshold             = paste0(str_A, ", ", str_N),
      plot                  = plot_out
    ))
  }

  # ── Always-takers only ────────────────────────────────────────────────────
  if (kc$at_exist) {

    critical_y0_at    <- (kc$p_c * late + kc$p_at * kc$y1_at - decision_threshold) / kc$p_at
    implied_effect_at <- kc$y1_at - critical_y0_at

    threshold_out <- list(
      theta_A                 = critical_y0_at,
      implied_effect_additive = implied_effect_at
    )
    if (binary && critical_y0_at != 0) {
      threshold_out$implied_effect_risk_ratio <- kc$y1_at / critical_y0_at
    }

    plot_out <- plot_sensitivity_one_sided(
      x_label        = expression(theta[A]),
      target_val     = critical_y0_at,
      known_vals     = c(kc$y0_c, kc$y1_at, kc$y1_c),
      known_labels   = c("E[Y(0)|C]", "E[Y(1)|A]", "E[Y(1)|C]"),
      known_lty      = c("dotted",    "dashed",    "dashed"),
      threshold      = decision_threshold,
      value_range    = plot_outcome_range,
      ate_fun        = function(y0_at) kc$p_c * late + kc$p_at * (kc$y1_at - y0_at),
      type_label     = "Always-Takers",
      show_title     = plot_show_title,
      show_subtitle  = plot_show_subtitle,
      label_tol      = plot_label_tol
    )

    return(list(
      identified_quantities = kc,
      threshold             = threshold_out,
      plot                  = plot_out
    ))
  }

  # ── Never-takers only ─────────────────────────────────────────────────────
  if (kc$nt_exist) {

    critical_y1_nt    <- (decision_threshold - kc$p_c * late + kc$p_nt * kc$y0_nt) / kc$p_nt
    implied_effect_nt <- critical_y1_nt - kc$y0_nt

    threshold_out <- list(
      theta_N                 = critical_y1_nt,
      implied_effect_additive = implied_effect_nt
    )
    if (binary && kc$y0_nt != 0) {
      threshold_out$implied_effect_risk_ratio <- critical_y1_nt / kc$y0_nt
    }

    plot_out <- plot_sensitivity_one_sided(
      x_label        = expression(theta[N]),
      target_val     = critical_y1_nt,
      known_vals     = c(kc$y0_c, kc$y1_c, kc$y0_nt),
      known_labels   = c("E[Y(0)|C]", "E[Y(1)|C]", "E[Y(0)|N]"),
      known_lty      = c("dotted",    "dashed",    "dotted"),
      threshold      = decision_threshold,
      value_range    = plot_outcome_range,
      ate_fun        = function(y1_nt) kc$p_c * late + kc$p_nt * (y1_nt - kc$y0_nt),
      type_label     = "Never-Takers",
      show_title     = plot_show_title,
      show_subtitle  = plot_show_subtitle,
      label_tol      = plot_label_tol
    )

    return(list(
      identified_quantities = kc,
      threshold             = threshold_out,
      plot                  = plot_out
    ))
  }

  # ── Full compliance: ATE is point identified ──────────────────────────────
  ate <- kc$p_c * late
  message("No non-compliance observed in data. ATE is point identified: ", ate)
  invisible(NULL)
}
