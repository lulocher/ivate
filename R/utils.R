# Internal helper functions — not exported.
# These are not intended to be called directly by users.

#' Check whether a vector is binary (0/1)
#' @keywords internal
is_binary <- function(y) all(y %in% c(0, 1))

#' Round a number to a given number of digits
#' @keywords internal
fmt <- function(x, digits = 3) round(x, digits)

#' Collapse nearby secondary-axis labels
#'
#' Values within \code{tol * diff(range)} of each other are grouped; each
#' group's break is placed at the group mean and labels are joined with
#' \code{"\\n"}.
#'
#' @param vals   Numeric vector of axis positions.
#' @param labels Character vector of axis labels (same length as \code{vals}).
#' @param range  Numeric vector of length 2 giving the axis range.
#' @param tol    Tolerance as a fraction of the axis range (default 0.05).
#' @return A named list with elements \code{vals} and \code{labels}.
#' @keywords internal
collapse_nearby_labels <- function(vals, labels, range, tol = 0.05) {
  tol_abs <- tol * diff(range)
  ord     <- order(vals)
  vals    <- vals[ord]
  labels  <- labels[ord]

  groups     <- cumsum(c(TRUE, diff(vals) > tol_abs))
  new_vals   <- as.numeric(tapply(vals,   groups, mean))
  new_labels <- as.character(tapply(labels, groups, paste, collapse = "\n"))

  list(vals = new_vals, labels = new_labels)
}

#' Extract outcome, treatment, and instrument vectors from a formula
#'
#' Parses a formula of the form \code{outcome ~ treatment | instrument} and
#' returns the corresponding columns from \code{data}.
#'
#' @param data    A data frame.
#' @param formula A formula with a \code{|} separator for the instrument, e.g.
#'   \code{y ~ a | z}.
#' @return A named list with elements \code{outcome}, \code{treatment}, and
#'   \code{instrument}.
#' @keywords internal
get_data_from_formula <- function(data, formula) {
  vars <- all.vars(as.formula(formula))
  list(
    outcome    = data[[vars[1]]],
    treatment  = data[[vars[2]]],
    instrument = data[[vars[3]]]
  )
}

#' Validate and coerce IV data inputs
#'
#' Checks that the treatment and instrument are binary (0/1) and that the
#' outcome is numeric. If treatment or instrument have the correct values but
#' wrong storage type they are silently coerced to numeric.
#'
#' @param data    A data frame.
#' @param formula A formula of the form \code{outcome ~ treatment | instrument}.
#' @return A named list with elements \code{outcome}, \code{treatment}, and
#'   \code{instrument}, all of class \code{numeric}.
#' @keywords internal
validate_iv_data <- function(data, formula) {
  d <- get_data_from_formula(data, formula)

  for (var in c("treatment", "instrument")) {
    x     <- d[[var]]
    x_num <- suppressWarnings(as.numeric(x))

    if (anyNA(x_num)) {
      stop(sprintf(
        "'%s' cannot be coerced to numeric. Please check your data.", var
      ))
    }
    if (!all(x_num %in% c(0, 1))) {
      stop(sprintf(
        "'%s' must be binary (0/1 only). Found values: %s.",
        var, paste(sort(unique(x_num)), collapse = ", ")
      ))
    }
    d[[var]] <- x_num
  }

  outcome_num <- suppressWarnings(as.numeric(d$outcome))
  if (anyNA(outcome_num)) {
    stop("'outcome' cannot be coerced to numeric. Please check your data.")
  }
  d$outcome <- outcome_num

  d
}

#' Compute identified type shares and potential outcome means
#'
#' Under the standard IV assumptions of Angrist, Imbens, and Rubin (1996), the
#' compliance-type shares and the conditional potential outcome means that are
#' identifiable from observed data are recovered here (see Appendix A of Locher
#' and Stensrud (2025) for derivations).
#'
#' @param y Numeric vector of outcomes.
#' @param a Numeric binary vector of treatments (0/1).
#' @param z Numeric binary vector of instruments (0/1).
#' @return A named list:
#'   \describe{
#'     \item{p_at}{Share of always-takers \eqn{\pi_A}.}
#'     \item{p_nt}{Share of never-takers \eqn{\pi_N}.}
#'     \item{p_c}{Share of compliers \eqn{\pi_C}.}
#'     \item{y1_c}{Complier mean under treatment \eqn{\mu_C^{(1)}}.}
#'     \item{y0_c}{Complier mean under control \eqn{\mu_C^{(0)}}.}
#'     \item{y1_at}{Always-taker mean under treatment \eqn{\mu_A^{(1)}}.}
#'     \item{y0_nt}{Never-taker mean under control \eqn{\mu_N^{(0)}}.}
#'     \item{at_exist}{Logical; \code{TRUE} if always-takers are present.}
#'     \item{nt_exist}{Logical; \code{TRUE} if never-takers are present.}
#'   }
#' @keywords internal
get_known_quantities <- function(y, a, z) {
  p_at <- mean(a[z == 0])
  p_nt <- 1 - mean(a[z == 1])
  p_c  <- mean(a[z == 1]) - mean(a[z == 0])

  at_exist <- p_at > 0
  nt_exist <- p_nt > 0

  y1_at <- if (at_exist) mean(y[a == 1 & z == 0]) else 0
  y0_nt <- if (nt_exist) mean(y[a == 0 & z == 1]) else 0

  y1_c <- ((p_c + p_at) * mean(y[a == 1 & z == 1]) - p_at * y1_at) / p_c
  y0_c <- ((p_c + p_nt) * mean(y[a == 0 & z == 0]) - p_nt * y0_nt) / p_c

  list(
    p_at     = p_at,  p_nt  = p_nt,  p_c   = p_c,
    y0_c     = y0_c,  y1_c  = y1_c,
    y1_at    = y1_at, y0_nt = y0_nt,
    at_exist = at_exist, nt_exist = nt_exist
  )
}
