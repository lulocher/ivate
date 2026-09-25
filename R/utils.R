# Internal helpers -- not exported.

#' Parse a formula of the form `outcome ~ treatment | instrument`
#'
#' @param formula A formula or a character string coercible to one.
#' @return A named character vector with elements `outcome`, `treatment`, and
#'   `instrument`.
#' @keywords internal
#' @noRd
parse_iv_formula <- function(formula) {
  if (is.character(formula)) formula <- stats::as.formula(formula)
  msg <- "`formula` must be of the form `outcome ~ treatment | instrument`."
  if (!inherits(formula, "formula") || length(formula) != 3L) stop(msg, call. = FALSE)

  rhs <- formula[[3L]]
  if (!is.call(rhs) || !identical(rhs[[1L]], as.name("|"))) stop(msg, call. = FALSE)

  vars <- list(outcome = formula[[2L]], treatment = rhs[[2L]], instrument = rhs[[3L]])
  if (!all(vapply(vars, is.name, logical(1)))) {
    stop("Each part of `formula` must be a single variable name.", call. = FALSE)
  }
  vapply(vars, as.character, character(1))
}

#' Validate and extract IV data
#'
#' Extracts outcome, treatment, and instrument from `data`, drops incomplete
#' observations (with a warning), coerces to numeric, and checks that
#' treatment and instrument are binary (0/1) and that both instrument arms are
#' observed.
#'
#' @param data A data frame.
#' @param formula A formula of the form `outcome ~ treatment | instrument`.
#' @return A list with numeric vectors `y`, `a`, `z` and the variable names
#'   `vars`.
#' @keywords internal
#' @noRd
prepare_iv_data <- function(data, formula) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  vars <- parse_iv_formula(formula)

  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0) {
    stop(sprintf("Variable(s) not found in `data`: %s.",
                 paste(missing_vars, collapse = ", ")), call. = FALSE)
  }

  d <- lapply(vars, function(v) {
    x <- data[[v]]
    if (is.factor(x)) x <- as.character(x)
    x_num <- suppressWarnings(as.numeric(x))
    if (any(is.na(x_num) & !is.na(x))) {
      stop(sprintf("'%s' cannot be coerced to numeric. Please check your data.", v),
           call. = FALSE)
    }
    x_num
  })

  complete <- stats::complete.cases(d$outcome, d$treatment, d$instrument)
  if (!all(complete)) {
    warning(sprintf("Removed %d observation(s) with missing values.", sum(!complete)),
            call. = FALSE)
    d <- lapply(d, `[`, complete)
  }

  for (part in c("treatment", "instrument")) {
    x <- d[[part]]
    if (!all(x %in% c(0, 1))) {
      vals <- sort(unique(x))
      stop(sprintf("The %s '%s' must be binary (0/1 only). Found values: %s.",
                   part, vars[[part]], paste(vals[seq_len(min(10, length(vals)))], collapse = ", ")),
           call. = FALSE)
    }
  }
  if (!all(c(0, 1) %in% d$instrument)) {
    stop("Both instrument arms (0 and 1) must be observed.", call. = FALSE)
  }

  list(y = d$outcome, a = d$treatment, z = d$instrument, vars = vars)
}

#' Check whether a vector is binary (0/1)
#' @keywords internal
#' @noRd
is_binary <- function(y) all(y %in% c(0, 1))

#' Resolve the outcome range \eqn{[y_{min}, y_{max}]}
#'
#' Defaults to `c(0, 1)` for binary outcomes and to the observed range
#' otherwise.
#' @keywords internal
#' @noRd
resolve_y_range <- function(y, y_range = NULL) {
  if (is.null(y_range)) {
    return(if (is_binary(y)) c(0, 1) else range(y))
  }
  if (!is.numeric(y_range) || length(y_range) != 2L || anyNA(y_range) ||
      y_range[1] >= y_range[2]) {
    stop("`y_range` must be a numeric vector c(y_min, y_max) with y_min < y_max.",
         call. = FALSE)
  }
  if (min(y) < y_range[1] || max(y) > y_range[2]) {
    stop("Observed outcomes lie outside `y_range`.", call. = FALSE)
  }
  y_range
}

#' Format numbers with three significant digits
#' @keywords internal
#' @noRd
fmt_num <- function(x, digits = 3) {
  out <- formatC(x, digits = digits, format = "fg")
  out[is.na(x)] <- "NA"
  trimws(out)
}

#' Format a confidence interval as "(lower, upper)" in square brackets
#' @keywords internal
#' @noRd
fmt_ci <- function(lower, upper) {
  ifelse(is.na(lower) & is.na(upper), "",
         paste0("[", fmt_num(lower), ", ", fmt_num(upper), "]"))
}
