#' Expand a table of counts into individual-level data
#'
#' With a binary instrument, treatment, and outcome, the data are fully
#' summarised by the counts in the \eqn{2 \times 2 \times 2} table, as often
#' reported in trial publications. `expand_counts()` turns such a table into
#' one row per individual, which can then be passed to [ate_sensitivity()],
#' [ate_bounds()], or [ate_bayes()].
#'
#' @param counts A data frame with one row per cell and a column of counts,
#'   or a contingency table (e.g. from [table()] or [xtabs()]).
#' @param count Name of the column holding the counts. Defaults to `"n"` if
#'   present, and to `"Freq"` otherwise (the name used by
#'   [as.data.frame.table()]).
#'
#' @return A data frame with the non-count columns of `counts`, each row
#'   repeated as often as its count.
#'
#' @examples
#' # Lipid data (Efron and Feldman, 1991), as analysed by Richardson, Evans,
#' # and Robins (2011): one-sided non-compliance
#' lipid <- data.frame(
#'   z = c(0, 0, 0, 0, 1, 1, 1, 1),
#'   a = c(0, 0, 1, 1, 0, 0, 1, 1),
#'   y = c(0, 1, 0, 1, 0, 1, 0, 1),
#'   n = c(158, 14, 0, 0, 52, 12, 23, 78)
#' )
#' dat <- expand_counts(lipid)
#' nrow(dat)
#'
#' @export
expand_counts <- function(counts, count = NULL) {
  if (inherits(counts, "table")) counts <- as.data.frame(counts, stringsAsFactors = FALSE)
  if (!is.data.frame(counts)) stop("`counts` must be a data frame or a table.", call. = FALSE)
  if (is.null(count)) {
    count <- if ("n" %in% names(counts)) "n" else if ("Freq" %in% names(counts)) "Freq" else
      stop("Specify the column holding the counts via `count`.", call. = FALSE)
  }
  if (!count %in% names(counts)) stop(sprintf("Column '%s' not found.", count), call. = FALSE)
  n <- counts[[count]]
  if (!is.numeric(n) || anyNA(n) || any(n < 0) || any(n != round(n))) {
    stop("Counts must be non-negative whole numbers.", call. = FALSE)
  }
  out <- counts[rep(seq_len(nrow(counts)), n), setdiff(names(counts), count), drop = FALSE]
  for (v in names(out)) {
    if (is.factor(out[[v]])) out[[v]] <- as.character(out[[v]])
    num <- suppressWarnings(as.numeric(out[[v]]))
    if (!anyNA(num[!is.na(out[[v]])])) out[[v]] <- num
  }
  rownames(out) <- NULL
  out
}
