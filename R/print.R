# Print methods.

has_ci <- function(tab) !is.null(tab) && "ci_lower" %in% names(tab)

setting_label <- function(setting) {
  switch(setting,
    "two-sided"          = "two-sided non-compliance",
    "never-takers only"  = "one-sided non-compliance (never-takers only)",
    "always-takers only" = "one-sided non-compliance (always-takers only)",
    setting
  )
}

ci_header <- function(level, inference = "bootstrap") {
  paste0(format(100 * level), if (identical(inference, "bayes")) "% CrI" else "% CI")
}

print_bayes_info <- function(x) {
  b <- x$bayes
  if (is.null(b)) return(invisible())
  cat(sprintf(paste0(
    "  Inference: Bayesian; Dirichlet(%s) prior on p(y, a | z) restricted to the IV model,\n",
    "             %d posterior draws; estimates are posterior medians\n"),
    fmt_num(b$prior_obs), b$n_draws))
  v <- b$violation
  cat(sprintf(paste0(
    "  Model check: P(IV inequalities violated | data) = %s  [(11): %s, (12): %s, (13): %s]\n"),
    fmt_num(v[["any"]]), fmt_num(v[["(11)"]]), fmt_num(v[["(12)"]]), fmt_num(v[["(13)"]])))
  invisible()
}

print_model_check <- function(x) {
  m <- x$model_check
  if (is.null(m) || identical(x$inference, "bayes")) return(invisible())
  rng <- paste0("[", fmt_num(m$y_range[1]), ", ", fmt_num(m$y_range[2]), "]")
  labels <- c(pi_C = "pi_C > 0", mu_C0 = paste("mu_C0 in", rng), mu_C1 = paste("mu_C1 in", rng))
  sample_txt <- if (any(m$sample)) {
    paste("VIOLATED in the sample:", paste(labels[m$sample], collapse = "; "))
  } else {
    "satisfied in the sample"
  }
  boot_txt <- ""
  if (!is.null(m$bootstrap)) {
    pct <- function(p) paste0(fmt_num(100 * p), "%")
    b <- m$bootstrap
    boot_txt <- sprintf("; violated in %s of %d bootstrap samples [pi_C: %s, mu_C0: %s, mu_C1: %s]",
                        pct(b[["any"]]), m$n_boot, pct(b[["pi_C"]]), pct(b[["mu_C0"]]), pct(b[["mu_C1"]]))
  }
  cat(sprintf("  Model check: IV restrictions %s%s\n", sample_txt, boot_txt))
  invisible()
}

print_header <- function(x, title) {
  cat(title, "\n\n", sep = "")
  cat(sprintf("  Model:    %s ~ %s | %s  (n = %d)\n",
              x$vars[["outcome"]], x$vars[["treatment"]], x$vars[["instrument"]], x$n))
  shares <- x$q[c("pi_C", "pi_A", "pi_N")]
  cat(sprintf("  Setting:  %s (pi_C = %s, pi_A = %s, pi_N = %s)\n",
              setting_label(x$setting), fmt_num(shares[1]), fmt_num(shares[2]), fmt_num(shares[3])))
  late <- x$identified[x$identified$quantity == "late", ]
  cat(sprintf("  LATE:     %s %s\n", fmt_num(late$estimate),
              if (has_ci(late)) paste(ci_header(x$level, x$inference), fmt_ci(late$ci_lower, late$ci_upper)) else ""))
}

#' @export
print.ivate_sensitivity <- function(x, ...) {
  print_header(x, "Sensitivity analysis for the ATE with an instrumental variable")
  rel <- if (x$direction == "greater") ">" else "<"
  cat(sprintf("  Criterion: ATE %s %s\n", rel, fmt_num(x$cutoff)))
  print_bayes_info(x)
  print_model_check(x)

  if (!is.null(x$boundary)) {
    cat("\nBoundary B (ATE = ", fmt_num(x$cutoff), "); the criterion holds iff:\n", sep = "")
    for (i in seq_len(nrow(x$boundary))) {
      b <- x$boundary[i, ]
      sign <- if (b$slope < 0) "-" else "+"
      cat(sprintf("  %-17s %s %s %s %s * (%s)\n", b$lhs, b$relation, fmt_num(b$intercept),
                  sign, fmt_num(abs(b$slope)), b$rhs))
    }
  }

  if (!is.null(x$fixed)) {
    f <- x$fixed
    cat(sprintf("\nRestriction: %s\n  least favourable value %s = %s %s\n",
                f$restriction, f$parameter, fmt_num(f$value),
                if (has_ci(f) && !isTRUE(f$ci_lower == f$ci_upper)) {
                  paste(ci_header(x$level, x$inference), fmt_ci(f$ci_lower, f$ci_upper))
                } else ""))
  }

  if (!is.null(x$threshold)) {
    t <- x$threshold
    cat(sprintf("\nThreshold: the criterion holds iff %s %s %s %s\n",
                t$parameter, t$relation, fmt_num(t$estimate),
                if (has_ci(t)) paste(ci_header(x$level, x$inference), fmt_ci(t$ci_lower, t$ci_upper)) else ""))
    if (t$estimate < x$y_range[1] || t$estimate > x$y_range[2]) {
      holds <- (t$relation == ">") == (t$estimate < x$y_range[1])
      cat(sprintf("  The threshold lies outside the outcome range [%s, %s]: the criterion %s for all feasible values.\n",
                  fmt_num(x$y_range[1]), fmt_num(x$y_range[2]), if (holds) "holds" else "fails"))
    }

    cat("\nImplied comparisons at the threshold (implied vs. identified benchmark):\n")
    im <- x$implied
    implied <- paste(im$implied, "=", fmt_num(im$implied_value))
    if (has_ci_cols(im, "implied")) {
      implied <- paste(implied, fmt_ci(im$implied_ci_lower, im$implied_ci_upper))
    }
    benchmark <- paste(im$benchmark, "=", fmt_num(im$benchmark_value))
    relative  <- ifelse(is.na(im$relative), "", paste0("(ratio ", fmt_num(im$relative), ")"))
    cat(paste0("  ", format(im$comparison), "  ", format(implied), "  vs  ",
               format(benchmark), "  ", relative, "\n"), sep = "")
  } else if (!is.null(x$boundary)) {
    cat("\nRestrict `theta_A0` or `theta_N1` to obtain a scalar threshold for the other parameter.\n")
  }
  invisible(x)
}

has_ci_cols <- function(tab, prefix) paste0(prefix, "_ci_lower") %in% names(tab)

#' @export
print.ivate_bounds <- function(x, ...) {
  print_header(x, "Bounds on the ATE with an instrumental variable")
  print_bayes_info(x)
  print_model_check(x)

  cat("\nFeasible sets:\n")
  for (i in seq_len(nrow(x$feasible))) {
    f <- x$feasible[i, ]
    cat(sprintf("  %s in [%s, %s]   restriction: %s\n", f$parameter,
                fmt_num(f$lower), fmt_num(f$upper), f$restriction))
  }

  show <- function(tab, label) {
    lo <- tab[tab$bound == "lower", ]
    up <- tab[tab$bound == "upper", ]
    cat(sprintf("  %-13s [%s, %s]\n", label, fmt_num(lo$estimate), fmt_num(up$estimate)))
    if (has_ci(tab)) {
      cat(sprintf("  %-13s lower %s, upper %s\n", ci_header(x$level, x$inference),
                  fmt_ci(lo$ci_lower, lo$ci_upper), fmt_ci(up$ci_lower, up$ci_upper)))
    }
  }
  cat("\nBounds on the ATE:\n")
  show(x$bounds, "Restricted:")
  show(x$unrestricted, "Unrestricted:")
  invisible(x)
}
