status_of <- function(res) {
  n <- comparison_nodes(res, ci = FALSE)
  stats::setNames(n$status, n$id)
}

test_that("node statuses reflect estimated, assumed, and implied quantities", {
  d <- sim_iv()
  unrestricted <- ate_sensitivity(d, y ~ a | z, ci = FALSE)
  s <- status_of(unrestricted)
  expect_equal(unname(s[c("C1", "C0", "A1", "N0")]), rep("estimated", 4))
  expect_equal(unname(s[c("A0", "N1")]), c("implied", "implied"))

  step1 <- ate_sensitivity(d, y ~ a | z, theta_A0 = 1, ci = FALSE)
  s1 <- status_of(step1)
  expect_equal(unname(s1[c("A0", "N1")]), c("assumed", "implied"))

  step2 <- ate_sensitivity(d, y ~ a | z, theta_N1 = restrict_effect(lower = 0), ci = FALSE)
  s2 <- status_of(step2)
  expect_equal(unname(s2[c("A0", "N1")]), c("implied", "assumed"))
})

test_that("values are filled in from the analysis", {
  d <- sim_iv()
  res <- ate_sensitivity(d, y ~ a | z, theta_A0 = 1, ci = FALSE)
  n <- comparison_nodes(res, ci = FALSE)
  expect_equal(n$value[n$id == "N1"], paste(">", fmt_digits(res$threshold$estimate, 2)))
  expect_equal(n$value[n$id == "A0"], "1.0")
  expect_equal(n$value[n$id == "C1"], fmt_digits(res$q[["mu_C1"]], 2))

  # Unrestricted two-sided: boundary formula, drawn as plotmath
  n0 <- comparison_nodes(ate_sensitivity(d, y ~ a | z, ci = FALSE), ci = FALSE)
  expect_true(all(n0$value_parse[n0$id %in% c("A0", "N1")]))
  expect_no_error(parse(text = n0$value[n0$id == "N1"]))
})

test_that("absent compliance types are dropped and arrows follow Figure 1", {
  one_sided <- ate_sensitivity(sim_iv(p_at = 0), y ~ a | z, ci = FALSE)
  n <- comparison_nodes(one_sided, ci = FALSE)
  expect_setequal(n$id, c("C1", "C0", "N1", "N0"))
  expect_equal(nrow(comparison_arrows(n)), 4)

  two_sided <- ate_sensitivity(sim_iv(), y ~ a | z, theta_A0 = 1, ci = FALSE)
  a <- comparison_arrows(comparison_nodes(two_sided, ci = FALSE))
  expect_equal(nrow(a), 9)
  # Arrows touching the implied node are red-status and dashed
  expect_true(all(a$linetype[a$status == "implied"] == "dashed"))
  expect_equal(sum(a$status == "estimated" & a$linetype == "solid"), 3)
})

test_that("plot_comparisons returns a ggplot, also via plot(type = 'comparisons')", {
  set.seed(1)
  res <- ate_sensitivity(sim_iv(), y ~ a | z, theta_A0 = restrict_effect(lower = 0), n_boot = 50)
  for (p in list(plot_comparisons(res), plot_comparisons(res, ci = TRUE),
                 plot(res, type = "comparisons"),
                 plot_comparisons(ate_sensitivity(sim_iv(), y ~ a | z, ci = FALSE)),
                 plot_comparisons(ate_sensitivity(sim_iv(p_nt = 0), y ~ a | z, ci = FALSE)))) {
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
  set.seed(2)
  res_b <- ate_sensitivity(expand_counts(lipid_counts), y ~ a | z, inference = "bayes", n_draws = 200)
  expect_no_error(ggplot2::ggplot_build(plot_comparisons(res_b, ci = TRUE)))

  expect_error(plot_comparisons(ate_bounds(sim_iv(), y ~ a | z, ci = FALSE)), "ate_sensitivity")
  expect_error(plot_comparisons(res, colours = c(a = "red")), "named")
})

test_that("fmt_digits shows significant digits", {
  expect_equal(fmt_digits(c(0.703, 0.011, 0.000836, 0.00505, 0.00488), 2),
               c("0.70", "0.011", "0.00084", "0.0050", "0.0049"))
  expect_equal(fmt_digits(c(-1.24, 4.46, 1, 12, 0), 2), c("-1.2", "4.5", "1.0", "12", "0"))
  expect_equal(fmt_digits(c(0.4567, 0.011), 3), c("0.457", "0.0110"))
  n <- comparison_nodes(ate_sensitivity(sim_iv(), y ~ a | z, theta_A0 = 1, ci = FALSE), ci = FALSE, digits = 3)
  expect_equal(n$value[n$id == "A0"], "1.00")
  expect_error(plot_comparisons(ate_sensitivity(sim_iv(), y ~ a | z, ci = FALSE), digits = -1), "digits")
})

