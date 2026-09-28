# Data with P(A = 1 | Z = 1) = P(A = 1 | Z = 0), i.e. no compliers (pi_C = 0)
no_complier_data <- function() {
  expand_counts(data.frame(
    z = c(0, 0, 0, 0, 1, 1, 1, 1), a = c(0, 0, 1, 1, 0, 0, 1, 1),
    y = c(0, 1, 0, 1, 0, 1, 0, 1), n = c(30, 20, 30, 20, 30, 20, 30, 20)
  ))
}

test_that("contributions pi_C * mu_C are identified and equal pi_C * mu_C when pi_C > 0", {
  d <- sim_iv()
  q <- identify_quantities(d$y, d$a, d$z)
  expect_equal(q[["pimu_C1"]], q[["pi_C"]] * q[["mu_C1"]])
  expect_equal(q[["pimu_C0"]], q[["pi_C"]] * q[["mu_C0"]])
  # pi_C * LATE equals the intention-to-treat effect on the outcome
  expect_equal(complier_term(q), mean(d$y[d$z == 1]) - mean(d$y[d$z == 0]))
})

test_that("thresholds and bounds are computed without compliers", {
  d <- no_complier_data()
  q <- identify_quantities(d$y, d$a, d$z)
  expect_equal(q[["pi_C"]], 0)
  expect_true(is.na(q[["mu_C0"]]) && is.na(q[["mu_C1"]]))
  expect_equal(unname(q[c("pimu_C0", "pimu_C1")]), c(0, 0))

  expect_warning(b <- ate_bounds(d, y ~ a | z, ci = FALSE), "No compliers")
  # Natural bounds: P(Y=1, A=1 | Z=1) - P(Y=1, A=0 | Z=0) - P(A=1 | Z=0), and
  #                 P(Y=1, A=1 | Z=1) + P(A=0 | Z=1) - P(Y=1, A=0 | Z=0)
  expect_equal(b$bounds$estimate, c(0.2 - 0.2 - 0.5, 0.2 + 0.5 - 0.2))
  expect_false(any(b$model_check$sample))
  expect_output(print(b), "LATE:     undefined \\(no compliers\\)")

  expect_warning(s <- ate_sensitivity(d, y ~ a | z, theta_A0 = 1, ci = FALSE), "No compliers")
  expect_equal(ate_at(s$q, 1, s$threshold$estimate), 0)
  expect_false(any(s$identified$quantity %in% c("mu_C0", "mu_C1", "late")))
  expect_true(all(is.na(s$implied$benchmark_value[s$implied$comparison != "ratio"][1])))
  expect_no_error(ggplot2::ggplot_build(suppressWarnings(plot_comparisons(s))))

  # Restrictions that refer to complier quantities are undefined
  expect_error(suppressWarnings(ate_bounds(d, y ~ a | z, theta_N1 = restrict_effect(upper = ~ late),
                                           ci = FALSE)), "undefined")
})

test_that("a negative share of compliers is reported, not fatal", {
  d <- sim_iv(1000)
  expect_warning(b <- ate_bounds(transform(d, z = 1 - z), y ~ a | z, ci = FALSE), "recode")
  expect_true(b$model_check$sample[["pi_C"]])
})
