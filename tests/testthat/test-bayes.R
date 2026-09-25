lipid <- expand_counts(lipid_counts)
flu   <- expand_counts(flu_counts)

test_that("constraint-violation probabilities reproduce Richardson, Evans, and Robins (2011)", {
  set.seed(11)
  # Flu vaccine data, uniform prior on the saturated model: Table 7, first row
  v_flu <- sample_dirichlet_q(arm_counts(flu$y, flu$a, flu$z), "two-sided", 20000, 1)$violation
  expect_equal(v_flu[["(11)"]], 0)
  expect_lt(abs(v_flu[["(12)"]] - 0.0603), 0.015)
  expect_lt(abs(v_flu[["(13)"]] - 0.5411), 0.015)

  # Lipid data: posterior probability 0.38, prior probability 0.5 (Section 3.1)
  v_lip <- sample_dirichlet_q(arm_counts(lipid$y, lipid$a, lipid$z), "never-takers only", 20000, 1)$violation
  expect_lt(abs(v_lip[["any"]] - 0.38), 0.015)
  v_prior <- sample_dirichlet_q(NULL, "never-takers only", 20000, 1)$violation
  expect_lt(abs(v_prior[["any"]] - 0.5), 0.015)
})

test_that("posterior draws satisfy the IV constraints and respect structural zeros", {
  set.seed(12)
  qd <- sample_dirichlet_q(arm_counts(flu$y, flu$a, flu$z), "two-sided", 2000, 1)$q
  expect_true(all(qd[, "pi_C"] > 0))
  expect_true(all(qd[, c("mu_C0", "mu_C1", "mu_A1", "mu_N0")] >= 0 &
                  qd[, c("mu_C0", "mu_C1", "mu_A1", "mu_N0")] <= 1))

  ql <- sample_dirichlet_q(arm_counts(lipid$y, lipid$a, lipid$z), "never-takers only", 2000, 1)$q
  expect_true(all(ql[, "pi_A"] == 0))
  expect_true(all(is.na(ql[, "mu_A1"])))
})

test_that("identified quantities from cell probabilities match the data-based formulas", {
  cnt <- arm_counts(flu$y, flu$a, flu$z)
  q_probs <- q_from_probs(matrix(cnt$z0 / sum(cnt$z0), 1), matrix(cnt$z1 / sum(cnt$z1), 1))
  expect_equal(q_probs[1, ], identify_quantities(flu$y, flu$a, flu$z))
})

test_that("Bayesian inference is restricted to binary outcomes", {
  d <- sim_iv(continuous = TRUE)
  expect_error(ate_bayes(d, y ~ a | z, theta_A0 = 0.5, theta_N1 = 0.5), "binary")
  expect_error(ate_sensitivity(d, y ~ a | z, inference = "bayes"), "binary")
  expect_error(ate_bounds(d, y ~ a | z, inference = "bayes"), "binary")
  expect_error(ate_bounds(lipid, y ~ a | z, inference = "bayes", y_range = c(0, 2)), "fixed at c\\(0, 1\\)")
})

test_that("ate_sensitivity and ate_bounds support Bayesian inference", {
  set.seed(13)
  res <- ate_sensitivity(lipid, y ~ a | z, inference = "bayes", n_draws = 500)
  expect_equal(res$inference, "bayes")
  expect_equal(nrow(res$draws), 500)
  expect_named(res$bayes$violation, c("(11)", "(12)", "(13)", "any"))
  # Estimates are posterior medians of the draws
  tau_draws <- apply(res$draws, 1, function(q) threshold_N1(q, 0))
  expect_equal(res$threshold$estimate, stats::median(tau_draws))
  expect_output(print(res), "CrI")

  set.seed(14)
  b <- suppressWarnings(ate_bounds(flu, y ~ a | z, theta_A0 = restrict_effect(lower = 0),
                                   inference = "bayes", n_draws = 500))
  expect_true(all(b$bounds$ci_lower <= b$bounds$estimate & b$bounds$estimate <= b$bounds$ci_upper))
  expect_output(print(b), "Model check")
})

test_that("Bayesian and bootstrap intervals agree in large samples", {
  d <- sim_iv(n = 20000, p_at = 0, seed = 3)
  set.seed(15)
  boot  <- ate_sensitivity(d, y ~ a | z, n_boot = 300)$threshold
  bayes <- ate_sensitivity(d, y ~ a | z, inference = "bayes", n_draws = 1000)$threshold
  expect_lt(abs(bayes$estimate - boot$estimate), 0.02)
  width <- function(t) t$ci_upper - t$ci_lower
  expect_lt(abs(width(bayes) / width(boot) - 1), 0.25)
})

test_that("ate_bayes sensitivity mode gives the posterior of the ATE at supplied values", {
  set.seed(16)
  res <- ate_bayes(lipid, y ~ a | z, theta_N1 = c(0, 0.5, 1), n_draws = 1000)
  expect_s3_class(res, "ivate_bayes")
  expect_equal(res$mode, "sensitivity")
  expect_equal(res$grid$theta_N1, c(0, 0.5, 1))
  expect_true(all(is.na(res$grid$theta_A0)))

  lin <- ate_linear(res$q_draws)
  ate_05 <- lin$intercept + lin$slope_N * 0.5
  expect_equal(res$grid$median[2], stats::quantile(ate_05, 0.5, names = FALSE))
  expect_equal(res$grid$prob[2], mean(ate_05 > 0))
  # Simultaneous bands contain pointwise intervals
  expect_true(all(res$grid$band_lower <= res$grid$ci_lower + 1e-12))
  expect_true(all(res$grid$band_upper >= res$grid$ci_upper - 1e-12))
  # At the extremes of theta_N1 the ATE equals the identified bounds
  expect_equal(res$grid$median[c(1, 3)], res$bounds$estimate, tolerance = 1e-8)

  two <- suppressWarnings(ate_bayes(flu, y ~ a | z, theta_A0 = c(0, 1), theta_N1 = c(0, 0.5, 1),
                                    n_draws = 500))
  expect_equal(nrow(two$grid), 6)
  one <- suppressWarnings(ate_bayes(flu, y ~ a | z, theta_A0 = 1, theta_N1 = 0.5, n_draws = 500))
  expect_equal(nrow(one$draws), 500)
})

test_that("ate_bayes full mode combines posterior draws with priors", {
  set.seed(17)
  # Point-mass priors reproduce the ATE at fixed values draw by draw
  res <- suppressWarnings(ate_bayes(flu, y ~ a | z, theta_A0 = prior_beta(), theta_N1 = 0.4,
                                    n_draws = 500))
  expect_equal(res$mode, "full")
  lin <- ate_linear(res$q_draws)
  expect_equal(res$draws$ate, lin$intercept + lin$slope_A * res$draws$theta_A0 + lin$slope_N * 0.4)
  expect_true(all(res$draws$theta_N1 == 0.4))
  expect_true(all(c("median", "ci_lower", "ci_upper", "prob") %in% names(res$posterior)))
  expect_false(is.null(res$prior_summary))

  # The posterior lies within the posterior range of the identified bounds
  res2 <- ate_bayes(lipid, y ~ a | z, theta_N1 = prior_uniform(), n_draws = 500, prior_draws = FALSE)
  expect_null(res2$prior_summary)
  expect_true(res2$posterior$ci_lower >= min(res2$bounds$ci_lower))
  expect_true(res2$posterior$ci_upper <= max(res2$bounds$ci_upper))
})

test_that("priors are sampled on the right scale and within [0, 1]", {
  set.seed(18)
  qd <- sample_dirichlet_q(arm_counts(flu$y, flu$a, flu$z), "two-sided", 300, 1)$q
  b <- sample_prior(prior_beta(2, 5), "theta_N1", qd)
  expect_true(all(b >= 0 & b <= 1))

  u <- sample_prior(prior_uniform(restrict_effect(lower = 0)), "theta_N1", qd)
  expect_true(all(u >= qd[, "mu_N0"] - 1e-12))
  uA <- sample_prior(prior_uniform(restrict_effect(lower = 0)), "theta_A0", qd)
  expect_true(all(uA <= qd[, "mu_A1"] + 1e-12))

  late <- qd[, "mu_C1"] - qd[, "mu_C0"]
  n_eff <- sample_prior(prior_normal(~ late, 1e-6, on = "effect"), "theta_N1", qd)
  target <- pmin(pmax(qd[, "mu_N0"] + late, 0), 1)
  expect_equal(n_eff, target, tolerance = 1e-4)
  a_eff <- sample_prior(prior_normal(0.05, 1e-6, on = "effect"), "theta_A0", qd)
  expect_equal(a_eff, qd[, "mu_A1"] - 0.05, tolerance = 1e-4)

  r <- sample_prior(prior_normal(1, 0.5, on = "ratio"), "theta_N1", qd)
  expect_true(all(r >= 0 & r <= 1))
  h <- sample_prior(prior_normal(0, 0.1, on = "heterogeneity"), "theta_A0", qd)
  expect_true(all(h >= 0 & h <= 1))
})

test_that("ate_bayes validates its sensitivity-parameter arguments", {
  expect_error(suppressWarnings(ate_bayes(flu, y ~ a | z, theta_N1 = 0.5)), "Supply `theta_A0`")
  expect_error(suppressWarnings(ate_bayes(flu, y ~ a | z, theta_A0 = c(0, 1), theta_N1 = prior_beta())),
               "single value or a prior")
  expect_error(ate_bayes(lipid, y ~ a | z, theta_N1 = restrict_effect(lower = 0)), "prior_uniform")
  expect_error(ate_bayes(lipid, y ~ a | z, theta_N1 = 1.5), "in \\[0, 1\\]")
  expect_warning(ate_bayes(lipid, y ~ a | z, theta_A0 = 0.5, theta_N1 = 0.5, n_draws = 100), "ignored")
  expect_error(prior_beta(0, 1), "positive")
  expect_error(prior_normal(0, -1), "positive")
})

test_that("print and plot methods work for ate_bayes results", {
  set.seed(19)
  sens <- ate_bayes(lipid, y ~ a | z, theta_N1 = seq(0, 1, 0.1), n_draws = 300)
  expect_output(print(sens), "Posterior of the ATE given")
  for (p in list(plot(sens), plot(sens, type = "probability"))) {
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
  surf <- suppressWarnings(ate_bayes(flu, y ~ a | z, theta_A0 = seq(0, 1, 0.25),
                                     theta_N1 = seq(0, 1, 0.25), n_draws = 300))
  expect_no_error(suppressWarnings(ggplot2::ggplot_build(plot(surf))))
  expect_no_error(suppressWarnings(ggplot2::ggplot_build(plot(surf, type = "ate"))))

  full <- ate_bayes(lipid, y ~ a | z, theta_N1 = prior_beta(), n_draws = 300)
  expect_output(print(full), "Priors on the sensitivity parameters")
  expect_no_error(ggplot2::ggplot_build(plot(full)))
  expect_no_error(ggplot2::ggplot_build(plot(full, show_prior = FALSE)))

  sens_b <- ate_sensitivity(lipid, y ~ a | z, inference = "bayes", n_draws = 300)
  expect_no_error(ggplot2::ggplot_build(plot(sens_b)))
})

test_that("expand_counts turns tables into individual-level data", {
  d <- expand_counts(lipid_counts)
  expect_equal(nrow(d), sum(lipid_counts$n))
  expect_named(d, c("z", "a", "y"))
  expect_equal(sum(d$z == 1 & d$a == 1 & d$y == 1), 78)

  tab <- xtabs(n ~ z + a + y, data = lipid_counts)
  d_tab <- expand_counts(tab)
  expect_equal(nrow(d_tab), nrow(d))
  expect_equal(ate_bounds(d_tab, y ~ a | z, ci = FALSE)$bounds, ate_bounds(d, y ~ a | z, ci = FALSE)$bounds)

  expect_error(expand_counts(data.frame(z = 1, k = 1)), "Specify the column")
  expect_error(expand_counts(data.frame(z = 1, n = -1)), "non-negative")
})
