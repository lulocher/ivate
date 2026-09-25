test_that("one-sided threshold makes the ATE equal the cutoff", {
  d <- sim_iv(p_at = 0)
  q <- manual_q(d)
  for (cutoff in c(0, 0.05)) {
    res <- ate_sensitivity(d, y ~ a | z, cutoff = cutoff, ci = FALSE)
    expect_s3_class(res, "ivate_sensitivity")
    expect_equal(res$setting, "never-takers only")
    expect_equal(res$threshold$parameter, "theta_N1")
    expect_equal(ate_manual(q, NA, res$threshold$estimate), cutoff)
    expect_null(res$boundary)
    expect_null(res$fixed)
  }

  d_at <- sim_iv(p_nt = 0)
  q_at <- manual_q(d_at)
  res_at <- ate_sensitivity(d_at, y ~ a | z, cutoff = 0.1, ci = FALSE)
  expect_equal(res_at$threshold$parameter, "theta_A0")
  expect_equal(ate_manual(q_at, res_at$threshold$estimate, NA), 0.1)
})

test_that("relations follow the monotonicity of the ATE in each parameter", {
  d <- sim_iv(p_at = 0)
  expect_equal(ate_sensitivity(d, y ~ a | z, ci = FALSE)$threshold$relation, ">")
  expect_equal(ate_sensitivity(d, y ~ a | z, direction = "less", ci = FALSE)$threshold$relation, "<")
  d_at <- sim_iv(p_nt = 0)
  expect_equal(ate_sensitivity(d_at, y ~ a | z, ci = FALSE)$threshold$relation, "<")
})

test_that("implied comparisons are computed at the threshold", {
  d <- sim_iv(p_at = 0)
  q <- manual_q(d)
  res <- ate_sensitivity(d, y ~ a | z, ci = FALSE)
  tau <- res$threshold$estimate
  im  <- res$implied
  expect_equal(im$implied_value, unname(c(tau - q["mu_N0"], tau / q["mu_N0"], tau - q["mu_C1"])))
  expect_equal(im$benchmark_value,
               unname(c(q["mu_C1"] - q["mu_C0"], q["mu_C1"] / q["mu_C0"], q["mu_N0"] - q["mu_C0"])))
  expect_true(is.na(im$relative[im$comparison == "ratio"]))

  # Ratios are not reported for outcomes that can be negative
  d_neg <- transform(sim_iv(p_at = 0, continuous = TRUE), y = y - 20)
  res_neg <- ate_sensitivity(d_neg, y ~ a | z, ci = FALSE)
  expect_false("ratio" %in% res_neg$implied$comparison)
})

test_that("two-sided boundary lies on ATE = cutoff in all parametrisations", {
  d <- sim_iv()
  q <- manual_q(d)
  res <- ate_sensitivity(d, y ~ a | z, cutoff = 0.05, ci = FALSE)
  expect_equal(res$setting, "two-sided")
  expect_null(res$threshold)
  b <- res$boundary
  expect_equal(nrow(b), 4)

  tA <- c(0.1, 0.6)
  lvl_N <- b[b$lhs == "theta_N1", ]
  expect_equal(ate_manual(q, tA[1], lvl_N$intercept + lvl_N$slope * tA[1]), 0.05)
  expect_equal(ate_manual(q, tA[2], lvl_N$intercept + lvl_N$slope * tA[2]), 0.05)

  tN <- 0.4
  lvl_A <- b[b$lhs == "theta_A0", ]
  expect_equal(ate_manual(q, lvl_A$intercept + lvl_A$slope * tN, tN), 0.05)

  dA <- 0.2
  eff_N <- b[b$lhs == "theta_N1 - mu_N0", ]
  tN_eff <- q[["mu_N0"]] + eff_N$intercept + eff_N$slope * dA
  expect_equal(ate_manual(q, q[["mu_A1"]] - dA, tN_eff), 0.05)

  expect_equal(b$relation, c(">", "<", ">", ">"))
})

test_that("restricting one parameter uses its least favourable value", {
  d <- sim_iv()
  q <- manual_q(d)

  res <- ate_sensitivity(d, y ~ a | z, theta_A0 = 1, ci = FALSE)
  expect_equal(res$fixed$value, 1)
  expect_equal(ate_manual(q, 1, res$threshold$estimate), 0)

  # Interval: least favourable theta_A0 is the upper end for ATE > x ...
  res_int <- ate_sensitivity(d, y ~ a | z, theta_A0 = c(0.2, 0.7), ci = FALSE)
  expect_equal(res_int$fixed$value, 0.7)
  # ... and the lower end for ATE < x
  res_less <- ate_sensitivity(d, y ~ a | z, theta_A0 = c(0.2, 0.7), direction = "less", ci = FALSE)
  expect_equal(res_less$fixed$value, 0.2)

  # Non-negative effect among always-takers: worst case theta_A0 = mu_A1
  res_eff <- ate_sensitivity(d, y ~ a | z, theta_A0 = restrict_effect(lower = 0), ci = FALSE)
  expect_equal(res_eff$fixed$value, q[["mu_A1"]])

  # Restricting theta_N1 gives a threshold on theta_A0
  res_N <- ate_sensitivity(d, y ~ a | z, theta_N1 = 0, ci = FALSE)
  expect_equal(res_N$threshold$parameter, "theta_A0")
  expect_equal(ate_manual(q, res_N$threshold$estimate, 0), 0)
})

test_that("invalid restriction combinations are rejected", {
  d <- sim_iv()
  expect_error(ate_sensitivity(d, y ~ a | z, theta_A0 = 1, theta_N1 = 0, ci = FALSE), "at most one")
  d_nt <- sim_iv(p_at = 0)
  expect_error(ate_sensitivity(d_nt, y ~ a | z, theta_N1 = 0, ci = FALSE), "sole sensitivity")
  expect_warning(ate_sensitivity(d_nt, y ~ a | z, theta_A0 = 1, ci = FALSE), "ignored")
  d_full <- transform(d, a = z)
  expect_error(ate_sensitivity(d_full, y ~ a | z, ci = FALSE), "point identified")
})

test_that("bootstrap confidence intervals are returned and cover the estimate", {
  set.seed(2)
  d <- sim_iv()
  res <- ate_sensitivity(d, y ~ a | z, theta_A0 = restrict_effect(lower = 0), n_boot = 200)
  expect_equal(dim(res$boot), c(200, 7))
  t <- res$threshold
  expect_true(t$ci_lower < t$estimate && t$estimate < t$ci_upper)
  # The least favourable value depends on mu_A1 and hence has a CI
  expect_true(res$fixed$ci_lower < res$fixed$ci_upper)
  expect_true(all(c("implied_ci_lower", "implied_ci_upper") %in% names(res$implied)))
  expect_true(all(c("intercept_ci_lower", "slope_ci_upper") %in% names(res$boundary)))
  expect_true(all(c("ci_lower", "ci_upper") %in% names(res$identified)))

  # Stratified bootstrap keeps the arm sizes fixed
  set.seed(3)
  res_s <- ate_sensitivity(d, y ~ a | z, n_boot = 50, boot_strata = "instrument")
  expect_equal(nrow(res_s$boot), 50)
})

test_that("bootstrap is reproducible with a seed and disabled with ci = FALSE", {
  d <- sim_iv(p_at = 0)
  set.seed(5); r1 <- ate_sensitivity(d, y ~ a | z, n_boot = 50)
  set.seed(5); r2 <- ate_sensitivity(d, y ~ a | z, n_boot = 50)
  expect_equal(r1$threshold, r2$threshold)
  r3 <- ate_sensitivity(d, y ~ a | z, ci = FALSE)
  expect_null(r3$boot)
  expect_false("ci_lower" %in% names(r3$threshold))
})

test_that("print methods run", {
  d <- sim_iv()
  expect_output(print(ate_sensitivity(d, y ~ a | z, ci = FALSE)), "Boundary")
  set.seed(1)
  expect_output(print(ate_sensitivity(d, y ~ a | z, theta_A0 = 1, n_boot = 50)), "Threshold")
})
