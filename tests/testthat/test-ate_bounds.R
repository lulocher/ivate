test_that("unrestricted bounds evaluate the ATE at the corners of the outcome range", {
  d <- sim_iv()
  q <- manual_q(d)
  res <- ate_bounds(d, y ~ a | z, ci = FALSE)
  expect_s3_class(res, "ivate_bounds")
  expect_equal(res$bounds$estimate, c(ate_manual(q, 1, 0), ate_manual(q, 0, 1)))
  expect_equal(res$bounds, res$unrestricted)

  d_c <- sim_iv(continuous = TRUE)
  q_c <- manual_q(d_c)
  r <- range(d_c$y)
  res_c <- ate_bounds(d_c, y ~ a | z, ci = FALSE)
  expect_equal(res_c$bounds$estimate, c(ate_manual(q_c, r[2], r[1]), ate_manual(q_c, r[1], r[2])))
  res_c2 <- ate_bounds(d_c, y ~ a | z, y_range = c(-20, 30), ci = FALSE)
  expect_equal(res_c2$bounds$estimate, c(ate_manual(q_c, 30, -20), ate_manual(q_c, -20, 30)))
})

test_that("fixing both parameters collapses the bounds to a point", {
  d <- sim_iv()
  q <- manual_q(d)
  res <- ate_bounds(d, y ~ a | z, theta_A0 = 0.4, theta_N1 = 0.3, ci = FALSE)
  expect_equal(res$bounds$estimate, rep(ate_manual(q, 0.4, 0.3), 2))
})

test_that("restrictions translate into the right intervals", {
  d <- sim_iv()
  q <- manual_q(d)
  late <- q[["mu_C1"]] - q[["mu_C0"]]
  feas <- function(...) ate_bounds(d, y ~ a | z, ..., ci = FALSE)$feasible

  f <- feas(theta_A0 = restrict_effect(0, 0.1), theta_N1 = restrict_effect(-0.1, 0))
  expect_equal(unlist(f[1, c("lower", "upper")], use.names = FALSE), q[["mu_A1"]] - c(0.1, 0))
  expect_equal(unlist(f[2, c("lower", "upper")], use.names = FALSE), q[["mu_N0"]] + c(-0.1, 0))

  f <- feas(theta_A0 = restrict_ratio(1, 2), theta_N1 = restrict_ratio(upper = 1.2))
  expect_equal(unlist(f[1, c("lower", "upper")], use.names = FALSE), q[["mu_A1"]] / c(2, 1))
  expect_equal(unlist(f[2, c("lower", "upper")], use.names = FALSE), c(0, 1.2 * q[["mu_N0"]]))

  f <- feas(theta_A0 = restrict_heterogeneity(-0.05, 0.05), theta_N1 = restrict_heterogeneity(upper = 0))
  expect_equal(unlist(f[1, c("lower", "upper")], use.names = FALSE), q[["mu_C0"]] + c(-0.05, 0.05))
  expect_equal(unlist(f[2, c("lower", "upper")], use.names = FALSE), c(0, q[["mu_C1"]]))

  # Formula and function limits are evaluated on the identified quantities
  f1 <- feas(theta_A0 = restrict_effect(lower = ~ late), theta_N1 = restrict_effect(upper = ~ late))
  f2 <- feas(theta_A0 = restrict_effect(lower = function(q) q$late),
             theta_N1 = restrict_effect(upper = function(q) q$late))
  expect_equal(f1[, c("lower", "upper")], f2[, c("lower", "upper")])
  expect_equal(f1$upper, c(q[["mu_A1"]] - late, q[["mu_N0"]] + late))

  # Lists of restrictions are intersected
  f <- feas(theta_N1 = list(c(0.1, 0.9), restrict_effect(upper = 0.05)))
  expect_equal(unlist(f[2, c("lower", "upper")], use.names = FALSE),
               c(0.1, min(0.9, q[["mu_N0"]] + 0.05)))
})

test_that("restricted bounds are nested in the unrestricted bounds", {
  d <- sim_iv()
  res <- ate_bounds(d, y ~ a | z, theta_A0 = restrict_effect(lower = 0),
                    theta_N1 = restrict_ratio(upper = 1.5), ci = FALSE)
  expect_gte(res$bounds$estimate[1], res$unrestricted$estimate[1])
  expect_lte(res$bounds$estimate[2], res$unrestricted$estimate[2])
})

test_that("incompatible or irrelevant restrictions are handled", {
  d <- sim_iv()
  expect_error(ate_bounds(d, y ~ a | z, theta_N1 = list(0.9, c(0, 0.5)), ci = FALSE), "incompatible")
  expect_error(restrict_effect(1, 0), "must not exceed")
  expect_error(restrict_ratio(-1, 2), ">= 0")
  expect_error(ate_bounds(d, y ~ a | z, theta_A0 = "high", ci = FALSE), "must be NULL")

  d_nt <- sim_iv(p_at = 0)
  expect_warning(res <- ate_bounds(d_nt, y ~ a | z, theta_A0 = 1, ci = FALSE), "ignored")
  expect_equal(res$feasible$parameter, "theta_N1")
})

test_that("bootstrap confidence intervals propagate data-dependent restrictions", {
  set.seed(4)
  d <- sim_iv()
  res <- ate_bounds(d, y ~ a | z, theta_A0 = restrict_effect(lower = ~ late), n_boot = 200)
  b <- res$bounds
  expect_true(all(b$ci_lower <= b$estimate & b$estimate <= b$ci_upper))
  expect_true(all(c("ci_lower", "ci_upper") %in% names(res$unrestricted)))
  expect_output(print(res), "Restricted")
})

test_that("restrictions print a readable description", {
  expect_equal(describe_restriction(restrict_ratio(upper = 1.2), "theta_N1"), "theta_N1 / mu_N0 <= 1.2")
  expect_equal(describe_restriction(restrict_effect(lower = ~ late), "theta_A0"), "late <= mu_A1 - theta_A0")
  expect_equal(describe_restriction(restrict_level(1, 1), "theta_A0"), "theta_A0 = 1")
  expect_output(print(restrict_effect(lower = 0)), "effect scale")
})
