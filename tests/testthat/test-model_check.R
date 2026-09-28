test_that("restriction violations are detected per restriction", {
  mk <- function(pi_C, mu_C0, mu_C1) {
    c(pi_C = pi_C, pi_A = 0.1, pi_N = 0.6, mu_C0 = mu_C0, mu_C1 = mu_C1, mu_A1 = 0.5, mu_N0 = 0.2,
      pimu_C0 = pi_C * mu_C0, pimu_C1 = pi_C * mu_C1)
  }
  q_ok  <- mk(0.3, 0.4, 0.6)
  q_bad <- mk(0.3, -0.1, 0.5)
  q_pic <- mk(-0.01, 0.4, 0.6)
  v <- iv_restriction_violations(rbind(q_ok, q_bad, q_pic), c(0, 1))
  expect_equal(unname(v[1, ]), c(FALSE, FALSE, FALSE))
  expect_equal(unname(v[2, ]), c(FALSE, TRUE, FALSE))
  expect_true(v[3, "pi_C"])
  # A vector is treated as a single sample; the range is respected
  expect_equal(unname(iv_restriction_violations(q_ok, c(0, 0.5))[1, ]), c(FALSE, FALSE, TRUE))
})

test_that("well-behaved data pass the model check in the sample and the bootstrap", {
  set.seed(1)
  res <- ate_sensitivity(sim_iv(5000), y ~ a | z, n_boot = 200)
  m <- res$model_check
  expect_false(any(m$sample))
  expect_equal(m$bootstrap[["any"]], 0)
  expect_equal(m$n_boot, 200)
  expect_output(print(res), "Model check: IV restrictions satisfied in the sample; violated in 0% of 200")

  b <- ate_bounds(sim_iv(5000), y ~ a | z, ci = FALSE)
  expect_null(b$model_check$bootstrap)
  expect_output(print(b), "satisfied in the sample\n")
})

test_that("violations in the sample and frequent bootstrap violations are reported", {
  # Flu vaccine data: the plug-in complier mean under treatment exceeds 1
  flu <- expand_counts(flu_counts)
  set.seed(2)
  expect_warning(
    expect_warning(res <- ate_sensitivity(flu, y ~ a | z, n_boot = 200), "violate the IV model restrictions"),
    "bootstrap samples violate")
  expect_true(res$model_check$sample[["mu_C1"]])
  expect_gt(res$model_check$bootstrap[["mu_C1"]], 0.5)
  expect_output(print(res), "VIOLATED in the sample: mu_C1 in \\[0, 1\\]")

  # Weak instrument and rare outcome: plug-in fine, but bootstrap violations frequent
  set.seed(3)
  n <- 400
  z <- rbinom(n, 1, 0.5)
  type <- sample(c("C", "N"), n, replace = TRUE, prob = c(0.08, 0.92))
  a <- ifelse(type == "C", z, 0)
  y <- rbinom(n, 1, 0.03)
  d <- data.frame(y, a, z)
  res_w <- suppressWarnings(ate_bounds(d, y ~ a | z, n_boot = 200))
  expect_gt(res_w$model_check$bootstrap[["any"]], 0.05)
})

test_that("Bayesian results keep their own model check", {
  set.seed(4)
  res <- ate_sensitivity(expand_counts(lipid_counts), y ~ a | z, inference = "bayes", n_draws = 200)
  expect_null(res$model_check$bootstrap)
  out <- capture.output(print(res))
  expect_true(any(grepl("P\\(IV inequalities violated \\| data\\)", out)))
  expect_false(any(grepl("bootstrap samples", out)))
})
