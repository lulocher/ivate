test_that("ate_sensitivity returns correct structure for two-sided non-compliance", {
  set.seed(10)
  n   <- 1000
  z   <- rbinom(n, 1, 0.5)
  a   <- ifelse(z == 1, rbinom(n, 1, 0.7), rbinom(n, 1, 0.2))
  y   <- rbinom(n, 1, 0.3 + 0.2 * a)
  dat <- data.frame(y = y, a = a, z = z)

  res <- ate_sensitivity(dat, y ~ a | z)

  expect_named(res, c("identified_quantities", "threshold", "plot"))
  expect_type(res$threshold, "character")
  expect_s3_class(res$plot, "ggplot")
})

test_that("ate_sensitivity returns correct structure for never-takers only", {
  set.seed(11)
  n   <- 1000
  z   <- rbinom(n, 1, 0.5)
  a   <- ifelse(z == 1, rbinom(n, 1, 0.6), 0L)   # no always-takers
  y   <- rbinom(n, 1, 0.3 + 0.2 * a)
  dat <- data.frame(y = y, a = a, z = z)

  res <- ate_sensitivity(dat, y ~ a | z)

  expect_named(res, c("identified_quantities", "threshold", "plot"))
  expect_type(res$threshold, "list")
  expect_true("theta_N" %in% names(res$threshold))
  expect_true("implied_effect_additive" %in% names(res$threshold))
  expect_s3_class(res$plot, "ggplot")
})

test_that("ate_sensitivity returns correct structure for always-takers only", {
  set.seed(12)
  n   <- 1000
  z   <- rbinom(n, 1, 0.5)
  a   <- ifelse(z == 1, 1L, rbinom(n, 1, 0.2))   # no never-takers
  y   <- rbinom(n, 1, 0.3 + 0.2 * a)
  dat <- data.frame(y = y, a = a, z = z)

  res <- ate_sensitivity(dat, y ~ a | z)

  expect_named(res, c("identified_quantities", "threshold", "plot"))
  expect_type(res$threshold, "list")
  expect_true("theta_A" %in% names(res$threshold))
  expect_true("implied_effect_additive" %in% names(res$threshold))
  expect_s3_class(res$plot, "ggplot")
})

test_that("ate_sensitivity always returns element named 'threshold'", {
  set.seed(13)
  n   <- 1000
  z   <- rbinom(n, 1, 0.5)
  dat_two  <- data.frame(
    y = rbinom(n, 1, 0.4),
    a = ifelse(z == 1, rbinom(n, 1, 0.7), rbinom(n, 1, 0.2)),
    z = z
  )
  dat_nt <- data.frame(
    y = rbinom(n, 1, 0.4),
    a = ifelse(z == 1, rbinom(n, 1, 0.6), 0L),
    z = z
  )
  dat_at <- data.frame(
    y = rbinom(n, 1, 0.4),
    a = ifelse(z == 1, 1L, rbinom(n, 1, 0.2)),
    z = z
  )

  for (dat in list(dat_two, dat_nt, dat_at)) {
    res <- ate_sensitivity(dat, y ~ a | z)
    expect_true("threshold" %in% names(res))
  }
})

test_that("ate_sensitivity respects decision_threshold argument", {
  set.seed(14)
  n   <- 1000
  z   <- rbinom(n, 1, 0.5)
  a   <- ifelse(z == 1, rbinom(n, 1, 0.6), 0L)
  y   <- rbinom(n, 1, 0.3 + 0.2 * a)
  dat <- data.frame(y = y, a = a, z = z)

  res0  <- ate_sensitivity(dat, y ~ a | z, decision_threshold = 0)
  res01 <- ate_sensitivity(dat, y ~ a | z, decision_threshold = 0.1)

  expect_false(identical(res0$threshold, res01$threshold))
})

test_that("ate_sensitivity returns NULL invisibly for full compliance", {
  set.seed(15)
  n   <- 500
  z   <- rbinom(n, 1, 0.5)
  a   <- z                        # perfect compliance, no AT or NT
  y   <- rbinom(n, 1, 0.3 + 0.2 * a)
  dat <- data.frame(y = y, a = a, z = z)

  expect_message(
    res <- ate_sensitivity(dat, y ~ a | z),
    "point identified"
  )
  expect_null(res)
})

test_that("ate_sensitivity errors on non-binary treatment or instrument", {
  dat <- data.frame(y = c(0, 1, 0, 1), a = c(0, 1, 2, 1), z = c(0, 1, 0, 1))
  expect_error(ate_sensitivity(dat, y ~ a | z), "'treatment' must be binary")
})

test_that("ate_sensitivity silently coerces integer inputs", {
  set.seed(16)
  n   <- 500
  z   <- as.integer(rbinom(n, 1, 0.5))
  a   <- as.integer(ifelse(z == 1, rbinom(n, 1, 0.6), 0L))
  y   <- rnorm(n)
  dat <- data.frame(y = y, a = a, z = z)

  expect_no_error(ate_sensitivity(dat, y ~ a | z))
})
