test_that("ate_bounds returns correct length-2 named vector", {
  set.seed(1)
  n   <- 1000
  z   <- rbinom(n, 1, 0.5)
  a   <- ifelse(z == 1, rbinom(n, 1, 0.7), rbinom(n, 1, 0.2))
  y   <- rbinom(n, 1, 0.3 + 0.2 * a)
  dat <- data.frame(y = y, a = a, z = z)

  res <- ate_bounds(dat, y ~ a | z)

  expect_length(res, 2)
  expect_named(res, c("Lower bound for ATE", "Upper bound for ATE"))
  expect_true(res["Lower bound for ATE"] <= res["Upper bound for ATE"])
})

test_that("ate_bounds lower <= upper for default [0,1] bounds", {
  set.seed(2)
  n   <- 500
  z   <- rbinom(n, 1, 0.5)
  a   <- ifelse(z == 1, rbinom(n, 1, 0.8), rbinom(n, 1, 0.1))
  y   <- rbinom(n, 1, 0.4 + 0.15 * a)
  dat <- data.frame(y = y, a = a, z = z)

  res <- ate_bounds(dat, y ~ a | z)
  expect_lte(res[1], res[2])
})

test_that("ate_bounds narrows with restricted sensitivity parameter ranges", {
  set.seed(3)
  n   <- 1000
  z   <- rbinom(n, 1, 0.5)
  a   <- ifelse(z == 1, rbinom(n, 1, 0.7), rbinom(n, 1, 0.2))
  y   <- rbinom(n, 1, 0.3 + 0.2 * a)
  dat <- data.frame(y = y, a = a, z = z)

  res_wide   <- ate_bounds(dat, y ~ a | z)
  res_narrow <- ate_bounds(dat, y ~ a | z,
                           bounds_y0_always_takers = c(0.2, 0.6),
                           bounds_y1_never_takers  = c(0.2, 0.6))

  width_wide   <- res_wide["Upper bound for ATE"]   - res_wide["Lower bound for ATE"]
  width_narrow <- res_narrow["Upper bound for ATE"] - res_narrow["Lower bound for ATE"]
  expect_lte(width_narrow, width_wide)
})

test_that("ate_bounds errors on non-binary treatment", {
  dat <- data.frame(y = c(0, 1, 0), a = c(0, 1, 2), z = c(0, 1, 0))
  expect_error(ate_bounds(dat, y ~ a | z), "'treatment' must be binary")
})

test_that("ate_bounds errors on non-binary instrument", {
  dat <- data.frame(y = c(0, 1, 0), a = c(0, 1, 0), z = c(0, 1, 2))
  expect_error(ate_bounds(dat, y ~ a | z), "'instrument' must be binary")
})

test_that("ate_bounds silently coerces integer treatment/instrument to numeric", {
  set.seed(4)
  n   <- 500
  z   <- as.integer(rbinom(n, 1, 0.5))
  a   <- as.integer(ifelse(z == 1, rbinom(n, 1, 0.7), rbinom(n, 1, 0.2)))
  y   <- rbinom(n, 1, 0.3 + 0.2 * a)
  dat <- data.frame(y = y, a = a, z = z)

  expect_no_error(ate_bounds(dat, y ~ a | z))
})
