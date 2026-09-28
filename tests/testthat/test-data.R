test_that("identified quantities match the definitions", {
  for (args in list(list(), list(p_at = 0), list(p_nt = 0), list(continuous = TRUE))) {
    d <- do.call(sim_iv, args)
    q <- identify_quantities(d$y, d$a, d$z)
    expect_equal(q, manual_q(d))
  }
})

test_that("iv_setting detects the compliance setting", {
  expect_equal(iv_setting(manual_q(sim_iv())), "two-sided")
  expect_equal(iv_setting(manual_q(sim_iv(p_at = 0))), "never-takers only")
  expect_equal(iv_setting(manual_q(sim_iv(p_nt = 0))), "always-takers only")
})

test_that("formula parsing accepts formulas and strings and rejects others", {
  expect_equal(unname(parse_iv_formula(y ~ a | z)), c("y", "a", "z"))
  expect_equal(unname(parse_iv_formula("y ~ a | z")), c("y", "a", "z"))
  expect_error(parse_iv_formula(y ~ a), "outcome ~ treatment \\| instrument")
  expect_error(parse_iv_formula(y ~ a + b | z), "single variable")
})

test_that("data validation catches common problems", {
  d <- sim_iv(500)
  expect_error(ate_bounds(d, y ~ a | w), "not found")
  expect_error(ate_bounds(transform(d, a = a * 2), y ~ a | z), "treatment 'a' must be binary")
  expect_error(ate_bounds(transform(d, z = z + 1), y ~ a | z), "instrument 'z' must be binary")
  expect_error(ate_bounds(transform(d, a = ifelse(a == 1, "yes", "no")), y ~ a | z),
               "cannot be coerced")
  expect_warning(ate_bounds(transform(d, z = 1 - z), y ~ a | z, ci = FALSE), "recode")
  expect_error(ate_bounds(d, y ~ a | z, y_range = c(0.2, 0.5), ci = FALSE), "outside")

  d_na <- d
  d_na$y[1:3] <- NA
  expect_warning(ate_bounds(d_na, y ~ a | z, ci = FALSE), "Removed 3")

  d_int <- transform(d, a = as.integer(a), z = factor(z), y = as.logical(y))
  expect_equal(ate_bounds(d_int, y ~ a | z, ci = FALSE)$bounds,
               ate_bounds(d, y ~ a | z, ci = FALSE)$bounds)
})

test_that("outcome range defaults to [0, 1] for binary and observed range otherwise", {
  expect_equal(resolve_y_range(c(0, 1, 1)), c(0, 1))
  expect_equal(resolve_y_range(c(2, 5, 3)), c(2, 5))
  expect_equal(resolve_y_range(c(2, 5, 3), c(0, 10)), c(0, 10))
})
