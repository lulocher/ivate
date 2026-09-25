build_ok <- function(p) {
  expect_s3_class(p, "ggplot")
  expect_no_error(suppressWarnings(ggplot2::ggplot_build(p)))
}

test_that("threshold plot is drawn for one-sided non-compliance", {
  set.seed(1)
  res <- ate_sensitivity(sim_iv(p_at = 0), y ~ a | z, n_boot = 50)
  p <- plot(res)
  build_ok(p)
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true(all(c("GeomRibbon", "GeomLine", "GeomHline", "GeomVline") %in% geoms))
  build_ok(plot(res, range = c(0, 0.5), ci = FALSE))
  build_ok(ggplot2::autoplot(res))
  expect_error(plot(res, type = "boundary"), "two-sided")
})

test_that("boundary and restricted threshold plots are drawn for two-sided non-compliance", {
  set.seed(1)
  d <- sim_iv()
  res <- ate_sensitivity(d, y ~ a | z, n_boot = 50)
  p <- plot(res)
  build_ok(p)
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true(all(c("GeomRaster", "GeomRibbon", "GeomLine") %in% geoms))
  expect_error(plot(res, type = "threshold"), "No scalar threshold")

  res_fixed <- ate_sensitivity(d, y ~ a | z, theta_A0 = restrict_effect(lower = 0), n_boot = 50)
  build_ok(plot(res_fixed))
  build_ok(plot(res_fixed, type = "boundary"))
})

test_that("threshold outside the plotted range does not break the plot", {
  res <- ate_sensitivity(sim_iv(p_at = 0), y ~ a | z, cutoff = 0.5, ci = FALSE)
  build_ok(plot(res))
})

test_that("nearby reference labels are merged but ticks are kept", {
  out <- collapse_nearby_labels(c(0.5, 0.1, 0.105), c("a", "b", "c"), c(0, 1), tol = 0.03)
  expect_equal(out$vals, c(0.1, 0.105, 0.5))
  expect_equal(out$labels, c("b\nc", "", "a"))
})
