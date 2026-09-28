# ivate 0.1.0

First public version, aligned with Locher, Richardson, and Stensrud (2026).

## Sensitivity analysis and bounds

* `ate_sensitivity()` returns thresholds (one-sided non-compliance), the
  boundary line on the level and effect scales (two-sided), and implied
  comparisons (effect, ratio, heterogeneity) with identified benchmarks.
  Restricting one sensitivity parameter (`theta_A0` / `theta_N1`) implements
  the two-step strategy of the manuscript. `direction = "less"` handles
  criteria of the form ATE < x.
* `ate_bounds()` accepts restrictions on the sensitivity parameters on the
  level, effect, ratio, or heterogeneity scale (`restrict_level()`,
  `restrict_effect()`, `restrict_ratio()`, `restrict_heterogeneity()`); limits
  may depend on identified quantities (e.g. `~ late`).
* The default outcome range is `c(0, 1)` for binary outcomes and the observed
  range otherwise.

## Inference

* Nonparametric bootstrap (optionally stratified by instrument) for all
  thresholds, boundary coefficients, implied comparisons, and bounds.
* Bayesian inference for binary outcomes (`inference = "bayes"`), following
  Richardson, Evans, and Robins (2011), with a posterior check of the IV
  inequalities.
* `ate_bayes()`: posterior of the ATE as a function of the sensitivity
  parameters, or under priors on them (`prior_beta()`, `prior_uniform()`,
  `prior_normal()`).

* Model check for plug-in inference: `ate_sensitivity()` and `ate_bounds()`
  check whether the plug-in estimates violate the IV model restrictions
  (pi_C > 0, complier means within the outcome range; inequalities (11)-(13)
  of Richardson, Evans, and Robins (2011) for binary outcomes) in the
  original sample and report the share of violating bootstrap samples
  (`$model_check`), warning if it exceeds 5%.

## Other

* `plot()` methods reproduce Figure 2 of the manuscript and Figures 3 and 5
  of Richardson, Evans, and Robins (2011).
* `plot_comparisons()` (or `plot(x, type = "comparisons")`) draws Figure 1
  of the manuscript with the values filled in: estimated quantities in
  green, assumed ones in orange, implied ones in red.
* `expand_counts()` converts published 2 x 2 x 2 tables into individual-level
  data.

## Breaking changes relative to the development version

* `decision_threshold` is now `cutoff`; `bounds_y0_always_takers` /
  `bounds_y1_never_takers` are now `theta_A0` / `theta_N1`; the `plot_*`
  arguments were replaced by a `plot()` method.
* Results are S3 objects with data-frame components instead of named vectors
  and character strings.
