# ivate

Informative inference on the average treatment effect (ATE) with instrumental
variables.

With a binary instrument and a binary treatment, the standard IV assumptions
point-identify the local ATE among compliers, but not the ATE. Assumption-free
(Balke–Pearl) bounds on the ATE are often too wide to answer the question a
decision maker cares about: *does the ATE exceed a policy-relevant cutoff
`x`?* Point identification of the ATE is precluded by exactly two unobserved
quantities:

- `theta_A0` = E[Y(0) | always-taker], and
- `theta_N1` = E[Y(1) | never-taker].

`ivate` implements the strategies of Locher, Richardson, and Stensrud (2026)
to reason about these two quantities directly:

| Function | Question it answers |
|---|---|
| `ate_sensitivity()` | How extreme must `theta_A0` / `theta_N1` be for the decision to flip? Thresholds, the boundary line in the two-sided case, and implied comparisons with identified quantities. |
| `ate_bounds()` | Which bounds on the ATE follow from restrictions I find credible (e.g. "no harm among always-takers")? |
| `ate_bayes()` | Posterior of the ATE as a function of, or under priors on, the sensitivity parameters (binary outcomes; Richardson, Evans, and Robins, 2011). |
| `plot_comparisons()` | Roster of estimated (green), assumed (orange), and implied (red) potential outcome means by compliance type, with the comparisons between them (as in the comparison figures of the manuscript). |

Uncertainty is quantified by the nonparametric bootstrap or, for binary
outcomes, by Bayesian inference (`inference = "bayes"`), which also reports a
posterior check of the IV inequalities.

## Installation

```r
# install.packages("remotes")
remotes::install_github("lulocher/ivate", build_vignettes = TRUE)
```

The vignette `vignette("replication", package = "ivate")` reproduces all
results in the examples of Locher, Richardson, and Stensrud (2026) from the
cell counts reported in the manuscript.

## Example

Simulated data from a trial with two-sided non-compliance: 50% compliers
(effect 0.15), 20% always-takers (effect 0.05), 30% never-takers (effect
0.02). The true ATE is 0.091; we ask whether it exceeds the cutoff 0.05.

```r
library(ivate)

set.seed(2026)
n    <- 5000
z    <- rbinom(n, 1, 0.5)
type <- sample(c("C", "A", "N"), n, replace = TRUE, prob = c(0.5, 0.2, 0.3))
a    <- ifelse(type == "A", 1, ifelse(type == "N", 0, z))
p0 <- c(C = 0.30, A = 0.55, N = 0.20)   # E[Y(0) | type]
p1 <- c(C = 0.45, A = 0.60, N = 0.22)   # E[Y(1) | type]
y0 <- rbinom(n, 1, p0[type])
y1 <- rbinom(n, 1, p1[type])
y  <- ifelse(a == 1, y1, y0)
dat <- data.frame(y, a, z)
```

**Without assumptions** the bounds contain the cutoff:

```r
ate_bounds(dat, y ~ a | z)
#>   Unrestricted: [-0.0673, 0.42]
```

**Sensitivity analysis.** Assume only that treatment does not harm
always-takers. Which values of `theta_N1` are compatible with ATE > 0.05?

```r
set.seed(1)
sens <- ate_sensitivity(dat, y ~ a | z, cutoff = 0.05,
                        theta_A0 = restrict_effect(lower = 0))
sens
#> Restriction: 0 <= mu_A1 - theta_A0
#>   least favourable value theta_A0 = 0.622 95% CI [0.577, 0.664]
#>
#> Threshold: the criterion holds iff theta_N1 > 0.146 95% CI [0.0549, 0.232]
#>
#> Implied comparisons at the threshold (implied vs. identified benchmark):
#>   effect         theta_N1 - mu_N0 = -0.0569 [-0.151, 0.0333]  vs  mu_C1 - mu_C0 = 0.129    (ratio -0.44)
#>   ratio          theta_N1 / mu_N0 = 0.719 [0.262, 1.18]       vs  mu_C1 / mu_C0 = 1.45
#>   heterogeneity  theta_N1 - mu_C1 = -0.272 [-0.393, -0.157]   vs  mu_N0 - mu_C0 = -0.0859  (ratio 3.17)
plot(sens)                    # ATE as a function of theta_N1
plot_comparisons(sens)        # roster of estimated / assumed / implied means (as in the example figures)
```

The ATE falls below 0.05 only if treatment *lowers* the outcome of
never-takers by more than 0.057, while it raises that of compliers by 0.129.
If such opposite-sign effects are implausible, the decision is resolved.

**Improved bounds** under restrictions stated on interpretable scales
(limits can depend on identified quantities such as the LATE):

```r
ate_bounds(dat, y ~ a | z,
           theta_A0 = restrict_effect(lower = 0),
           theta_N1 = restrict_effect(lower = 0, upper = ~ late))
#> Feasible sets:
#>   theta_A0 in [0, 0.622]       restriction: 0 <= mu_A1 - theta_A0
#>   theta_N1 in [0.202, 0.332]   restriction: 0 <= theta_N1 - mu_N0 <= late
#>
#> Bounds on the ATE:
#>   Restricted:   [0.0664, 0.228]
```

**Bayesian analysis** (binary outcomes): the posterior probability that the
criterion holds, as a function of the sensitivity parameters.

```r
set.seed(1)
ate_bayes(dat, y ~ a | z, cutoff = 0.05, theta_A0 = 0.6,
          theta_N1 = seq(0, 0.4, by = 0.05))
#>  theta_A0 theta_N1 median 95% CrI          P(ATE > 0.05)
#>  0.6      0.10     0.0408 [0.0174, 0.0652] 0.23
#>  0.6      0.15     0.0551 [0.0321, 0.0795] 0.67
#>  0.6      0.20     0.0695 [0.0464, 0.0936] 0.952
#>  ...
```

Passing priors instead of numbers (`prior_beta()`, `prior_uniform()`,
`prior_normal()`) gives a full posterior of the ATE. Since the data carry no
information about the sensitivity parameters, that posterior depends on
their prior even in large samples; `ate_bayes()` therefore reports the prior
of the ATE and the posterior of the identified bounds alongside it.

Published 2 × 2 × 2 tables can be analysed via `expand_counts()`.

## References

Locher, L., Richardson, T., and Stensrud, M. (2026). Informative inference on
the average treatment effect with instrumental variables. Working paper.

Richardson, T. S., Evans, R. J., and Robins, J. M. (2011). Transparent
parametrizations of models for potential outcomes. *Bayesian Statistics 9*,
569–610.

Balke, A. and Pearl, J. (1997). Bounds on treatment effects from studies with
imperfect compliance. *Journal of the American Statistical Association*, 92,
1171–1176.

Please cite the package with `citation("ivate")`.
