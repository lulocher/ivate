# Simulate IV data with always-takers (share p_at), never-takers (p_nt), and
# compliers. The instrument is randomised; outcomes depend on type and
# treatment so that all identified means differ.
sim_iv <- function(n = 3000, p_at = 0.2, p_nt = 0.3, continuous = FALSE, seed = 1) {
  set.seed(seed)
  z    <- rbinom(n, 1, 0.5)
  type <- sample(c("C", "A", "N"), n, replace = TRUE, prob = c(1 - p_at - p_nt, p_at, p_nt))
  a    <- ifelse(type == "A", 1, ifelse(type == "N", 0, z))
  base <- c(C = 0.3, A = 0.5, N = 0.2)[type]
  eff  <- c(C = 0.2, A = 0.1, N = 0.05)[type]
  mu   <- base + eff * a
  y    <- if (continuous) rnorm(n, 10 * mu) else rbinom(n, 1, mu)
  data.frame(y = y, a = a, z = z)
}

# Identified quantities computed directly from the definitions in Appendix A.
manual_q <- function(d) {
  with(d, {
    p_at <- mean(a[z == 0]); p_nt <- 1 - mean(a[z == 1]); p_c <- 1 - p_at - p_nt
    y1_at <- if (p_at > 0) mean(y[a == 1 & z == 0]) else NA
    y0_nt <- if (p_nt > 0) mean(y[a == 0 & z == 1]) else NA
    y1_c <- (mean(a[z == 1]) * mean(y[a == 1 & z == 1]) - ifelse(p_at > 0, p_at * y1_at, 0)) / p_c
    y0_c <- ((1 - mean(a[z == 0])) * mean(y[a == 0 & z == 0]) - ifelse(p_nt > 0, p_nt * y0_nt, 0)) / p_c
    c(pi_C = p_c, pi_A = p_at, pi_N = p_nt, mu_C0 = y0_c, mu_C1 = y1_c, mu_A1 = y1_at, mu_N0 = y0_nt)
  })
}

ate_manual <- function(q, tA, tN) {
  unname(q["pi_C"] * (q["mu_C1"] - q["mu_C0"]) +
           ifelse(q["pi_A"] > 0, q["pi_A"] * (q["mu_A1"] - tA), 0) +
           ifelse(q["pi_N"] > 0, q["pi_N"] * (tN - q["mu_N0"]), 0))
}
