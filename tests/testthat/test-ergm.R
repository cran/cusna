# E6 smoke tests: the native ERGM path (R -> cpp11 -> vendored CPU host) is wired
# correctly. The heavy bit-identity / distributional validation lives in
# native/test (the C ABI is the same); these check shapes and hand values.

test_that("ergm_term validates term names", {
  expect_s3_class(ergm_term("edges"), "ergm_term")
  expect_s3_class(ergm_term("gwesp", 0.25), "ergm_term")
  expect_error(ergm_term("nope"), "unknown ERGM term")
  expect_error(ergm_term(c("a", "b")), "single ERGM term")
})

test_that("ergm_stats matches hand-computed edges/mutual (directed)", {
  # directed ties 1->2, 2->1, 2->3
  a <- matrix(0L, 3, 3)
  a[1, 2] <- 1L; a[2, 1] <- 1L; a[2, 3] <- 1L
  s <- ergm_stats(a, list(ergm_term("edges"), ergm_term("mutual")),
                  directed = TRUE)
  expect_equal(unname(s["edges"]), 3)    # 3 directed ties
  expect_equal(unname(s["mutual"]), 1)   # one mutual dyad (1<->2)
})

test_that("ergm_stats edges is half the entry count when undirected", {
  a <- matrix(0L, 4, 4)
  a[1, 2] <- a[2, 1] <- a[2, 3] <- a[3, 2] <- 1L   # 2 undirected edges
  s <- ergm_stats(a, list(ergm_term("edges")))
  expect_equal(unname(s["edges"]), 2)
})

test_that("ergm_simulate returns an nsim x k matrix of statistics", {
  a <- matrix(0L, 6, 6)
  a[1, 2] <- a[2, 1] <- a[2, 3] <- a[3, 2] <- 1L
  sf <- ergm_simulate(a, coef = -1, terms = list(ergm_term("edges")),
                      nsim = 8, burnin = 1000L, seed = 1L)
  expect_equal(dim(sf), c(8L, 1L))
  expect_equal(colnames(sf), "edges")
  expect_true(all(sf[, 1] >= 0))
})

test_that("ergm_simulate return_nets: stats recomputed from nets match (directed)", {
  # consistency: statistics recomputed from each returned network equal the row of
  # `stats` for that simulation, to machine precision.
  set.seed(101)
  n <- 8L
  a <- matrix(as.integer(runif(n * n) < 0.25), n, n); diag(a) <- 0L
  terms <- list(ergm_term("edges"), ergm_term("mutual"))
  sim <- ergm_simulate(a, coef = c(-1, 0.5), terms = terms, nsim = 6,
                       directed = TRUE, burnin = 500L, seed = 7L,
                       return_nets = TRUE)
  expect_named(sim, c("stats", "nets"))
  expect_length(sim$nets, 6L)
  expect_equal(dim(sim$nets[[1]]), c(n, n))
  for (c in seq_len(6L)) {
    rec <- ergm_stats(sim$nets[[c]], terms, directed = TRUE)
    expect_equal(unname(rec), unname(sim$stats[c, ]))
  }
})

test_that("ergm_simulate seed-in-seed: stats identical with/without return_nets", {
  # byte-for-byte regression: the default (scalar) call is unchanged, and the
  # stats returned alongside nets are the same as the bare stats matrix.
  a <- matrix(0L, 6, 6)
  a[1, 2] <- a[2, 1] <- a[2, 3] <- a[3, 2] <- 1L
  s_old <- ergm_simulate(a, coef = -1, terms = list(ergm_term("edges")),
                         nsim = 8, burnin = 1000L, seed = 1L)
  s_new <- ergm_simulate(a, coef = -1, terms = list(ergm_term("edges")),
                         nsim = 8, burnin = 1000L, seed = 1L, return_nets = TRUE)
  expect_identical(s_old, s_new$stats)
})

test_that("temporal terms are rejected in the cross-sectional path", {
  a <- matrix(0L, 4, 4); a[1, 2] <- a[2, 1] <- 1L
  expect_error(ergm_stats(a, list(ergm_term("memory")), directed = TRUE),
               "temporal term")
  expect_error(ergm_simulate(a, coef = 0, terms = list(ergm_term("memory")),
                             directed = TRUE),
               "temporal term")
})

test_that("ergm_mple returns 3 named pseudo-MLE coefficients", {
  set.seed(1)
  a <- matrix(as.integer(runif(18 * 18) < 0.25), 18, 18)
  diag(a) <- 0L
  co <- ergm_mple(a)
  expect_named(co, c("edges", "mutual", "transitive"))
  expect_length(co, 3)
  expect_true(all(is.finite(co)))
})

test_that("ergm_mcmle recovers the analytic edges-only MLE", {
  # for an edges-only ERGM (independent dyads) the MLE is exactly the logit of
  # the observed density -- an analytic check that needs no ergm reference.
  set.seed(2)
  n <- 12L
  a <- matrix(as.integer(runif(n * n) < 0.3), n, n)
  diag(a) <- 0L
  fit <- ergm_mcmle(a, list(ergm_term("edges")), directed = TRUE,
                    nsim = 300L, burnin = 3000L, maxit = 25L, seed = 3L)
  expect_s3_class(fit, "ergm_mcmle")
  expect_true(fit$converged)
  dens <- sum(a) / (n * (n - 1))
  expect_equal(unname(fit$coef[["edges"]]), stats::qlogis(dens), tolerance = 0.2)
})

test_that("tergm_mple fits a network sequence and recovers persistence", {
  # a persistent sequence: each wave keeps most ties of the previous one, so
  # the memory (lagged-tie) coefficient must be strongly positive.
  set.seed(4)
  n <- 15L
  base <- matrix(as.integer(runif(n * n) < 0.2), n, n)
  diag(base) <- 0L
  nets <- list(base)
  for (t in 2:3) {
    nx <- base
    flip <- matrix(runif(n * n) < 0.05, n, n)      # small perturbation
    nx[flip] <- 1L - nx[flip]
    diag(nx) <- 0L
    nets[[t]] <- nx
    base <- nx
  }
  fit <- tergm_mple(nets, ttriple = FALSE, R = 50L, seed = 5L)
  expect_s3_class(fit, "tergm_mple")
  expect_true("memory" %in% names(fit$coef))
  expect_gt(unname(fit$coef[["memory"]]), 1)       # strong persistence
  expect_equal(nrow(fit$ci), length(fit$coef))
})

# recompute a tergm term list's statistics from a network and the lag (cross terms
# via ergm_stats(directed = TRUE); memory = sum(x*lag), delrecip = sum(x*t(lag)))
.tergm_recompute <- function(x, lag, terms) {
  vapply(seq_along(terms), function(i) {
    nm <- terms[[i]]$name
    switch(nm,
           "memory"   = sum(x * lag),
           "delrecip" = sum(x * t(lag)),
           as.numeric(ergm_stats(x, terms[i], directed = TRUE)))
  }, numeric(1))
}

test_that("tergm_simulate seed-in-seed: legacy call is byte-for-byte unchanged", {
  set.seed(201)
  n <- 7L
  lag <- matrix(as.integer(runif(n * n) < 0.3), n, n); diag(lag) <- 0L
  x0  <- matrix(as.integer(runif(n * n) < 0.3), n, n); diag(x0) <- 0L
  # legacy: fixed flags, no `terms`/`attr`; must equal the historical assembly
  m_legacy <- tergm_simulate(x0, lag, coef = c(-1, 0.5, 1.0, 0.2),
                             mutual = TRUE, memory = TRUE, delrecip = TRUE,
                             nsim = 12, burnin = 800L, seed = 3L)
  expect_equal(colnames(m_legacy), c("edges", "mutual", "memory", "delrecip"))
  # same run through the explicit `terms` API must give identical stats
  terms <- list(ergm_term("edges"), ergm_term("mutual"),
                ergm_term("memory"), ergm_term("delrecip"))
  m_terms <- tergm_simulate(x0, lag, coef = c(-1, 0.5, 1.0, 0.2), terms = terms,
                            nsim = 12, burnin = 800L, seed = 3L)
  expect_identical(m_legacy, m_terms)
})

test_that("tergm_simulate return_nets: stats recomputed from nets match", {
  # consistency for a mixed cross-sectional + temporal term list
  set.seed(202)
  n <- 8L
  lag <- matrix(as.integer(runif(n * n) < 0.3), n, n); diag(lag) <- 0L
  x0  <- matrix(as.integer(runif(n * n) < 0.3), n, n); diag(x0) <- 0L
  terms <- list(ergm_term("edges"), ergm_term("mutual"),
                ergm_term("memory"), ergm_term("delrecip"))
  sim <- tergm_simulate(x0, lag, coef = c(-1, 0.4, 0.8, 0.3), terms = terms,
                        nsim = 6, burnin = 800L, seed = 5L, return_nets = TRUE)
  expect_named(sim, c("stats", "nets"))
  expect_length(sim$nets, 6L)
  for (c in seq_len(6L)) {
    rec <- .tergm_recompute(sim$nets[[c]], lag, terms)
    expect_equal(rec, unname(sim$stats[c, ]))
  }
})

test_that("tergm_simulate supports the nodematch cross-sectional term with attr", {
  set.seed(203)
  n <- 10L
  lag <- matrix(as.integer(runif(n * n) < 0.25), n, n); diag(lag) <- 0L
  x0  <- matrix(as.integer(runif(n * n) < 0.25), n, n); diag(x0) <- 0L
  grp <- rep(c(0, 1), length.out = n)
  terms <- list(ergm_term("edges"), ergm_term("nodematch"), ergm_term("memory"))
  sim <- tergm_simulate(x0, lag, coef = c(-1, 0.5, 1.0), terms = terms,
                        attr = grp, nsim = 5, burnin = 600L, seed = 6L,
                        return_nets = TRUE)
  expect_equal(colnames(sim$stats), c("edges", "nodematch", "memory"))
  for (c in seq_len(5L)) {
    x <- sim$nets[[c]]
    cross <- as.numeric(ergm_stats(x, terms[1:2], attr = grp, directed = TRUE))
    expect_equal(cross[1], unname(sim$stats[c, "edges"]))
    expect_equal(cross[2], unname(sim$stats[c, "nodematch"]))
    expect_equal(sum(x * lag), unname(sim$stats[c, "memory"]))
  }
})

test_that("alaam_mple equals glm on the hand-computed contagion (undirected)", {
  # the native contagion change statistic is A %*% y, so the ALAAM MPLE must be
  # exactly a logistic regression of y on that count + covariates -- an exact
  # self-consistency check that needs no external ALAAM package.
  set.seed(6)
  n <- 40L
  net <- matrix(as.integer(runif(n * n) < 0.12), n, n)
  net <- ((net + t(net)) > 0) * 1L                 # undirected
  diag(net) <- 0L
  x <- rnorm(n)
  y <- as.integer(runif(n) < plogis(-1 + 0.4 * (net %*% rep(1, n)) / 5 + 0.5 * x))
  fit <- alaam_mple(y, net, cov = x)
  con <- as.numeric(net %*% y)                     # hand-computed contagion
  ref <- glm(y ~ con + x, family = binomial())
  expect_equal(unname(fit$coef), unname(coef(ref)), tolerance = 1e-8)
  expect_named(fit$coef, c("density", "contagion", "cov"))
})

test_that("alaam_simulate returns attribute vectors and statistics", {
  set.seed(7)
  n <- 25L
  net <- matrix(as.integer(runif(n * n) < 0.15), n, n)
  net <- ((net + t(net)) > 0) * 1L
  diag(net) <- 0L
  sim <- alaam_simulate(net, coef = c(density = -0.5, contagion = 0.4),
                        nsim = 10L, sweeps = 50L, seed = 8L)
  expect_equal(dim(sim$y), c(10L, n))
  expect_true(all(sim$y %in% c(0L, 1L)))
  expect_equal(colnames(sim$stats), c("density", "contagion"))
  expect_true(all(sim$stats[, "density"] >= 0 & sim$stats[, "density"] <= n))
})

test_that("alaam_mcmle matches the moments of an attribute drawn from the model", {
  # generate an attribute from a known ALAAM, then check the MCMC-MLE drives the
  # simulated sufficient statistics onto the observed ones (moment matching).
  set.seed(9)
  n <- 30L
  net <- matrix(as.integer(runif(n * n) < 0.15), n, n)
  net <- ((net + t(net)) > 0) * 1L
  diag(net) <- 0L
  y <- alaam_simulate(net, c(density = -1.2, contagion = 0.5),
                      nsim = 1L, sweeps = 300L, seed = 11L)$y[1, ]
  fit <- alaam_mcmle(y, net, nsim = 400L, maxit = 30L, seed = 12L)
  expect_s3_class(fit, "alaam_mcmle")
  expect_true(fit$converged)
  expect_lt(max(abs(fit$tstat)), 0.3)     # observed moments recovered at the MLE
})

test_that("stergm_cmle recovers the analytic edges-only formation/persistence rates", {
  # for an edges-only model the constrained ERGM is dyad-independent, so the CMLE
  # of the edges coefficient is exactly the logit of the free-dyad on-rate
  # (formation: absent-at-t-1 dyads that form; persistence: present-at-t-1 dyads
  # that survive) -- an analytic check of the permitted-mask sampler + CMLE loop.
  set.seed(10)
  n <- 20L
  Y1 <- matrix(as.integer(runif(n * n) < 0.2), n, n); diag(Y1) <- 0L
  Y2 <- matrix(as.integer(runif(n * n) < 0.25), n, n); diag(Y2) <- 0L
  fit <- stergm_cmle(list(Y1, Y2), formation = list(ergm_term("edges")),
                     nsim = 400L, maxit = 25L, seed = 1L)
  expect_s3_class(fit, "stergm_cmle")
  free_f <- Y1 == 0; diag(free_f) <- FALSE
  free_d <- Y1 == 1; diag(free_d) <- FALSE
  rate_f <- sum(Y2[free_f]) / sum(free_f)
  rate_d <- sum(Y2[free_d]) / sum(free_d)
  expect_equal(unname(fit$formation$coef[["edges"]]),
               stats::qlogis(rate_f), tolerance = 0.3)
  expect_equal(unname(fit$dissolution$coef[["edges"]]),
               stats::qlogis(rate_d), tolerance = 0.3)
})
