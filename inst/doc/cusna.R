library(cusna)
cusna_abi_version()
cusna_has_cuda()      # FALSE on the CPU-only build

# Keep the vignette build well within CRAN's two-core policy: the native CPU
# backend is OpenMP-parallel over simulation chains, so run it single-threaded
# here. cusna_set_threads() calls the OpenMP runtime directly, so it takes
# effect immediately (unlike Sys.setenv(OMP_NUM_THREADS = ...), whose effect
# on an already-initialised OpenMP runtime is platform-dependent). The demo
# data is tiny, so a single thread is still fast.
old_threads <- cusna_set_threads(1L)

effects <- list(
  cusna_effect("density"),
  cusna_effect("recip"),
  cusna_effect("transTrip")
)
effects[[1]]

alcohol <- c(3, 2, 2, 1, 3)         # one value per actor (toy)
cusna_effect("egoX", covariate = alcohol)

list(
  cusna_effect("egoX", covariate = alcohol),
  cusna_effect("altX", covariate = alcohol),
  cusna_interaction(c(1, 2))          # egoX x altX
)

set.seed(7)
n <- 20
w1 <- matrix(as.integer(runif(n * n) < 0.12), n, n); diag(w1) <- 0L
w2 <- w1; flip <- sample(n * n, 40); w2[flip] <- 1L - w2[flip]; diag(w2) <- 0L
w3 <- w2; flip <- sample(n * n, 40); w3[flip] <- 1L - w3[flip]; diag(w3) <- 0L

dat <- saom_data(list(w1, w2, w3))
dat

fit <- mom_estimate(
  dat,
  effects = list(cusna_effect("density"), cusna_effect("recip")),
  conditional = TRUE,                 # RSiena's default for one network
  control = mom_control(n1 = 100, nsub = 1, n2 = 10, batch2 = 50, n3 = 200))
fit

summary(fit)          # RSiena-style table (rates, estimates, s.e., t-conv)
coef(fit)             # named parameter vector
round(vcov(fit), 4)   # delta-method covariance

dat <- saom_data(s50$waves, behavior = s50$alcohol)    # the bundled s50 panel

fit <- mom_estimate(
  dat,
  effects = list(cusna_effect("density"), cusna_effect("recip"),
                 cusna_effect("transTrip"),
                 cusna_effect("altX", dyn = TRUE),      # covariate = the behavior
                 cusna_effect("egoX", dyn = TRUE)),
  beh_effects = list(cusna_beh_effect("linear"), cusna_beh_effect("quad"),
                     cusna_beh_effect("avSim")),
  conditional = FALSE)

fit_gof <- mom_estimate(
  dat,
  effects = list(cusna_effect("density"), cusna_effect("recip")),
  conditional = TRUE, return_nets = TRUE,
  control = mom_control(n1 = 100, nsub = 1, n2 = 10, batch2 = 50, n3 = 200))
sims <- fit_gof$sims[[2]]                       # period 2 ends at wave 3
dim(sims)
sim_indeg <- cusna_gof_distribution(sims, kind = 0L, levls = 0:4)
obs_indeg <- cusna_gof_distribution(array(w3, c(1, n, n)), kind = 0L, levls = 0:4)
rbind(observed = obs_indeg[1, ], simulated = colMeans(sim_indeg))

str(mom_control())

x <- w1                                       # any 0/1 adjacency
terms <- list(ergm_term("edges"), ergm_term("mutual"))
ergm_stats(x, terms, directed = TRUE)

ergm_simulate(x, coef = c(-2, 1), terms = terms, nsim = 5,
              directed = TRUE, seed = 1)

ergm_mcmle(x, terms, directed = TRUE)       # MCMC-MLE       (matches ergm::ergm)
ergm_mple(x)                                # pseudo-likelihood demo
tergm_mple(list(w1, w2, w3))                # TERGM           (matches btergm)
tergm_simulate(w2, lag = w1, coef = c(-4, 3, 2, 0.5))
stergm_cmle(list(w1, w2, w3),               # separable TERGM (matches tergm CMLE)
            formation = list(ergm_term("edges"), ergm_term("mutual")))
alaam_mple(y, net = x)                      # ALAAM MPLE      (matches glm)
alaam_mcmle(y, net = x)                     # ALAAM MCMC-MLE

