# Test bug fixes: predict RE indexing, GP reshaping, type case,
# ignore.RE, residuals with multiple chains, ppcOcc single-visit,
# single-factor w.0.samples dimensions.

skip_on_cran()

# predict with >= 3 unstructured random effects ---------------------------
# The RE column offset into beta.star.samples must be cumulative across
# all preceding REs. Predicting at the fitted sites with the fitted
# covariates and RE levels must reproduce the fitted psi.samples exactly.
set.seed(100)
J.x <- 8
J.y <- 8
J <- J.x * J.y
n.rep <- rep(3, J)
beta <- c(0.3, 0.5)
alpha <- c(-0.2)
psi.RE <- list(levels = c(10, 5, 8),
	       sigma.sq.psi = c(1.2, 0.8, 1.5))
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
	      psi.RE = psi.RE, p.RE = list(), sp = FALSE)
data.list <- list(y = dat$y,
		  occ.covs = data.frame(occ.cov = dat$X[, 2],
					f.1 = dat$X.re[, 1],
					f.2 = dat$X.re[, 2],
					f.3 = dat$X.re[, 3]),
		  det.covs = list(int = matrix(1, J, max(n.rep))))
out <- PGOcc(occ.formula = ~ occ.cov + (1 | f.1) + (1 | f.2) + (1 | f.3),
	     det.formula = ~ 1,
	     data = data.list,
	     n.samples = 400,
	     n.burn = 200,
	     n.thin = 2,
	     n.chains = 1,
	     verbose = FALSE)
X.0 <- cbind(1, dat$X[, 2], dat$X.re)
colnames(X.0) <- c('(Intercept)', 'occ.cov', 'f.1', 'f.2', 'f.3')
test_that("predict recovers fitted psi with three occurrence REs", {
  pred.out <- predict(out, X.0)
  expect_equal(pred.out$psi.0.samples, out$psi.samples,
	       check.attributes = FALSE, tolerance = 1e-10)
})

# residuals with multiple chains ------------------------------------------
test_that("residuals works when n.post.samples exceeds the per-chain count", {
  out.4 <- PGOcc(occ.formula = ~ occ.cov + (1 | f.1) + (1 | f.2) + (1 | f.3),
		 det.formula = ~ 1,
		 data = data.list,
		 n.samples = 300,
		 n.burn = 220,
		 n.thin = 1,
		 n.chains = 4,
		 verbose = FALSE)
  # 80 saved samples per chain, 320 total.
  resid.out <- residuals(out.4, n.post.samples = 100)
  expect_equal(dim(resid.out$occ.resids), c(100, J))
})

# spMsPGOcc GP predictions -------------------------------------------------
# The GP (NNGP = FALSE) branch must not reshape the prediction output
# twice. Predicting at the fitted coordinates must reproduce psi.samples.
set.seed(101)
N <- 3
beta <- cbind(rnorm(N, 0.3, 0.5), rnorm(N, 0.4, 0.5))
alpha <- cbind(rnorm(N, 0.2, 0.4))
dat.ms <- simMsOcc(J.x = 6, J.y = 6, n.rep = rep(3, 36), N = N,
		   beta = beta, alpha = alpha, sp = TRUE,
		   cov.model = 'exponential', sigma.sq = rep(1, N),
		   phi = rep(3 / 0.5, N))
data.list.ms <- list(y = dat.ms$y,
		     occ.covs = data.frame(occ.cov = dat.ms$X[, 2]),
		     det.covs = list(int = matrix(1, 36, 3)),
		     coords = dat.ms$coords)
out.gp <- spMsPGOcc(occ.formula = ~ occ.cov,
		    det.formula = ~ 1,
		    data = data.list.ms,
		    n.batch = 20,
		    batch.length = 25,
		    n.burn = 400,
		    n.thin = 1,
		    n.chains = 1,
		    NNGP = FALSE,
		    verbose = FALSE)
X.0.ms <- cbind(1, data.list.ms$occ.covs$occ.cov)
colnames(X.0.ms) <- c('(Intercept)', 'occ.cov')
test_that("predict.spMsPGOcc GP branch recovers fitted psi", {
  pred.gp <- predict(out.gp, X.0.ms, coords.0 = data.list.ms$coords,
		     verbose = FALSE)
  expect_equal(pred.gp$psi.0.samples, out.gp$psi.samples,
	       check.attributes = FALSE, tolerance = 1e-10)
})

test_that("predict accepts case-insensitive type", {
  pred.occ <- predict(out.gp, X.0.ms, coords.0 = data.list.ms$coords,
		      type = 'Occupancy', verbose = FALSE)
  expect_equal(dim(pred.occ$psi.0.samples), dim(out.gp$psi.samples))
})

# lfMsPGOcc ignore.RE ------------------------------------------------------
set.seed(102)
dat.lf <- simMsOcc(J.x = 6, J.y = 6, n.rep = rep(3, 36), N = N,
		   beta = beta, alpha = alpha, sp = FALSE,
		   psi.RE = list(levels = c(8), sigma.sq.psi = c(1)))
data.list.lf <- list(y = dat.lf$y,
		     occ.covs = data.frame(occ.cov = dat.lf$X[, 2],
					   site.re = dat.lf$X.re[, 1]),
		     det.covs = list(int = matrix(1, 36, 3)),
		     coords = dat.lf$coords)
out.lf <- lfMsPGOcc(occ.formula = ~ occ.cov + (1 | site.re),
		    det.formula = ~ 1,
		    data = data.list.lf,
		    n.factors = 2,
		    n.samples = 400,
		    n.burn = 300,
		    n.thin = 1,
		    n.chains = 1,
		    verbose = FALSE)
test_that("predict.lfMsPGOcc ignore.RE = TRUE works without RE columns", {
  X.0.lf <- cbind(1, data.list.lf$occ.covs$occ.cov)
  colnames(X.0.lf) <- c('(Intercept)', 'occ.cov')
  pred.lf <- predict(out.lf, X.0.lf, coords.0 = data.list.lf$coords + 0.013,
		     ignore.RE = TRUE)
  expect_equal(dim(pred.lf$psi.0.samples), c(out.lf$n.post * out.lf$n.chains,
					     N, 36))
})

# ppcOcc with a single replicate per season --------------------------------
set.seed(103)
J.t <- 25
dat.t <- simTOcc(J.x = 5, J.y = 5, n.time = rep(4, J.t),
		 n.rep = matrix(1, J.t, 4),
		 beta = c(0.4, 0.5), alpha = c(0.1), trend = TRUE, sp = FALSE)
data.list.t <- list(y = dat.t$y,
		    occ.covs = list(trend = dat.t$X[, , 2]),
		    det.covs = list(int = dat.t$X.p[, , , 1]))
out.t <- tPGOcc(occ.formula = ~ trend,
		det.formula = ~ 1,
		data = data.list.t,
		n.batch = 20,
		batch.length = 25,
		n.burn = 250,
		n.thin = 1,
		n.chains = 1,
		verbose = FALSE)
test_that("ppcOcc works on multi-season fits with one visit per season", {
  ppc.1 <- suppressMessages(ppcOcc(out.t, fit.stat = 'chi-squared', group = 1))
  ppc.2 <- suppressMessages(ppcOcc(out.t, fit.stat = 'freeman-tukey', group = 2))
  expect_equal(dim(ppc.1$fit.y), c(out.t$n.post * out.t$n.chains, 4))
  expect_equal(dim(ppc.2$fit.y), c(out.t$n.post * out.t$n.chains, 4))
})

# stMsPGOcc predictions with a single factor --------------------------------
set.seed(104)
n.time.ms <- 3
dat.st <- simTMsOcc(J.x = 5, J.y = 5, n.time = rep(n.time.ms, 25),
		    n.rep = matrix(2, 25, n.time.ms), N = N,
		    beta = beta, alpha = alpha, trend = TRUE, sp = TRUE,
		    cov.model = 'exponential', sigma.sq = rep(1, N),
		    phi = 3 / 0.5, factor.model = TRUE, n.factors = 1)
data.list.st <- list(y = dat.st$y,
		     occ.covs = list(trend = dat.st$X[, , 2]),
		     det.covs = list(int = array(1, dim = c(25, n.time.ms, 2))),
		     coords = dat.st$coords)
out.st <- stMsPGOcc(occ.formula = ~ trend,
		    det.formula = ~ 1,
		    data = data.list.st,
		    n.factors = 1,
		    n.batch = 20,
		    batch.length = 25,
		    n.burn = 400,
		    n.thin = 1,
		    n.chains = 1,
		    n.neighbors = 5,
		    verbose = FALSE)
test_that("predict.stMsPGOcc keeps the factor dimension when q = 1", {
  X.0.st <- array(NA, dim = c(25, n.time.ms, 2))
  X.0.st[, , 1] <- 1
  X.0.st[, , 2] <- data.list.st$occ.covs$trend
  pred.st <- predict(out.st, X.0.st, coords.0 = dat.st$coords + 0.013,
		     t.cols = 1:n.time.ms, verbose = FALSE)
  expect_equal(length(dim(pred.st$w.0.samples)), 3)
  expect_equal(dim(pred.st$w.0.samples)[2], 1)
})
