# Test fpPGOcc.R  --------------------------------------------------------

skip_on_cran()

# Intercept Only ----------------------------------------------------------
set.seed(100)
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
gamma <- c(1.5)
p.tp.det <- length(gamma)
delta <- c(0.5)
p.fp.det <- length(delta)
psi.RE <- list()
mu.RE <- list()
omega.RE <- list()
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = FALSE, prop.val = 0.2)
y <- dat$y
X <- dat$X
X.mu <- dat$X.mu
X.omega <- dat$X.omega

data.list <- list(y = y, n.val = dat$n.val, d = dat$d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0),
                   fix = TRUE)
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3)

n.batch <- 40
batch.length <- 25
n.report <- 100

out <- fpPGOcc(occ.formula = ~ 1,
               tp.det.formula = ~ 1,
               fp.det.formula = ~ 1,
               data = data.list,
               inits = inits.list,
               priors = prior.list,
               tuning = tuning.list,
               n.batch = n.batch,
               batch.length = batch.length,
               n.omp.threads = 1,
               verbose = FALSE,
               n.report = n.report,
               n.burn = 400,
               n.thin = 6,
               n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class fpPGOcc", {
  expect_s3_class(out, "fpPGOcc")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, FALSE)
  expect_equal(out$muRE, FALSE)
  expect_equal(out$omegaRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
  expect_equal(out$n.val, dat$n.val)
  expect_equal(out$d, dat$d)
})

# Check output dimensions -------------
test_that("posterior sample dimensions are correct", {
  n.post.samples <- out$n.post * out$n.chains
  expect_equal(n.post.samples, 200)
  expect_equal(dim(out$beta.samples), c(n.post.samples, p.occ))
  expect_equal(dim(out$gamma.samples), c(n.post.samples, p.tp.det))
  expect_equal(dim(out$delta.samples), c(n.post.samples, p.fp.det))
  expect_equal(dim(out$z.samples), c(n.post.samples, J))
  expect_equal(dim(out$psi.samples), c(n.post.samples, J))
  expect_equal(length(out$rhat$beta), p.occ)
  expect_equal(length(out$rhat$gamma), p.tp.det)
  expect_equal(length(out$rhat$delta), p.fp.det)
})

# Check fixed initial values ----------
test_that("validated true positives force z = 1", {
  z.fixed <- which(apply(dat$d, 1, sum, na.rm = TRUE) > 0)
  expect_true(all(out$z.samples[, z.fixed] == 1))
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- fpPGOcc(occ.formula = ~ 1,
                 tp.det.formula = ~ 1,
                 fp.det.formula = ~ 1,
                 data = data.list,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 n.omp.threads = 1,
                 verbose = FALSE)
  expect_s3_class(out, "fpPGOcc")
})
# Check non-integer n.post -------------
test_that("non-integer n.post", {
  expect_error(out <- fpPGOcc(occ.formula = ~ 1,
                              tp.det.formula = ~ 1,
                              fp.det.formula = ~ 1,
                              data = data.list,
                              n.thin = 13,
                              n.batch = n.batch,
                              batch.length = batch.length,
                              n.omp.threads = 1,
                              verbose = FALSE))
})
# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(fpPGOcc(occ.formula = ~ 1,
                        tp.det.formula = ~ 1,
                        fp.det.formula = ~ 1,
                        data = data.list,
                        n.batch = 4,
                        batch.length = 25,
                        n.omp.threads = 1,
                        verbose = TRUE,
                        n.report = n.report,
                        n.burn = 1,
                        n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for fpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$mu.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$omega.samples), c(n.post.samples, dim(y)))
  # Missing cells in y stay missing
  expect_true(all(is.na(fitted.out$y.rep.samples[1, , ][is.na(y)])))
  # Replicate validated true positives cannot exceed validated signals
  expect_true(all(sweep(fitted.out$d.rep.samples, c(2, 3), dat$n.val, '<=') |
                  is.na(fitted.out$d.rep.samples)))
})

# Check predictions -------------------
test_that("predict works for fpPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("true positive prediction works", {
  J.str <- 100
  X.mu.0 <- matrix(1, nrow = J.str, ncol = p.tp.det)
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J.str))
  expect_true(all(pred.out$mu.0.samples > 0))
})
test_that("false positive prediction works", {
  J.str <- 100
  X.omega.0 <- matrix(1, nrow = J.str, ncol = p.fp.det)
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, J.str))
  expect_true(all(pred.out$omega.0.samples > 0))
})
test_that("predict errors on bad input", {
  expect_error(predict(out, dat$X, type = 'detection'))
  expect_error(predict(out, matrix(1, 10, 3)))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(length(ppc.out$fit.d.rep), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.rep.group.quants), c(5, dim(y)))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(length(ppc.out$fit.d.rep), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
    expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, J))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J))
    expect_equal(dim(ppc.out$fit.d.rep.group.quants), c(5, J))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.y.rep), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(length(ppc.out$fit.d.rep), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.y.rep.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.rep.group.quants), c(5, max(n.rep)))
  }
  expect_output(summary(ppc.out))
  expect_error(ppcOcc(out, 'chi-square', 3))
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
  expect_equal(dim(out.resids$occ.resids), c(10, J))
  expect_equal(dim(out.resids$det.resids), c(10, dim(y)))
  expect_equal(dim(out.resids$val.resids), c(10, dim(y)))
  for (type in c('pearson', 'response')) {
    out.resids <- residuals(out, n.post.samples = 10, type = type)
    expect_equal(dim(out.resids$det.resids), c(10, dim(y)))
  }
  expect_error(residuals(out, n.post.samples = 10, type = 'deviance'))
  expect_error(residuals(out, n.post.samples = out$n.post * out$n.chains + 1))
})

# Check plot --------------------------
test_that("plot works for fpPGOcc", {
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(out, 'beta', density = FALSE))
  expect_no_error(plot(out, 'gamma', density = FALSE))
  expect_no_error(plot(out, 'delta', density = FALSE))
  expect_error(plot(out, 'beta.comm'))
})

# Occurrence covariate only -----------------------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.5)
p.occ <- length(beta)
gamma <- c(1.5)
p.tp.det <- length(gamma)
delta <- c(0.5)
p.fp.det <- length(delta)
psi.RE <- list()
mu.RE <- list()
omega.RE <- list()
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = FALSE, prop.val = 0.2)
y <- dat$y
X <- dat$X
X.mu <- dat$X.mu
X.omega <- dat$X.omega

colnames(X) <- c('int', 'occ.cov.1')

data.list <- list(y = y, occ.covs = X, n.val = dat$n.val, d = dat$d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3)

out <- fpPGOcc(occ.formula = ~ occ.cov.1,
               tp.det.formula = ~ 1,
               fp.det.formula = ~ 1,
               data = data.list,
               inits = inits.list,
               priors = prior.list,
               tuning = tuning.list,
               n.batch = n.batch,
               batch.length = batch.length,
               n.omp.threads = 1,
               verbose = FALSE,
               n.report = n.report,
               n.burn = 400,
               n.thin = 6,
               n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class fpPGOcc", {
  expect_s3_class(out, "fpPGOcc")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, FALSE)
  expect_equal(out$muRE, FALSE)
  expect_equal(out$omegaRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- fpPGOcc(occ.formula = ~ occ.cov.1,
                 tp.det.formula = ~ 1,
                 fp.det.formula = ~ 1,
                 data = data.list,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 n.chains = 1,
                 n.omp.threads = 1,
                 verbose = FALSE)
  expect_s3_class(out, "fpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(fpPGOcc(occ.formula = ~ occ.cov.1,
                       tp.det.formula = ~ 1,
                       fp.det.formula = ~ 1,
                       data = tmp.data,
                       n.batch = n.batch,
                       batch.length = batch.length,
                       n.omp.threads = 1,
                       verbose = FALSE))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(fpPGOcc(occ.formula = ~ occ.cov.1,
                        tp.det.formula = ~ 1,
                        fp.det.formula = ~ 1,
                        data = data.list,
                        n.batch = 4,
                        batch.length = 25,
                        n.omp.threads = 1,
                        verbose = TRUE,
                        n.report = n.report,
                        n.burn = 1,
                        n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for fpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
})

# Check predictions -------------------
test_that("predict works for fpPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("true positive prediction works", {
  J.str <- 100
  X.mu.0 <- matrix(1, nrow = J.str, ncol = p.tp.det)
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J.str))
})
test_that("false positive prediction works", {
  J.str <- 100
  X.omega.0 <- matrix(1, nrow = J.str, ncol = p.fp.det)
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, J.str))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, dim(y)))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, max(n.rep)))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
})

# True positive covariate only --------------------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
gamma <- c(1.5, -0.4)
p.tp.det <- length(gamma)
delta <- c(0.5)
p.fp.det <- length(delta)
psi.RE <- list()
mu.RE <- list()
omega.RE <- list()
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = FALSE, prop.val = 0.2)
y <- dat$y
X <- dat$X
X.mu <- dat$X.mu
X.omega <- dat$X.omega

data.list <- list(y = y, tp.det.covs = list(tp.cov.1 = X.mu[, , 2]),
                  n.val = dat$n.val, d = dat$d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3)
occ.formula <- ~ 1
tp.det.formula <- ~ tp.cov.1
fp.det.formula <- ~ 1

out <- fpPGOcc(occ.formula = occ.formula,
               tp.det.formula = tp.det.formula,
               fp.det.formula = fp.det.formula,
               data = data.list,
               inits = inits.list,
               priors = prior.list,
               tuning = tuning.list,
               n.batch = n.batch,
               batch.length = batch.length,
               n.omp.threads = 1,
               verbose = FALSE,
               n.report = n.report,
               n.burn = 400,
               n.thin = 6,
               n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class fpPGOcc", {
  expect_s3_class(out, "fpPGOcc")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, FALSE)
  expect_equal(out$muRE, FALSE)
  expect_equal(out$omegaRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- fpPGOcc(occ.formula = occ.formula,
                 tp.det.formula = tp.det.formula,
                 fp.det.formula = fp.det.formula,
                 data = data.list,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 n.omp.threads = 1,
                 verbose = FALSE)
  expect_s3_class(out, "fpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  # Missing covariate where y is observed
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1
  tmp.data$n.val[1, 1] <- 0
  tmp.data$d[1, 1] <- 0
  tmp.data$tp.det.covs$tp.cov.1[1, 1] <- NA
  expect_error(fpPGOcc(occ.formula = occ.formula,
                       tp.det.formula = tp.det.formula,
                       fp.det.formula = fp.det.formula,
                       data = tmp.data,
                       n.batch = n.batch,
                       batch.length = batch.length,
                       n.omp.threads = 1,
                       verbose = FALSE))
  # Missing y without matching missing n.val and d
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  expect_error(fpPGOcc(occ.formula = occ.formula,
                       tp.det.formula = tp.det.formula,
                       fp.det.formula = fp.det.formula,
                       data = tmp.data,
                       n.batch = n.batch,
                       batch.length = batch.length,
                       n.omp.threads = 1,
                       verbose = FALSE))
  # Missing y, n.val, and d in the same cell
  tmp.data$n.val[1, 1] <- NA
  tmp.data$d[1, 1] <- NA
  out <- fpPGOcc(occ.formula = occ.formula,
                 tp.det.formula = tp.det.formula,
                 fp.det.formula = fp.det.formula,
                 data = tmp.data,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 n.omp.threads = 1,
                 verbose = FALSE)
  expect_s3_class(out, "fpPGOcc")
  expect_equal(nrow(out$X.mu), sum(!is.na(tmp.data$y)))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(fpPGOcc(occ.formula = occ.formula,
                        tp.det.formula = tp.det.formula,
                        fp.det.formula = fp.det.formula,
                        data = data.list,
                        n.batch = 4,
                        batch.length = 25,
                        n.omp.threads = 1,
                        verbose = TRUE,
                        n.report = n.report,
                        n.burn = 1,
                        n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for fpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
})

# Check predictions -------------------
test_that("predict works for fpPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("true positive prediction works", {
  X.mu.0 <- dat$X.mu[, 1, ]
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("false positive prediction works", {
  J.str <- 100
  X.omega.0 <- matrix(1, nrow = J.str, ncol = p.fp.det)
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, J.str))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, dim(y)))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, max(n.rep)))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
})

# False positive covariate only -------------------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
gamma <- c(1.5)
p.tp.det <- length(gamma)
delta <- c(0.5, 0.6)
p.fp.det <- length(delta)
psi.RE <- list()
mu.RE <- list()
omega.RE <- list()
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = FALSE, prop.val = 0.2)
y <- dat$y
X <- dat$X
X.mu <- dat$X.mu
X.omega <- dat$X.omega

data.list <- list(y = y, fp.det.covs = list(fp.cov.1 = X.omega[, , 2]),
                  n.val = dat$n.val, d = dat$d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3)
occ.formula <- ~ 1
tp.det.formula <- ~ 1
fp.det.formula <- ~ fp.cov.1

out <- fpPGOcc(occ.formula = occ.formula,
               tp.det.formula = tp.det.formula,
               fp.det.formula = fp.det.formula,
               data = data.list,
               inits = inits.list,
               priors = prior.list,
               tuning = tuning.list,
               n.batch = n.batch,
               batch.length = batch.length,
               n.omp.threads = 1,
               verbose = FALSE,
               n.report = n.report,
               n.burn = 400,
               n.thin = 6,
               n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class fpPGOcc", {
  expect_s3_class(out, "fpPGOcc")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, FALSE)
  expect_equal(out$muRE, FALSE)
  expect_equal(out$omegaRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- fpPGOcc(occ.formula = occ.formula,
                 tp.det.formula = tp.det.formula,
                 fp.det.formula = fp.det.formula,
                 data = data.list,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 n.omp.threads = 1,
                 verbose = FALSE)
  expect_s3_class(out, "fpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1
  tmp.data$n.val[1, 1] <- 0
  tmp.data$d[1, 1] <- 0
  tmp.data$fp.det.covs$fp.cov.1[1, 1] <- NA
  expect_error(fpPGOcc(occ.formula = occ.formula,
                       tp.det.formula = tp.det.formula,
                       fp.det.formula = fp.det.formula,
                       data = tmp.data,
                       n.batch = n.batch,
                       batch.length = batch.length,
                       n.omp.threads = 1,
                       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  tmp.data$n.val[1, 1] <- NA
  tmp.data$d[1, 1] <- NA
  out <- fpPGOcc(occ.formula = occ.formula,
                 tp.det.formula = tp.det.formula,
                 fp.det.formula = fp.det.formula,
                 data = tmp.data,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 n.omp.threads = 1,
                 verbose = FALSE)
  expect_s3_class(out, "fpPGOcc")
  expect_equal(nrow(out$X.omega), sum(!is.na(tmp.data$y)))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(fpPGOcc(occ.formula = occ.formula,
                        tp.det.formula = tp.det.formula,
                        fp.det.formula = fp.det.formula,
                        data = data.list,
                        n.batch = 4,
                        batch.length = 25,
                        n.omp.threads = 1,
                        verbose = TRUE,
                        n.report = n.report,
                        n.burn = 1,
                        n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for fpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
})

# Check predictions -------------------
test_that("predict works for fpPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("true positive prediction works", {
  J.str <- 100
  X.mu.0 <- matrix(1, nrow = J.str, ncol = p.tp.det)
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J.str))
})
test_that("false positive prediction works", {
  X.omega.0 <- dat$X.omega[, 1, ]
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, dim(y)))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, max(n.rep)))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
})

# Covariates on all -------------------------------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.5)
p.occ <- length(beta)
gamma <- c(1.5, -0.4)
p.tp.det <- length(gamma)
delta <- c(0.5, 0.6, -0.3)
p.fp.det <- length(delta)
psi.RE <- list()
mu.RE <- list()
omega.RE <- list()
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = FALSE, prop.val = 0.2)
y <- dat$y
X <- dat$X
X.mu <- dat$X.mu
X.omega <- dat$X.omega

colnames(X) <- c('int', 'occ.cov.1')
data.list <- list(y = y, occ.covs = X,
                  tp.det.covs = list(tp.cov.1 = X.mu[, , 2]),
                  fp.det.covs = list(fp.cov.1 = X.omega[, , 2],
                                     fp.cov.2 = X.omega[, , 3]),
                  n.val = dat$n.val, d = dat$d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3)
occ.formula <- ~ occ.cov.1
tp.det.formula <- ~ tp.cov.1
fp.det.formula <- ~ fp.cov.1 + fp.cov.2

out <- fpPGOcc(occ.formula = occ.formula,
               tp.det.formula = tp.det.formula,
               fp.det.formula = fp.det.formula,
               data = data.list,
               inits = inits.list,
               priors = prior.list,
               tuning = tuning.list,
               n.batch = n.batch,
               batch.length = batch.length,
               n.omp.threads = 1,
               verbose = FALSE,
               n.report = n.report,
               n.burn = 400,
               n.thin = 6,
               n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class fpPGOcc", {
  expect_s3_class(out, "fpPGOcc")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, FALSE)
  expect_equal(out$muRE, FALSE)
  expect_equal(out$omegaRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check parameter names ---------------
test_that("parameter names are correct", {
  expect_equal(colnames(out$beta.samples), c('(Intercept)', 'occ.cov.1'))
  expect_equal(colnames(out$gamma.samples), c('(Intercept)', 'tp.cov.1'))
  expect_equal(colnames(out$delta.samples), c('(Intercept)', 'fp.cov.1', 'fp.cov.2'))
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- fpPGOcc(occ.formula = occ.formula,
                 tp.det.formula = tp.det.formula,
                 fp.det.formula = fp.det.formula,
                 data = data.list,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 n.omp.threads = 1,
                 verbose = FALSE)
  expect_s3_class(out, "fpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$occ.covs[3, ] <- NA
  expect_error(fpPGOcc(occ.formula = occ.formula,
                       tp.det.formula = tp.det.formula,
                       fp.det.formula = fp.det.formula,
                       data = tmp.data,
                       n.batch = n.batch,
                       batch.length = batch.length,
                       n.omp.threads = 1,
                       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- 1
  tmp.data$n.val[1, 1] <- 0
  tmp.data$d[1, 1] <- 0
  tmp.data$fp.det.covs$fp.cov.2[1, 1] <- NA
  expect_error(fpPGOcc(occ.formula = occ.formula,
                       tp.det.formula = tp.det.formula,
                       fp.det.formula = fp.det.formula,
                       data = tmp.data,
                       n.batch = n.batch,
                       batch.length = batch.length,
                       n.omp.threads = 1,
                       verbose = FALSE))
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  tmp.data$n.val[1, 1] <- NA
  tmp.data$d[1, 1] <- NA
  out <- fpPGOcc(occ.formula = occ.formula,
                 tp.det.formula = tp.det.formula,
                 fp.det.formula = fp.det.formula,
                 data = tmp.data,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 n.omp.threads = 1,
                 verbose = FALSE)
  expect_s3_class(out, "fpPGOcc")
})

# Check validation data ---------------
test_that("validation data error handling works", {
  fit.fp <- function(data, ...) {
    fpPGOcc(occ.formula = occ.formula,
            tp.det.formula = tp.det.formula,
            fp.det.formula = fp.det.formula,
            data = data,
            n.batch = n.batch,
            batch.length = batch.length,
            n.omp.threads = 1,
            verbose = FALSE, ...)
  }
  # n.val and d must be supplied
  tmp.data <- data.list
  tmp.data$n.val <- NULL
  expect_error(fit.fp(tmp.data), "n.val must be specified")
  tmp.data <- data.list
  tmp.data$d <- NULL
  expect_error(fit.fp(tmp.data), "d must be specified")
  # Dimensions must match y
  tmp.data <- data.list
  tmp.data$n.val <- tmp.data$n.val[, -1]
  expect_error(fit.fp(tmp.data), "dimension mismatch")
  tmp.data <- data.list
  tmp.data$d <- tmp.data$d[-1, ]
  expect_error(fit.fp(tmp.data), "dimension mismatch")
  # Missingness pattern must match y
  tmp.data <- data.list
  obs.indx <- which(!is.na(tmp.data$y))[1]
  tmp.data$d[obs.indx] <- NA
  expect_error(fit.fp(tmp.data), "mismatch between missing values")
  # Counts must be non-negative integers
  tmp.data <- data.list
  tmp.data$y[obs.indx] <- 1.5
  expect_error(fit.fp(tmp.data), "non-negative integer")
  tmp.data <- data.list
  tmp.data$n.val[obs.indx] <- -1
  expect_error(fit.fp(tmp.data), "non-negative integer")
  # d <= n.val <= y
  tmp.data <- data.list
  tmp.data$n.val[obs.indx] <- tmp.data$y[obs.indx] + 1
  expect_error(fit.fp(tmp.data), "n.val cannot be larger than y")
  tmp.data <- data.list
  tmp.data$d[obs.indx] <- tmp.data$n.val[obs.indx] + 1
  expect_error(fit.fp(tmp.data), "d cannot be larger than n.val")
  # z inits must be 1 where a validated true positive exists
  z.inits <- rep(0, J)
  expect_error(fit.fp(data.list, inits = list(z = z.inits)),
               "initial values for latent occurrence")
  # Missing formulas and data
  expect_error(fpPGOcc(occ.formula = occ.formula, fp.det.formula = fp.det.formula,
                       data = data.list, n.batch = n.batch,
                       batch.length = batch.length, verbose = FALSE),
               "tp.det.formula must be specified")
  expect_error(fpPGOcc(occ.formula = occ.formula, tp.det.formula = tp.det.formula,
                       data = data.list, n.batch = n.batch,
                       batch.length = batch.length, verbose = FALSE),
               "fp.det.formula must be specified")
  tmp.data <- data.list
  tmp.data$fp.det.covs <- NULL
  expect_error(fit.fp(tmp.data), "fp.det.covs must be specified")
})

# Check priors, inits, and tuning -----
test_that("prior, initial value, and tuning error handling works", {
  fit.fp <- function(...) {
    fpPGOcc(occ.formula = occ.formula,
            tp.det.formula = tp.det.formula,
            fp.det.formula = fp.det.formula,
            data = data.list,
            n.batch = n.batch,
            batch.length = batch.length,
            n.omp.threads = 1,
            verbose = FALSE, ...)
  }
  expect_error(fit.fp(priors = list(gamma.normal = list(0, rep(10, 3)))))
  expect_error(fit.fp(priors = list(delta.normal = list(rep(0, 2), 10))))
  expect_error(fit.fp(priors = list(delta.normal = 10)))
  expect_error(fit.fp(inits = list(gamma = c(0, 0, 0))))
  expect_error(fit.fp(inits = list(delta = c(0, 0))))
  expect_error(fit.fp(inits = list(z = rep(2, J))))
  expect_error(fit.fp(tuning = list(gamma = 0.3)))
  expect_error(fit.fp(tuning = list(delta = 0.3)))
  expect_error(fit.fp(tuning = list(gamma = c(0.3, 0.3, 0.3), delta = 0.3)))
  # Vector-valued priors, inits, and tuning work
  out <- fit.fp(priors = list(gamma.normal = list(rep(0, p.tp.det), rep(10, p.tp.det)),
                              delta.normal = list(rep(0, p.fp.det), rep(10, p.fp.det))),
                inits = list(gamma = rep(0, p.tp.det), delta = rep(0, p.fp.det)),
                tuning = list(gamma = rep(0.3, p.tp.det), delta = rep(0.3, p.fp.det)))
  expect_s3_class(out, "fpPGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(fpPGOcc(occ.formula = occ.formula,
                        tp.det.formula = tp.det.formula,
                        fp.det.formula = fp.det.formula,
                        data = data.list,
                        n.batch = 4,
                        batch.length = 25,
                        n.omp.threads = 1,
                        verbose = TRUE,
                        n.report = n.report,
                        n.burn = 1,
                        n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for fpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
})

# Check predictions -------------------
test_that("predict works for fpPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("true positive prediction works", {
  X.mu.0 <- dat$X.mu[, 1, ]
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("false positive prediction works", {
  X.omega.0 <- dat$X.omega[, 1, ]
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, J))
  expect_error(predict(out, dat$X.omega[, 1, 1:2], type = 'false positive'))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, dim(y)))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, max(n.rep)))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
})

# Interactions in all -----------------------------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.5, -0.4)
p.occ <- length(beta)
gamma <- c(1.5, -0.4, 0.2)
p.tp.det <- length(gamma)
delta <- c(0.5, 0.6, -0.3)
p.fp.det <- length(delta)
psi.RE <- list()
mu.RE <- list()
omega.RE <- list()
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = FALSE, prop.val = 0.2)
y <- dat$y
X <- dat$X
X.mu <- dat$X.mu
X.omega <- dat$X.omega

colnames(X) <- c('int', 'occ.cov.1', 'occ.cov.2')
data.list <- list(y = y, occ.covs = X,
                  tp.det.covs = list(tp.cov.1 = X.mu[, , 2],
                                     tp.cov.2 = X.mu[, , 3]),
                  fp.det.covs = list(fp.cov.1 = X.omega[, , 2],
                                     fp.cov.2 = X.omega[, , 3]),
                  n.val = dat$n.val, d = dat$d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3)
occ.formula <- ~ occ.cov.1 * occ.cov.2
tp.det.formula <- ~ tp.cov.1 * tp.cov.2
fp.det.formula <- ~ fp.cov.1 * fp.cov.2

out <- fpPGOcc(occ.formula = occ.formula,
               tp.det.formula = tp.det.formula,
               fp.det.formula = fp.det.formula,
               data = data.list,
               inits = inits.list,
               priors = prior.list,
               tuning = tuning.list,
               n.batch = n.batch,
               batch.length = batch.length,
               n.omp.threads = 1,
               verbose = FALSE,
               n.report = n.report,
               n.burn = 400,
               n.thin = 6,
               n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class fpPGOcc", {
  expect_s3_class(out, "fpPGOcc")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, FALSE)
  expect_equal(out$muRE, FALSE)
  expect_equal(out$omegaRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check interaction terms -------------
test_that("interaction terms are included", {
  expect_equal(ncol(out$beta.samples), 4)
  expect_equal(ncol(out$gamma.samples), 4)
  expect_equal(ncol(out$delta.samples), 4)
  expect_equal(colnames(out$gamma.samples)[4], 'tp.cov.1:tp.cov.2')
  expect_equal(colnames(out$delta.samples)[4], 'fp.cov.1:fp.cov.2')
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- fpPGOcc(occ.formula = occ.formula,
                 tp.det.formula = tp.det.formula,
                 fp.det.formula = fp.det.formula,
                 data = data.list,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 n.omp.threads = 1,
                 verbose = FALSE)
  expect_s3_class(out, "fpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(fpPGOcc(occ.formula = occ.formula,
                        tp.det.formula = tp.det.formula,
                        fp.det.formula = fp.det.formula,
                        data = data.list,
                        n.batch = 4,
                        batch.length = 25,
                        n.omp.threads = 1,
                        verbose = TRUE,
                        n.report = n.report,
                        n.burn = 1,
                        n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for fpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
})

# Check predictions -------------------
test_that("predict works for fpPGOcc", {
  X.0 <- cbind(dat$X, dat$X[, 2] * dat$X[, 3])
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("true positive prediction works", {
  X.mu.0 <- cbind(dat$X.mu[, 1, ], dat$X.mu[, 1, 2] * dat$X.mu[, 1, 3])
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("false positive prediction works", {
  X.omega.0 <- cbind(dat$X.omega[, 1, ], dat$X.omega[, 1, 2] * dat$X.omega[, 1, 3])
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, dim(y)))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, max(n.rep)))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
})

# Site covariate on detection ---------------------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.5)
p.occ <- length(beta)
gamma <- c(1.5)
p.tp.det <- 2
delta <- c(0.5)
p.fp.det <- 2
psi.RE <- list()
mu.RE <- list()
omega.RE <- list()
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = FALSE, prop.val = 0.2)
y <- dat$y
X <- dat$X
X.mu <- dat$X.mu
X.omega <- dat$X.omega

colnames(X) <- c('int', 'occ.cov.1')
# Site-level detection covariates supplied as length J vectors
data.list <- list(y = y, occ.covs = X,
                  tp.det.covs = list(site.cov = X[, 2]),
                  fp.det.covs = list(site.cov = X[, 2]),
                  n.val = dat$n.val, d = dat$d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3)
occ.formula <- ~ occ.cov.1
tp.det.formula <- ~ site.cov
fp.det.formula <- ~ site.cov

out <- fpPGOcc(occ.formula = occ.formula,
               tp.det.formula = tp.det.formula,
               fp.det.formula = fp.det.formula,
               data = data.list,
               inits = inits.list,
               priors = prior.list,
               tuning = tuning.list,
               n.batch = n.batch,
               batch.length = batch.length,
               n.omp.threads = 1,
               verbose = FALSE,
               n.report = n.report,
               n.burn = 400,
               n.thin = 6,
               n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class fpPGOcc", {
  expect_s3_class(out, "fpPGOcc")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, FALSE)
  expect_equal(out$muRE, FALSE)
  expect_equal(out$omegaRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check site covariate expansion ------
test_that("site-level detection covariates are expanded to all surveys", {
  n.obs <- sum(!is.na(y))
  expect_equal(nrow(out$X.mu), n.obs)
  expect_equal(nrow(out$X.omega), n.obs)
  site.long <- rep(1:J, ncol(y))[!is.na(c(y))]
  expect_equal(unname(out$X.mu[, 2]), unname(X[site.long, 2]))
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- fpPGOcc(occ.formula = occ.formula,
                 tp.det.formula = tp.det.formula,
                 fp.det.formula = fp.det.formula,
                 data = data.list,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 n.omp.threads = 1,
                 verbose = FALSE)
  expect_s3_class(out, "fpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(fpPGOcc(occ.formula = occ.formula,
                        tp.det.formula = tp.det.formula,
                        fp.det.formula = fp.det.formula,
                        data = data.list,
                        n.batch = 4,
                        batch.length = 25,
                        n.omp.threads = 1,
                        verbose = TRUE,
                        n.report = n.report,
                        n.burn = 1,
                        n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for fpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
})

# Check predictions -------------------
test_that("predict works for fpPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("true positive prediction works", {
  X.mu.0 <- dat$X
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("false positive prediction works", {
  X.omega.0 <- dat$X
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, dim(y)))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, max(n.rep)))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
})

# Random intercept on occurrence ------------------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
gamma <- c(1.5)
p.tp.det <- length(gamma)
delta <- c(0.5)
p.fp.det <- length(delta)
psi.RE <- list(levels = c(45),
               sigma.sq.psi = c(1.3))
mu.RE <- list()
omega.RE <- list()
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = FALSE, prop.val = 0.2)
y <- dat$y
X <- dat$X
X.re <- dat$X.re + 15
X.mu <- dat$X.mu
X.omega <- dat$X.omega

occ.covs <- cbind(X, X.re)
colnames(occ.covs) <- c('int', 'occ.factor.1')
data.list <- list(y = y, occ.covs = occ.covs, n.val = dat$n.val, d = dat$d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10),
                   sigma.sq.psi.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0, sigma.sq.psi = 1,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3)
occ.formula <- ~ (1 | occ.factor.1)
tp.det.formula <- ~ 1
fp.det.formula <- ~ 1

out <- fpPGOcc(occ.formula = occ.formula,
               tp.det.formula = tp.det.formula,
               fp.det.formula = fp.det.formula,
               data = data.list,
               inits = inits.list,
               priors = prior.list,
               tuning = tuning.list,
               n.batch = n.batch,
               batch.length = batch.length,
               n.omp.threads = 1,
               verbose = FALSE,
               n.report = n.report,
               n.burn = 400,
               n.thin = 6,
               n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class fpPGOcc", {
  expect_s3_class(out, "fpPGOcc")
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$psiRE, TRUE)
  expect_equal(out$muRE, FALSE)
  expect_equal(out$omegaRE, FALSE)
  expect_equal(dim(out$sigma.sq.psi.samples), c(out$n.post * out$n.chains, 1))
  expect_equal(ncol(out$beta.star.samples), length(unique(c(X.re))))
  expect_equal(length(out$rhat$sigma.sq.psi), 1)
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$re.level.names))),
               sort(unique(c(X.re))))
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$occ.covs <- as.data.frame(data.list$occ.covs)
  data.list$occ.covs$occ.factor.1 <- factor(data.list$occ.covs$occ.factor.1)
  expect_error(out <- fpPGOcc(occ.formula = occ.formula,
                              tp.det.formula = tp.det.formula,
                              fp.det.formula = fp.det.formula,
                              data = data.list,
                              n.batch = n.batch,
                              batch.length = batch.length,
                              n.omp.threads = 1,
                              verbose = FALSE))
  data.list$occ.covs$occ.factor.1 <- as.character(factor(data.list$occ.covs$occ.factor.1))
  expect_error(out <- fpPGOcc(occ.formula = occ.formula,
                              tp.det.formula = tp.det.formula,
                              fp.det.formula = fp.det.formula,
                              data = data.list,
                              n.batch = n.batch,
                              batch.length = batch.length,
                              n.omp.threads = 1,
                              verbose = FALSE))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- fpPGOcc(occ.formula = occ.formula,
                 tp.det.formula = tp.det.formula,
                 fp.det.formula = fp.det.formula,
                 data = data.list,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 n.omp.threads = 1,
                 verbose = FALSE)
  expect_s3_class(out, "fpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(fpPGOcc(occ.formula = occ.formula,
                        tp.det.formula = tp.det.formula,
                        fp.det.formula = fp.det.formula,
                        data = data.list,
                        n.batch = 4,
                        batch.length = 25,
                        n.omp.threads = 1,
                        verbose = TRUE,
                        n.report = n.report,
                        n.burn = 1,
                        n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for fpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
})

# Check predictions -------------------
test_that("predict works for fpPGOcc", {
  pred.out <- predict(out, occ.covs)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
  # Ignoring the random effects
  pred.out <- predict(out, occ.covs[, 1, drop = FALSE], ignore.RE = TRUE)
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  # New random effect levels
  X.0.new <- occ.covs
  X.0.new[, 'occ.factor.1'] <- max(X.re) + 1
  pred.out <- predict(out, X.0.new)
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  # Column names must match
  X.0.bad <- occ.covs
  colnames(X.0.bad) <- c('int', 'wrong.name')
  expect_error(predict(out, X.0.bad))
})
test_that("true positive prediction works", {
  J.str <- 100
  X.mu.0 <- matrix(1, nrow = J.str, ncol = p.tp.det)
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J.str))
})
test_that("false positive prediction works", {
  J.str <- 100
  X.omega.0 <- matrix(1, nrow = J.str, ncol = p.fp.det)
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, J.str))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, dim(y)))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, max(n.rep)))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
})

# Check plot --------------------------
test_that("plot works for random effect variances", {
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(out, 'sigma.sq.psi', density = FALSE))
})

# Multiple random intercepts on occurrence with covariate -----------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.5)
p.occ <- length(beta)
gamma <- c(1.5)
p.tp.det <- length(gamma)
delta <- c(0.5)
p.fp.det <- length(delta)
psi.RE <- list(levels = c(25, 10),
               sigma.sq.psi = c(1.3, 0.5))
mu.RE <- list()
omega.RE <- list()
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = FALSE, prop.val = 0.2)
y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.mu <- dat$X.mu
X.omega <- dat$X.omega

occ.covs <- cbind(X, X.re)
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.factor.1', 'occ.factor.2')
data.list <- list(y = y, occ.covs = occ.covs, n.val = dat$n.val, d = dat$d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10),
                   sigma.sq.psi.ig = list(a = c(0.1, 0.1), b = c(0.1, 0.1)))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0, sigma.sq.psi = c(1, 1),
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3)
occ.formula <- ~ occ.cov.1 + (1 | occ.factor.1) + (1 | occ.factor.2)
tp.det.formula <- ~ 1
fp.det.formula <- ~ 1

out <- fpPGOcc(occ.formula = occ.formula,
               tp.det.formula = tp.det.formula,
               fp.det.formula = fp.det.formula,
               data = data.list,
               inits = inits.list,
               priors = prior.list,
               tuning = tuning.list,
               n.batch = n.batch,
               batch.length = batch.length,
               n.omp.threads = 1,
               verbose = FALSE,
               n.report = n.report,
               n.burn = 400,
               n.thin = 6,
               n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class fpPGOcc", {
  expect_s3_class(out, "fpPGOcc")
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$psiRE, TRUE)
  expect_equal(out$muRE, FALSE)
  expect_equal(out$omegaRE, FALSE)
  expect_equal(dim(out$sigma.sq.psi.samples), c(out$n.post * out$n.chains, 2))
  expect_equal(ncol(out$beta.star.samples),
               length(unique(X.re[, 1])) + length(unique(X.re[, 2])))
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$re.level.names))),
               sort(unique(c(X.re))))
})

# Check RE priors and inits -----------
test_that("random effect prior and initial value errors work", {
  expect_error(fpPGOcc(occ.formula = occ.formula,
                       tp.det.formula = tp.det.formula,
                       fp.det.formula = fp.det.formula,
                       data = data.list,
                       priors = list(sigma.sq.psi.ig = list(a = rep(0.1, 3), b = 0.1)),
                       n.batch = n.batch,
                       batch.length = batch.length,
                       verbose = FALSE))
  expect_error(fpPGOcc(occ.formula = occ.formula,
                       tp.det.formula = tp.det.formula,
                       fp.det.formula = fp.det.formula,
                       data = data.list,
                       inits = list(sigma.sq.psi = c(1, 1, 1)),
                       n.batch = n.batch,
                       batch.length = batch.length,
                       verbose = FALSE))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- fpPGOcc(occ.formula = occ.formula,
                 tp.det.formula = tp.det.formula,
                 fp.det.formula = fp.det.formula,
                 data = data.list,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 n.omp.threads = 1,
                 verbose = FALSE)
  expect_s3_class(out, "fpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(fpPGOcc(occ.formula = occ.formula,
                        tp.det.formula = tp.det.formula,
                        fp.det.formula = fp.det.formula,
                        data = data.list,
                        n.batch = 4,
                        batch.length = 25,
                        n.omp.threads = 1,
                        verbose = TRUE,
                        n.report = n.report,
                        n.burn = 1,
                        n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for fpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
})

# Check predictions -------------------
test_that("predict works for fpPGOcc", {
  pred.out <- predict(out, occ.covs)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
  pred.out <- predict(out, occ.covs[, 1:2], ignore.RE = TRUE)
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("true positive prediction works", {
  J.str <- 100
  X.mu.0 <- matrix(1, nrow = J.str, ncol = p.tp.det)
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J.str))
})
test_that("false positive prediction works", {
  J.str <- 100
  X.omega.0 <- matrix(1, nrow = J.str, ncol = p.fp.det)
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, J.str))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, dim(y)))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, max(n.rep)))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
})

# Random intercept on true positive detection -----------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
gamma <- c(1.5, -0.4)
p.tp.det <- length(gamma)
delta <- c(0.5)
p.fp.det <- length(delta)
psi.RE <- list()
mu.RE <- list(levels = c(12),
              sigma.sq.mu = c(0.75))
omega.RE <- list()
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = FALSE, prop.val = 0.2)
y <- dat$y
X <- dat$X
X.mu <- dat$X.mu
X.mu.re <- dat$X.mu.re
X.omega <- dat$X.omega

data.list <- list(y = y,
                  tp.det.covs = list(tp.cov.1 = X.mu[, , 2],
                                     tp.factor.1 = X.mu.re[, , 1]),
                  n.val = dat$n.val, d = dat$d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10),
                   sigma.sq.mu.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0, sigma.sq.mu = 0.5,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3, gamma.star = 0.3)
occ.formula <- ~ 1
tp.det.formula <- ~ tp.cov.1 + (1 | tp.factor.1)
fp.det.formula <- ~ 1

out <- fpPGOcc(occ.formula = occ.formula,
               tp.det.formula = tp.det.formula,
               fp.det.formula = fp.det.formula,
               data = data.list,
               inits = inits.list,
               priors = prior.list,
               tuning = tuning.list,
               n.batch = n.batch,
               batch.length = batch.length,
               n.omp.threads = 1,
               verbose = FALSE,
               n.report = n.report,
               n.burn = 400,
               n.thin = 6,
               n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class fpPGOcc", {
  expect_s3_class(out, "fpPGOcc")
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$psiRE, FALSE)
  expect_equal(out$muRE, TRUE)
  expect_equal(out$omegaRE, FALSE)
  expect_equal(dim(out$sigma.sq.mu.samples), c(out$n.post * out$n.chains, 1))
  expect_equal(ncol(out$gamma.star.samples), length(unique(c(na.omit(c(X.mu.re))))))
  expect_equal(length(out$rhat$sigma.sq.mu), 1)
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$mu.re.level.names))),
               sort(unique(c(X.mu.re))))
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$tp.det.covs$tp.factor.1 <- factor(data.list$tp.det.covs$tp.factor.1)
  expect_error(out <- fpPGOcc(occ.formula = occ.formula,
                              tp.det.formula = tp.det.formula,
                              fp.det.formula = fp.det.formula,
                              data = data.list,
                              n.batch = n.batch,
                              batch.length = batch.length,
                              n.omp.threads = 1,
                              verbose = FALSE))
  data.list$tp.det.covs$tp.factor.1 <- as.character(factor(data.list$tp.det.covs$tp.factor.1))
  expect_error(out <- fpPGOcc(occ.formula = occ.formula,
                              tp.det.formula = tp.det.formula,
                              fp.det.formula = fp.det.formula,
                              data = data.list,
                              n.batch = n.batch,
                              batch.length = batch.length,
                              n.omp.threads = 1,
                              verbose = FALSE))
})

# Check RE tuning ---------------------
test_that("gamma.star tuning is required with true positive random effects", {
  expect_error(fpPGOcc(occ.formula = occ.formula,
                       tp.det.formula = tp.det.formula,
                       fp.det.formula = fp.det.formula,
                       data = data.list,
                       tuning = list(gamma = 0.3, delta = 0.3),
                       n.batch = n.batch,
                       batch.length = batch.length,
                       verbose = FALSE),
               "gamma.star must be specified")
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- fpPGOcc(occ.formula = occ.formula,
                 tp.det.formula = tp.det.formula,
                 fp.det.formula = fp.det.formula,
                 data = data.list,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 n.omp.threads = 1,
                 verbose = FALSE)
  expect_s3_class(out, "fpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  tmp.data$n.val[1, 1] <- NA
  tmp.data$d[1, 1] <- NA
  out <- fpPGOcc(occ.formula = occ.formula,
                 tp.det.formula = tp.det.formula,
                 fp.det.formula = fp.det.formula,
                 data = tmp.data,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 n.omp.threads = 1,
                 verbose = FALSE)
  expect_s3_class(out, "fpPGOcc")
  expect_equal(nrow(out$X.mu.re), sum(!is.na(tmp.data$y)))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(fpPGOcc(occ.formula = occ.formula,
                        tp.det.formula = tp.det.formula,
                        fp.det.formula = fp.det.formula,
                        data = data.list,
                        n.batch = 4,
                        batch.length = 25,
                        n.omp.threads = 1,
                        verbose = TRUE,
                        n.report = n.report,
                        n.burn = 1,
                        n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for fpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
})

# Check predictions -------------------
test_that("predict works for fpPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("true positive prediction works", {
  X.mu.0 <- cbind(dat$X.mu[, 1, ], dat$X.mu.re[, 1, 1])
  colnames(X.mu.0) <- c('(Intercept)', 'tp.cov.1', 'tp.factor.1')
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
  # Ignoring the random effects
  pred.out <- predict(out, dat$X.mu[, 1, ], type = 'true positive', ignore.RE = TRUE)
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
  # New random effect level
  X.mu.0[, 'tp.factor.1'] <- max(X.mu.re, na.rm = TRUE) + 1
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("false positive prediction works", {
  J.str <- 100
  X.omega.0 <- matrix(1, nrow = J.str, ncol = p.fp.det)
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, J.str))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, dim(y)))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, max(n.rep)))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
})

# Check plot --------------------------
test_that("plot works for random effect variances", {
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(out, 'sigma.sq.mu', density = FALSE))
})

# Random intercept on false positive detection ----------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3)
p.occ <- length(beta)
gamma <- c(1.5)
p.tp.det <- length(gamma)
delta <- c(0.5, 0.6)
p.fp.det <- length(delta)
psi.RE <- list()
mu.RE <- list()
omega.RE <- list(levels = c(10),
                 sigma.sq.omega = c(0.5))
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = FALSE, prop.val = 0.2)
y <- dat$y
X <- dat$X
X.mu <- dat$X.mu
X.omega <- dat$X.omega
X.omega.re <- dat$X.omega.re

data.list <- list(y = y,
                  fp.det.covs = list(fp.cov.1 = X.omega[, , 2],
                                     fp.factor.1 = X.omega.re[, , 1]),
                  n.val = dat$n.val, d = dat$d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10),
                   sigma.sq.omega.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0, sigma.sq.omega = 0.5,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3, delta.star = 0.3)
occ.formula <- ~ 1
tp.det.formula <- ~ 1
fp.det.formula <- ~ fp.cov.1 + (1 | fp.factor.1)

out <- fpPGOcc(occ.formula = occ.formula,
               tp.det.formula = tp.det.formula,
               fp.det.formula = fp.det.formula,
               data = data.list,
               inits = inits.list,
               priors = prior.list,
               tuning = tuning.list,
               n.batch = n.batch,
               batch.length = batch.length,
               n.omp.threads = 1,
               verbose = FALSE,
               n.report = n.report,
               n.burn = 400,
               n.thin = 6,
               n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class fpPGOcc", {
  expect_s3_class(out, "fpPGOcc")
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$psiRE, FALSE)
  expect_equal(out$muRE, FALSE)
  expect_equal(out$omegaRE, TRUE)
  expect_equal(dim(out$sigma.sq.omega.samples), c(out$n.post * out$n.chains, 1))
  expect_equal(ncol(out$delta.star.samples), length(unique(c(na.omit(c(X.omega.re))))))
  expect_equal(length(out$rhat$sigma.sq.omega), 1)
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$omega.re.level.names))),
               sort(unique(c(X.omega.re))))
})

# Check RE error ----------------------
test_that("random effect gives error when non-numeric", {
  data.list$fp.det.covs$fp.factor.1 <- factor(data.list$fp.det.covs$fp.factor.1)
  expect_error(out <- fpPGOcc(occ.formula = occ.formula,
                              tp.det.formula = tp.det.formula,
                              fp.det.formula = fp.det.formula,
                              data = data.list,
                              n.batch = n.batch,
                              batch.length = batch.length,
                              n.omp.threads = 1,
                              verbose = FALSE))
  data.list$fp.det.covs$fp.factor.1 <- as.character(factor(data.list$fp.det.covs$fp.factor.1))
  expect_error(out <- fpPGOcc(occ.formula = occ.formula,
                              tp.det.formula = tp.det.formula,
                              fp.det.formula = fp.det.formula,
                              data = data.list,
                              n.batch = n.batch,
                              batch.length = batch.length,
                              n.omp.threads = 1,
                              verbose = FALSE))
})

# Check RE tuning ---------------------
test_that("delta.star tuning is required with false positive random effects", {
  expect_error(fpPGOcc(occ.formula = occ.formula,
                       tp.det.formula = tp.det.formula,
                       fp.det.formula = fp.det.formula,
                       data = data.list,
                       tuning = list(gamma = 0.3, delta = 0.3),
                       n.batch = n.batch,
                       batch.length = batch.length,
                       verbose = FALSE),
               "delta.star must be specified")
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- fpPGOcc(occ.formula = occ.formula,
                 tp.det.formula = tp.det.formula,
                 fp.det.formula = fp.det.formula,
                 data = data.list,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 n.omp.threads = 1,
                 verbose = FALSE)
  expect_s3_class(out, "fpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(fpPGOcc(occ.formula = occ.formula,
                        tp.det.formula = tp.det.formula,
                        fp.det.formula = fp.det.formula,
                        data = data.list,
                        n.batch = 4,
                        batch.length = 25,
                        n.omp.threads = 1,
                        verbose = TRUE,
                        n.report = n.report,
                        n.burn = 1,
                        n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for fpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
})

# Check predictions -------------------
test_that("predict works for fpPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("true positive prediction works", {
  J.str <- 100
  X.mu.0 <- matrix(1, nrow = J.str, ncol = p.tp.det)
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J.str))
})
test_that("false positive prediction works", {
  X.omega.0 <- cbind(dat$X.omega[, 1, ], dat$X.omega.re[, 1, 1])
  colnames(X.omega.0) <- c('(Intercept)', 'fp.cov.1', 'fp.factor.1')
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, J))
  pred.out <- predict(out, dat$X.omega[, 1, ], type = 'false positive', ignore.RE = TRUE)
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, dim(y)))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, max(n.rep)))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
})

# Check plot --------------------------
test_that("plot works for random effect variances", {
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(out, 'sigma.sq.omega', density = FALSE))
})

# Random intercepts on all with covariates --------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
beta <- c(0.3, 0.5)
p.occ <- length(beta)
gamma <- c(1.5, -0.4)
p.tp.det <- length(gamma)
delta <- c(0.5, 0.6, -0.3)
p.fp.det <- length(delta)
psi.RE <- list(levels = c(20, 15),
               sigma.sq.psi = c(1.3, 0.2))
mu.RE <- list(levels = c(12),
              sigma.sq.mu = c(0.75))
omega.RE <- list(levels = c(10),
                 sigma.sq.omega = c(0.5))
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = FALSE, prop.val = 0.2)
y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.mu <- dat$X.mu
X.mu.re <- dat$X.mu.re
X.omega <- dat$X.omega
X.omega.re <- dat$X.omega.re

occ.covs <- cbind(X, X.re)
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.factor.1', 'occ.factor.2')
data.list <- list(y = y, occ.covs = occ.covs,
                  tp.det.covs = list(tp.cov.1 = X.mu[, , 2],
                                     tp.factor.1 = X.mu.re[, , 1]),
                  fp.det.covs = list(fp.cov.1 = X.omega[, , 2],
                                     fp.cov.2 = X.omega[, , 3],
                                     fp.factor.1 = X.omega.re[, , 1]),
                  n.val = dat$n.val, d = dat$d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10),
                   sigma.sq.psi.ig = list(a = 0.1, b = 0.1),
                   sigma.sq.mu.ig = list(a = 0.1, b = 0.1),
                   sigma.sq.omega.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0, sigma.sq.psi = 1,
                   sigma.sq.mu = 0.5, sigma.sq.omega = 0.5,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3, gamma.star = 0.3, delta.star = 0.3)
occ.formula <- ~ occ.cov.1 + (1 | occ.factor.1) + (1 | occ.factor.2)
tp.det.formula <- ~ tp.cov.1 + (1 | tp.factor.1)
fp.det.formula <- ~ fp.cov.1 + fp.cov.2 + (1 | fp.factor.1)

out <- fpPGOcc(occ.formula = occ.formula,
               tp.det.formula = tp.det.formula,
               fp.det.formula = fp.det.formula,
               data = data.list,
               inits = inits.list,
               priors = prior.list,
               tuning = tuning.list,
               n.batch = n.batch,
               batch.length = batch.length,
               n.omp.threads = 1,
               verbose = FALSE,
               n.report = n.report,
               n.burn = 400,
               n.thin = 6,
               n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class fpPGOcc", {
  expect_s3_class(out, "fpPGOcc")
})

# Check random effects ----------------
test_that("random effects are correct", {
  expect_equal(out$psiRE, TRUE)
  expect_equal(out$muRE, TRUE)
  expect_equal(out$omegaRE, TRUE)
  n.post.samples <- out$n.post * out$n.chains
  expect_equal(dim(out$sigma.sq.psi.samples), c(n.post.samples, 2))
  expect_equal(dim(out$sigma.sq.mu.samples), c(n.post.samples, 1))
  expect_equal(dim(out$sigma.sq.omega.samples), c(n.post.samples, 1))
})

# Check RE levels ---------------------
test_that("random effect levels are correct", {
  expect_equal(sort(unique(unlist(out$re.level.names))),
               sort(unique(c(X.re))))
  expect_equal(sort(unique(unlist(out$mu.re.level.names))),
               sort(unique(c(X.mu.re))))
  expect_equal(sort(unique(unlist(out$omega.re.level.names))),
               sort(unique(c(X.omega.re))))
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- fpPGOcc(occ.formula = occ.formula,
                 tp.det.formula = tp.det.formula,
                 fp.det.formula = fp.det.formula,
                 data = data.list,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 n.omp.threads = 1,
                 verbose = FALSE)
  expect_s3_class(out, "fpPGOcc")
})

# Check fixed initial values ----------
test_that("fixed initial values across chains work", {
  tmp.inits <- inits.list
  tmp.inits$fix <- TRUE
  out <- fpPGOcc(occ.formula = occ.formula,
                 tp.det.formula = tp.det.formula,
                 fp.det.formula = fp.det.formula,
                 data = data.list,
                 inits = tmp.inits,
                 tuning = tuning.list,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 n.chains = 2,
                 n.omp.threads = 1,
                 verbose = FALSE)
  expect_s3_class(out, "fpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check missing values ----------------
test_that("missing value error handling works", {
  tmp.data <- data.list
  tmp.data$y[1, 1] <- NA
  tmp.data$n.val[1, 1] <- NA
  tmp.data$d[1, 1] <- NA
  out <- fpPGOcc(occ.formula = occ.formula,
                 tp.det.formula = tp.det.formula,
                 fp.det.formula = fp.det.formula,
                 data = tmp.data,
                 tuning = tuning.list,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 n.omp.threads = 1,
                 verbose = FALSE)
  expect_s3_class(out, "fpPGOcc")
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(fpPGOcc(occ.formula = occ.formula,
                        tp.det.formula = tp.det.formula,
                        fp.det.formula = fp.det.formula,
                        data = data.list,
                        n.batch = 4,
                        batch.length = 25,
                        n.omp.threads = 1,
                        verbose = TRUE,
                        n.report = n.report,
                        n.burn = 1,
                        n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for fpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, dim(y)))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, dim(y)))
})

# Check predictions -------------------
test_that("predict works for fpPGOcc", {
  pred.out <- predict(out, occ.covs)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("true positive prediction works", {
  X.mu.0 <- cbind(dat$X.mu[, 1, ], dat$X.mu.re[, 1, 1])
  colnames(X.mu.0) <- c('(Intercept)', 'tp.cov.1', 'tp.factor.1')
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("false positive prediction works", {
  X.omega.0 <- cbind(dat$X.omega[, 1, ], dat$X.omega.re[, 1, 1])
  colnames(X.omega.0) <- c('(Intercept)', 'fp.cov.1', 'fp.cov.2', 'fp.factor.1')
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, dim(y)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, dim(y)))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, max(n.rep)))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, max(n.rep)))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
})

# Dimension of y is not equal to max(n.rep) -------------------------------
J.x <- 15
J.y <- 15
J <- J.x * J.y
n.rep <- sample(2:4, J, replace = TRUE)
n.rep.max <- 7
beta <- c(0.3, 0.5)
p.occ <- length(beta)
gamma <- c(1.5, -0.4)
p.tp.det <- length(gamma)
delta <- c(0.5, 0.6)
p.fp.det <- length(delta)
psi.RE <- list()
mu.RE <- list()
omega.RE <- list()
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, n.rep.max = n.rep.max,
                beta = beta, gamma = gamma, delta = delta, psi.RE = psi.RE,
                mu.RE = mu.RE, omega.RE = omega.RE, sp = FALSE, prop.val = 0.2)
y <- dat$y
X <- dat$X
X.mu <- dat$X.mu
X.omega <- dat$X.omega

colnames(X) <- c('int', 'occ.cov.1')
data.list <- list(y = y, occ.covs = X,
                  tp.det.covs = list(tp.cov.1 = X.mu[, , 2]),
                  fp.det.covs = list(fp.cov.1 = X.omega[, , 2]),
                  n.val = dat$n.val, d = dat$d)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   gamma.normal = list(mean = 0, var = 10),
                   delta.normal = list(mean = 0, var = 10))
# Starting values
inits.list <- list(beta = 0, gamma = 0, delta = 0,
                   z = as.numeric(apply(y, 1, max, na.rm = TRUE) > 0))
# Tuning
tuning.list <- list(gamma = 0.3, delta = 0.3)
occ.formula <- ~ occ.cov.1
tp.det.formula <- ~ tp.cov.1
fp.det.formula <- ~ fp.cov.1

out <- fpPGOcc(occ.formula = occ.formula,
               tp.det.formula = tp.det.formula,
               fp.det.formula = fp.det.formula,
               data = data.list,
               inits = inits.list,
               priors = prior.list,
               tuning = tuning.list,
               n.batch = n.batch,
               batch.length = batch.length,
               n.omp.threads = 1,
               verbose = FALSE,
               n.report = n.report,
               n.burn = 400,
               n.thin = 6,
               n.chains = 2)

# Test to make sure it worked ---------
test_that("out is of class fpPGOcc", {
  expect_s3_class(out, "fpPGOcc")
})

# Check random effects ----------------
test_that("random effects are empty", {
  expect_equal(out$psiRE, FALSE)
  expect_equal(out$muRE, FALSE)
  expect_equal(out$omegaRE, FALSE)
})

# Check output data output is correct -
test_that("out$y == y", {
  expect_equal(out$y, dat$y)
  expect_equal(ncol(out$y), n.rep.max)
})

# Check default priors ----------------
test_that("default priors, inits, burn, thin, tuning work", {
  out <- fpPGOcc(occ.formula = occ.formula,
                 tp.det.formula = tp.det.formula,
                 fp.det.formula = fp.det.formula,
                 data = data.list,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 n.omp.threads = 1,
                 verbose = FALSE)
  expect_s3_class(out, "fpPGOcc")
})

# Check summary -----------------------
test_that("summary works", {
  expect_output(summary(out))
})

# Check verbose -----------------------
test_that("verbose prints to the screen", {
  expect_output(fpPGOcc(occ.formula = occ.formula,
                        tp.det.formula = tp.det.formula,
                        fp.det.formula = fp.det.formula,
                        data = data.list,
                        n.batch = 4,
                        batch.length = 25,
                        n.omp.threads = 1,
                        verbose = TRUE,
                        n.report = n.report,
                        n.burn = 1,
                        n.thin = 1))
})

# Check waicOcc -----------------------
test_that("waicOCC works for fpPGOcc", {
  # as.vector gets rid of names
  waic.out <- as.vector(waicOcc(out))
  expect_equal(length(waic.out), 3)
  expect_equal(waic.out[3], -2 * (waic.out[1] - waic.out[2]))
})

# Check fitted ------------------------
test_that("fitted works for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  fitted.out <- fitted(out)
  expect_equal(length(fitted.out), 4)
  expect_equal(dim(fitted.out$y.rep.samples), c(n.post.samples, J, n.rep.max))
  expect_equal(dim(fitted.out$d.rep.samples), c(n.post.samples, J, n.rep.max))
})

# Check predictions -------------------
test_that("predict works for fpPGOcc", {
  X.0 <- dat$X
  pred.out <- predict(out, X.0)
  expect_type(pred.out, "list")
  expect_equal(dim(pred.out$psi.0.samples), c(out$n.post * out$n.chains, J))
  expect_equal(dim(pred.out$z.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("true positive prediction works", {
  X.mu.0 <- dat$X.mu[, 1, ]
  pred.out <- predict(out, X.mu.0, type = 'true positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$mu.0.samples), c(out$n.post * out$n.chains, J))
})
test_that("false positive prediction works", {
  X.omega.0 <- dat$X.omega[, 1, ]
  pred.out <- predict(out, X.omega.0, type = 'false positive')
  expect_type(pred.out, 'list')
  expect_equal(dim(pred.out$omega.0.samples), c(out$n.post * out$n.chains, J))
})

# Check PPCs --------------------------
test_that("posterior predictive checks work for fpPGOcc", {
  n.post.samples <- out$n.post * out$n.chains
  # Replicate groups are the columns of y with any observed surveys. simFPOcc
  # spreads surveys across all n.rep.max columns, so this is not max(n.rep).
  n.rep.obs <- sum(colSums(!is.na(y)) > 0)
  for (fit.stat in c('chi-square', 'freeman-tukey')) {
    ppc.out <- ppcOcc(out, fit.stat, 0)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J, n.rep.max))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J, n.rep.max))

    ppc.out <- ppcOcc(out, fit.stat, 1)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, J))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, J))

    ppc.out <- ppcOcc(out, fit.stat, 2)
    expect_type(ppc.out, "list")
    expect_equal(length(ppc.out$fit.y), n.post.samples)
    expect_equal(length(ppc.out$fit.d), n.post.samples)
    expect_equal(dim(ppc.out$fit.y.group.quants), c(5, n.rep.obs))
    expect_equal(dim(ppc.out$fit.d.group.quants), c(5, n.rep.obs))
  }
})

# Test residuals ----------------------
test_that("residuals works", {
  out.resids <- residuals(out, n.post.samples = 10)
  expect_type(out.resids, 'list')
  expect_equal(length(out.resids), 3)
  expect_equal(dim(out.resids$det.resids), c(10, J, n.rep.max))
})
