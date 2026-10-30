\name{predict.fpPGOcc}
\alias{predict.fpPGOcc}
\title{Function for prediction at new locations for single-species false positive occupancy models}

\description{
  The function \code{predict} collects posterior predictive samples for a set of new locations given an object of class `fpPGOcc`. Prediction is possible for both the latent occupancy state as well as the true positive detection rate and false positive detection rate. 
}

\usage{
\method{predict}{fpPGOcc}(object, X.0, ignore.RE = FALSE, type = 'occupancy', ...)
}

\arguments{

  \item{object}{an object of class fpPGOcc}

  \item{X.0}{the design matrix of covariates at the prediction locations. This should include a column of 1s for the intercept if an intercept is included in the model. If random effects are included in the portion of the model that is being predicted, the levels of the random effects at the new locations should be included as a column in the design matrix. The ordering of the levels should match the ordering used to fit the data in \code{fpPGOcc}. Columns should correspond to the order of how covariates were specified in the corresponding formula argument of \code{fpPGOcc}. Column names of the random effects must match the name of the random effects, if specified in the corresponding formula argument of \code{fpPGOcc}.}

  \item{ignore.RE}{logical value that specifies whether or not to remove random effects from the subsequent predictions. If \code{TRUE}, random effects will be included. If \code{FALSE}, random effects will be set to 0 and predictions will only be generated from the fixed effects.}

  \item{type}{a quoted keyword indicating what type of prediction to produce. Valid keywords are 'occupancy' to predict latent occupancy probability and latent occupancy values (this is the default), 'true positive' to predict the true positive detection rates given new values of covariates, and 'false positive' to predict the false positive rates given new values of covariates.}

  \item{...}{currently no additional arguments}
}

\note{
  When \code{ignore.RE = FALSE}, both sampled levels and non-sampled levels of random effects are supported for prediction. For sampled levels, the posterior distribution for the random intercept corresponding to that level of the random effect will be used in the prediction. For non-sampled levels, random values are drawn from a normal distribution using the posterior samples of the random effect variance, which results in fully propagated uncertainty in predictions with models that incorporate random effects. 
}

\author{
  Jeffrey W. Doser \email{jwdoser@ncsu.edu}
}

\value{
  A list object of class \code{predict.fpPGOcc}. When \code{type = 'occupancy'}, the list consists of: 

  \item{psi.0.samples}{a \code{coda} object of posterior predictive samples for the
    latent occupancy probability values.}

  \item{z.0.samples}{a \code{coda} object of posterior predictive samples for the
    latent occupancy values.}

  When \code{type = 'true positive'}, the list consists of: 

  \item{mu.0.samples}{a \code{coda} object of posterior predictive samples for the 
    true positive rate values.}

  When \code{type = 'false positive'}, the list consists of: 

  \item{omega.0.samples}{a \code{coda} object of posterior predictive samples for the 
    false positive rate values.}

  The return object will include additional objects used for standard 
  extractor functions. 
}

\examples{

set.seed(8372)
J.x <- 10
J.y <- 10
J <- J.x * J.y
n.rep <- sample(3, J, replace = TRUE)
n.rep.max <- max(n.rep)
beta <- c(0, 1.3)
p.occ <- length(beta)
# Covariates on true positive detection rate
gamma <- c(1, -0.5)
p.tp <- length(gamma)
# False positive rate
delta <- c(2, 0.3, -0.8)
p.fp <- length(delta)

# Proportion of samples manually validated
prop.val <- 0.05

# Simulate the data
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, sp = FALSE, prop.val = prop.val, family = "Poisson")

# Split into fitting and prediction data set
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
y <- dat$y[-pred.indx, , drop = FALSE]
# Occupancy covariates
X <- dat$X[-pred.indx, , drop = FALSE]
# Prediction covariates
X.0 <- dat$X[pred.indx, , drop = FALSE]
# True positive detection covariates
X.omega <- dat$X.omega[-pred.indx, , , drop = FALSE]
# False positive detection covariates
X.mu <- dat$X.mu[-pred.indx, , , drop = FALSE]
# Human validation data
n.val <- dat$n.val[-pred.indx, , drop = FALSE]
d <- dat$d[-pred.indx, , drop = FALSE]

# Package all data into a list
occ.covs <- cbind(X)
colnames(occ.covs) <- c('int', 'occ.cov.1')
# True positive covariates
tp.det.covs <- list(tp.cov.1 = X.mu[, , 2])
# False positive covariates
fp.det.covs <- list(fp.cov.1 = X.omega[, , 2], fp.cov.2 = X.omega[, , 3])
data.list <- list(y = y,
                  occ.covs = occ.covs,
                  tp.det.covs = tp.det.covs,
                  fp.det.covs = fp.det.covs,
                  n.val = n.val,
                  d = d)

# Priors
prior.list <- list(beta.normal = list(0, 2.72),
                   gamma.normal = list(mean = rep(0, p.tp),
                                       var = rep(10, p.tp)),
                   delta.normal = list(mean = rep(0, p.fp),
                                       var = rep(10, p.fp)))

# Initial values
inits.list <- list(gamma = 0, beta = 0, delta = 0)

# Tuning
tuning.list <- list(delta = 0.5, gamma = 0.5)

# Run a short test run of the model
out <- fpPGOcc(occ.formula = ~ occ.cov.1,
               tp.det.formula = ~ tp.cov.1,
               fp.det.formula = ~ fp.cov.1 + fp.cov.2,
               data = data.list,
               inits = inits.list,
               priors = prior.list,
               n.batch = 5,
               tuning = tuning.list,
               batch.length = 25,
               accept.rate = 0.43,
               n.omp.threads = 1,
               verbose = TRUE,
               n.report = 1,
               n.chains = 1,
               n.burn = 0,
               n.thin = 1)

# Predict at new locations ------------------------------------------------
colnames(X.0) <- c('intercept', 'occ.cov.1')
out.pred <- predict(out, X.0)
psi.0.quants <- apply(out.pred$psi.0.samples, 2, quantile, c(0.025, 0.5, 0.975))
plot(dat$psi[pred.indx], psi.0.quants[2, ], pch = 19, xlab = 'True', 
     ylab = 'Fitted', ylim = c(min(psi.0.quants), max(psi.0.quants)))
segments(dat$psi[pred.indx], psi.0.quants[1, ], dat$psi[pred.indx], psi.0.quants[3, ])
lines(dat$psi[pred.indx], dat$psi[pred.indx])
}
