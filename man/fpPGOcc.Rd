\name{fpPGOcc}
\alias{fpPGOcc}
\title{Function for Fitting Single-Species False Positive Occupancy Models}

\usage{
fpPGOcc(occ.formula, tp.det.formula, fp.det.formula, data, inits, priors, tuning, 
        n.batch, batch.length, accept.rate = 0.43, n.omp.threads = 1, verbose = TRUE, 
        n.report = 100, n.burn = round(.10 * n.batch * batch.length), 
        n.thin = 1, n.chains = 1, ...)
}

\description{
  Function for fitting single-species occupancy models using Polya-Gamma latent variables. 
  The function assumes data arise from some automated algorithm that processed remote sensing
  data (e.g., acoustic recording units, camera traps) to generate a number of observations 
  for each site x replicate combination (i.e., count data). A human observer is assumed to 
  verify a subset of the observations to help separately identify the false positive and 
  true positive detection rates. 
}

\arguments{
\item{occ.formula}{a symbolic description of the model to be fit
  for the occurrence portion of the model using R's model syntax. Only
  right-hand side of formula is specified. See example below. Random intercepts
  are allowed using lme4 syntax (Bates et al. 2015).}

\item{tp.det.formula}{a symbolic description of the model to be fit
  for the true positive rate detection portion of the model using R's model syntax. Only
  right-hand side of formula is specified. See example below. Random intercepts
  are allowed using lme4 syntax (Bates et al. 2015).}

\item{fp.det.formula}{a symbolic description of the model to be fit
  for the false positive rate detection portion of the model using R's model syntax. Only
  right-hand side of formula is specified. See example below. Random intercepts
  are allowed using lme4 syntax (Bates et al. 2015).}

\item{data}{a list containing data necessary for model fitting.
  Valid tags are \code{y}, \code{occ.covs}, \code{tp.det.covs}, \code{fp.det.covs}, 
  \code{n.val}, and \code{d}. \code{y} is a matrix or data frame with 
  first dimension equal to the number of sites (\eqn{J}{J}) and second dimension equal to the maximum number of 
  replicates at a given site. \code{y} contains the number of automated detections
  of the target species at each site x replicate combination. \code{occ.covs} is a matrix or data frame
  containing the variables used in the occurrence portion of the model, with
  \eqn{J}{J} rows for each column (variable). \code{tp.det.covs} is a list of
  variables included in the true positive rate detection portion of the model. Each list element
  is a different detection covariate, which can be site-level or 
  observational-level. Site-level covariates are specified as a vector 
  of length \eqn{J}{J} while observation-level covariates are specified 
  as a matrix or data frame with the number of rows equal to \eqn{J}{J} and 
  number of columns equal to the maximum number of replicates at a given site. 
  \code{fp.det.covs} takes the same form as \code{tp.det.covs} but contains covarites for the 
  false positive rate detection portion of the model. \code{n.val} is a two-dimensional 
  matrix with rows corresponding to sites and columns corresponding to replicates that contains
  the number of observations for each site x replicate cell that were validated by a human 
  observer. \code{d} is a two-dimensional matrix with rows corresponding to sites and columns
  corresponding to replicates that contains the number of the \code{n.val} observations that 
  were validated to be true during the manual validation process.}

\item{inits}{a list with each tag corresponding to a parameter name.
  Valid tags are \code{z}, \code{beta}, \code{gamma}, \code{sigma.sq.psi}, 
  \code{sigma.sq.mu}, and \code{sigma.sq.omega}. The value portion of each tag is the 
  parameter's initial value. \code{sigma.sq.psi}, \code{sigma.sq.mu}, and \code{sigma.sq.omega} are 
  only relevant when including random effects in the occurrence, true positive rate, and 
  false positive rate portions of the occupancy model, respectively. See \code{priors} 
  description for definition of each parameter name.
  Additionally, the tag \code{fix} can be set to \code{TRUE} 
  to fix the starting values across all chains. If \code{fix} is not specified
  (the default), starting values are varied randomly across chains.}

\item{priors}{a list with each tag corresponding to a parameter name. 
  Valid tags are \code{beta.normal}, \code{gamma.normal}, \code{delta.normal}, 
  \code{sigma.sq.psi.ig}, \code{sigma.sq.mu.ig}, and \code{sigma.sq.omega.ig}. 
  Occupancy (\code{beta}), true positive rate (\code{omega}), and false positive rate (\code{delta}) 
  regression coefficients are assumed to follow a normal distribution. 
  The hyperparameters of the normal distribution are passed as a list of 
  length two with the first and second elements corresponding to the mean 
  and variance of the normal distribution, which are each specified as vectors of 
  length equal to the number of coefficients to be estimated or of length one
  if priors are the same for all coefficients. If not specified, prior means 
  are set to 0 and prior variances set to 2.72 for the occupancy coefficients and 10 for 
  the true positive and false positive detection rate coefficients. \code{sigma.sq.psi}, \code{sigma.sq.mu}, and  
  \code{sigma.sq.omega} are the random effect variances for any occurrence, true positive rate, 
  or false positive rate random effects, respectively, and are assumed to follow an inverse
  Gamma distribution. The hyperparameters of the inverse-Gamma distribution
  are passed as a list of length two with first and second elements corresponding
  to the shape and scale parameters, respectively, which are each specified as
  vectors of length equal to the number of random intercepts or of length one 
  if priors are the same for all random effect variances.}

\item{tuning}{a list with each tag corresponding to a parameter
  name. Valid tags are \code{gamma}, \code{delta}, \code{sigma.sq.mu}, 
  \code{sigma.sq.omega}, \code{gamma.star} (the true-positive rate random intercept values), 
  and \code{delta.star} (the false-positive rate random intercept values).
  The value portion of each tag defines the initial variance of the Adaptive sampler. Note 
  that the introduction of count data through \code{y} and \code{d} eliminates the use of 
  Polya-Gamma data augmentation for the detection model parameter updates, and so all detection 
  model parameters are updated with Metropolis (or Metropolis-Hastings) updates using an 
  Adaptive sampler. See Roberts and Rosenthal (2009) for details.}

\item{n.batch}{the number of MCMC batches in each chain to run for the Adaptive MCMC
  sampler. See Roberts and Rosenthal (2009) for details.}

\item{batch.length}{the length of each MCMC batch in each chain to run for the Adaptive
  MCMC sampler. See Roberts and Rosenthal (2009) for details.}

\item{accept.rate}{target acceptance rate for Adaptive MCMC. Default is
  0.43. See Roberts and Rosenthal (2009) for details.}

\item{n.omp.threads}{a positive integer indicating the number of threads
  to use for SMP parallel processing within-chains. This will have no impact
  on model run time for non-spatial models. The package must be compiled for 
  OpenMP support. For most Intel-based machines, we recommend setting 
  \code{n.omp.threads} up to the number of hypterthreaded cores. Note, 
  \code{n.omp.threads} > 1 might not work on some systems. Currently only 
  relevant for spatial models.}

\item{verbose}{if \code{TRUE}, messages about data preparation, 
  model specification, and progress of the sampler are printed to the screen. 
  Otherwise, no messages are printed.}

\item{n.report}{the interval to report MCMC progress.}

\item{n.burn}{the number of samples out of the total \code{n.samples} to 
  discard as burn-in for each chain. By default, the first 10\% of samples is discarded.}

\item{n.thin}{the thinning interval for collection of MCMC samples. The
  thinning occurs after the \code{n.burn} samples are discarded. Default 
  value is set to 1.}

\item{n.chains}{the number of chains to run.}

\item{...}{currently no additional arguments}
}

\note{
 Some of the underlying code used for generating random numbers from the
 Polya-Gamma distribution is taken from the \pkg{pgdraw} package written
 by Daniel F. Schmidt and Enes Makalic. Their code implements Algorithm
 6 in PhD thesis of Jesse Bennett Windle (2013) \url{https://repositories.lib.utexas.edu/handle/2152/21842}.
}

\references{

  Chambert, T., Waddle, J. H., Miller, D. A., Walls, S. C., & Nichols, J. D. (2018). 
  A new framework for analysing automated acoustic species detection data: 
  Occupancy estimation and optimization of recordings post‐processing. 
  Methods in Ecology and Evolution, 9(3), 560-570.
  
  MacKenzie, D. I., J. D. Nichols, G. B. Lachman, S. Droege, 
  J. Andrew Royle, and C. A. Langtimm. 2002. Estimating Site Occupancy 
  Rates When Detection Probabilities Are Less Than One. Ecology 83: 2248-2255.

  Polson, N.G., J.G. Scott, and J. Windle. (2013) Bayesian Inference for
  Logistic Models Using Polya-Gamma Latent Variables. 
  \emph{Journal of the American Statistical Association}, 108:1339-1349.

  Bates, Douglas, Martin Maechler, Ben Bolker, Steve Walker (2015).
  Fitting Linear Mixed-Effects Models Using lme4. Journal of
  Statistical Software, 67(1), 1-48. \doi{10.18637/jss.v067.i01}.

  Roberts, G.O. and Rosenthal J.S. (2009) Examples  of adaptive MCMC. 
  \emph{Journal of Computational and Graphical Statistics}, 18(2):349-367.

}

\author{
  Jeffrey W. Doser \email{jwdoser@ncsu.edu}
}

\value{
  An object of class \code{fpPGOcc} that is a list comprised of: 

  \item{beta.samples}{a \code{coda} object of posterior samples
    for the occupancy regression coefficients.}

  \item{gamma.samples}{a \code{coda} object of posterior samples
    for the true positive rate regression coefficients.}
  
  \item{delta.samples}{a \code{coda} object of posterior samples
    for the false positive rate regression coefficients.}

  \item{z.samples}{a \code{coda} object of posterior samples 
    for the latent occupancy values}

  \item{psi.samples}{a \code{coda} object of posterior samples
    for the latent occupancy probability values}

  \item{sigma.sq.psi.samples}{a \code{coda} object of posterior samples
    for variances of random intercepts included in the occupancy portion
    of the model. Only included if random intercepts are specified in 
    \code{occ.formula}.}

  \item{sigma.sq.mu.samples}{a \code{coda} object of posterior samples
    for variances of random intercpets included in the true positive rate portion 
    of the model. Only included if random intercepts are specified in 
    \code{tp.det.formula}.}

  \item{sigma.sq.omega.samples}{a \code{coda} object of posterior samples
    for variances of random intercpets included in the false positive rate portion 
    of the model. Only included if random intercepts are specified in 
    \code{fp.det.formula}.}

  \item{beta.star.samples}{a \code{coda} object of posterior samples
    for the occurrence random effects. Only included if random intercepts
    are specified in \code{occ.formula}.}

  \item{gamma.star.samples}{a \code{coda} object of posterior samples
    for the true positive rate random effects. Only included if random intercepts 
    are specified in \code{tp.det.formula}.}
  
  \item{delta.star.samples}{a \code{coda} object of posterior samples
    for the false positive rate random effects. Only included if random intercepts 
    are specified in \code{fp.det.formula}.}

  \item{like.samples}{a \code{coda} object of posterior samples
    for the likelihood value associated with each site. Used for calculating
    WAIC.}

  \item{rhat}{a list of Gelman-Rubin diagnostic values for some of the model
    parameters.}

  \item{ESS}{a list of effective sample sizes for some of the model parameters.}

  \item{run.time}{execution time reported using \code{proc.time()}.}

  The return object will include additional objects used for 
  subsequent prediction and/or model fit evaluation. 
}

\examples{
set.seed(8372)
J.x <- 5
J.y <- 5
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

# Extract necessary data from the data list
y <- dat$y
# Occupancy covariates
X <- dat$X
# True positive rate covariates
X.mu <- dat$X.mu
# False positive rate covariates
X.omega <- dat$X.omega
# Human validation data
n.val <- dat$n.val
# True positive validated data
d <- dat$d

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
                  coords = dat$coords,
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

summary(out)

}
