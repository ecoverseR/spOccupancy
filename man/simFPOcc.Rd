\name{simFPOcc}
\alias{simFPOcc}
\title{Simulate Single-Species Data from Automated Algorithms with False Positive Classifications}

\description{
  The function \code{simFPOcc} simulates single-species data assumed to arise from some automated algorithm that processed remote sensing data (e.g., acoustic recording units, camera traps). These data are assumed to contain both false positive and false negative observations. This function is designed to be used for simulation studies, power assessments, or function testing. Data can be optionally simulated with a spatial Gaussian Process in the occurrence portion of the model. Non-spatial random intercepts can also be included in the observation models (i.e., true positive detection rate, false positive rate) and the occurrence model. 
}

\usage{
simFPOcc(J.x, J.y, n.rep, n.rep.max, beta, gamma, delta, psi.RE = list(), 
       mu.RE = list(), omega.RE = omega.RE, sp = FALSE, svc.cols = 1, cov.model, 
       sigma.sq, phi, nu, x.positive = FALSE, grid, prop.val, kappa, 
       family = "Poisson", ...)
}

\arguments{
\item{J.x}{a single numeric value indicating the number of sites to simulate detection-nondetection data along the horizontal axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{J.y}{a single numeric value indicating the number of sites to simulate detection-nondetection data along the vertical axis. Total number of sites with simulated data is \eqn{J.x \times J.y}{J.x x J.y}. }

\item{n.rep}{a numeric vector of length \eqn{J = J.x \times J.y}{J = J.x x J.y} indicating the number of repeat visits at each of the \eqn{J}{J} sites.}

\item{n.rep.max}{a single numeric value indicating the maximum number of replicate surveys. This is an optional argument, with its default value set to \code{max(n.rep)}. This can be used to generate data sets with different types of missingness (e.g., simulate data across 20 days (replicate surveys) but sites are only sampled a maximum of ten times each).}

\item{beta}{a numeric vector containing the intercept and regression coefficient parameters for the occupancy portion of the single-species occupancy model.}

\item{gamma}{a numeric vector containing the intercept and regression coefficient parameters for the true positive rate observation model of the single-species occupancy model.}

\item{delta}{a numeric vector containing the intercept and regression coefficient parameters for the false positive rate observation model of the single-species occupancy model.}

\item{psi.RE}{a list used to specify the non-spatial random intercepts included in the occupancy portion of the model. The list must have two tags: \code{levels} and \code{sigma.sq.psi}. \code{levels} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the number of levels there are in each intercept. \code{sigma.sq.psi} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the variances for each random effect. If not specified, no random effects are included in the occupancy portion of the model.}

\item{mu.RE}{a list used to specify the non-spatial random intercepts included in the true positive rate observation model. The list must have two tags: \code{levels} and \code{sigma.sq.mu}. \code{levels} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the number of levels there are in each intercept. \code{sigma.sq.mu} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the variances for each random effect. If not specified, no random effects are included in the true positive rate observation submodel.}

\item{omega.RE}{a list used to specify the non-spatial random intercepts included in the false positive rate observation model. The list must have two tags: \code{levels} and \code{sigma.sq.omega}. \code{levels} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the number of levels there are in each intercept. \code{sigma.sq.omega} is a vector of length equal to the number of distinct random intercepts to include in the model and contains the variances for each random effect. If not specified, no random effects are included in the false positive rate observation submodel.}

\item{sp}{a logical value indicating whether to simulate a spatially-explicit occupancy model with a Gaussian process. By default set to \code{FALSE}.}

\item{svc.cols}{a vector indicating the variables whose effects will be
  estimated as spatially-varying coefficients. \code{svc.cols} is an
  integer vector with values indicating the order of covariates specified
  in the model formula (with 1 being the intercept if specified).}

\item{cov.model}{a quoted keyword that specifies the covariance function used to model the spatial dependence structure among the latent occurrence values.  Supported covariance model key words are: \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and \code{"gaussian"}.}

\item{sigma.sq}{a numeric value indicating the spatial variance parameter. Ignored when \code{sp = FALSE}.}

\item{phi}{a numeric value indicating the spatial decay parameter. Ignored when \code{sp = FALSE}. }

\item{nu}{a numeric value indicating the spatial smoothness parameter. Only used when \code{sp = TRUE} and \code{cov.model = "matern"}.} 

\item{x.positive}{a logical value indicating whether the simulated covariates should be simulated as random standard normal covariates (\code{x.positive = FALSE}) or restricted to positive values using a uniform distribution with lower bound 0 and upper bound 1 (\code{x.positive = TRUE}).}

\item{grid}{an atomic vector used to specify the grid across which to simulate the latent spatial processes. This argument is used to simulate the underlying spatial processes at a different resolution than the coordinates (e.g., if coordinates are distributed across a grid).}

\item{prop.val}{a single numeric value ranging between 0 and 1 that indicates the proportion of the total signals that are manually vetted by a human observer.}

\item{kappa}{a single numeric value containing the dispersion parameter for the Poisson model of the detected signals. Only relevant when \code{family = 'NB'}.}

\item{family}{the distribution to use for the observation model. Currently supports \code{"Poisson"} (default) and \code{"NB"} (negative binomial).}

\item{...}{currently no additional arguments}
}

\author{
  Jeffrey W. Doser \email{jwdoser@ncsu.edu}
}

\references{
  Chambert, T., Waddle, J. H., Miller, D. A., Walls, S. C., & Nichols, J. D. (2018). 
  A new framework for analysing automated acoustic species detection data: 
  Occupancy estimation and optimization of recordings post‐processing. 
  Methods in Ecology and Evolution, 9(3), 560-570.
}


\value{
  A list comprised of: 

  \item{X}{a \eqn{J \times p.occ}{J x p.occ} numeric design matrix for the occupancy portion of the model.}
  \item{X.mu}{a three-dimensional numeric array with dimensions corresponding to sites, repeat visits, and number of true positive rate regression coefficients. This is the design matrix used for the true positive rate sub-model of the occupancy model.}
  \item{X.omega}{a three-dimensional numeric array with dimensions corresponding to sites, repeat visits, and number of false positive rate regression coefficients. This is the design matrix used for the false positive rate sub-model of the occupancy model.}
  \item{coords}{a matrix of coordiantes of each occupancy site, with the rows corresponding to a unique spatial coordinate and the columns corresponding to the easting and northing, respectively. This will be identical to coords.full unless \code{grid} is specified as an argument.}
  \item{coords.full}{a \eqn{J \times 2}{J x 2} numeric matrix of coordinates of each occupancy site. Required for spatial models. Note that this may contain repeated coordinates if setting \code{grid} as an function argument.}
  \item{w}{a matrix of the spatial random effect values for each site. The number of columns is determined by the \code{svc.cols} argument (the number of spatially-varying coefficients).}
  \item{psi}{a \eqn{J \times 1}{J x 1} matrix of the occupancy probabilities for each site.}
  \item{z}{a length \eqn{J}{J} vector of the latent occupancy states at each site.}
  \item{mu}{a \code{J x max(n.rep)} matrix of the true positive rates for each site and replicate combination. Sites with fewer than \code{max(n.rep)} replicates will contain \code{NA} values.}
  \item{omega}{a \code{J x max(n.rep)} matrix of the false positive rates for each site and replicate combination. Sites with fewer than \code{max(n.rep)} replicates will contain \code{NA} values.}
  \item{y}{a \code{J x max(n.rep)} matrix of the raw detection-nondetection data for each site and replicate combination.}
  \item{n.val}{a \code{J x max(n.rep)} matrix of the number of human-validated signals at each site and replicate combination.}
  \item{d}{a \code{J x max(n.rep)} matrix of the number of human-validated signals at each site and replicate combination that were validated to be true.}
  \item{D}{a \code{J x max(n.rep)} matrix of the total number of true positive signals at each site and replicate combination.}
  \item{Q}{a \code{J x max(n.rep)} matrix of the total number of false positive signalst at each site and replicate combination.}
  \item{X.mu.re}{a three-dimensional numeric array containing the levels of any true positive rate random effect included in the model. Only relevant when true positive rate random effects are specified in \code{mu.RE}.}
  \item{X.omega.re}{a three-dimensional numeric array containing the levels of any false positive rate random effect included in the model. Only relevant when false positive rate random effects are specified in \code{omega.RE}.}
  \item{X.re}{a numeric matrix containing the levels of any occurrence random effect included in the model. Only relevant when occurrence random effects are specified in \code{psi.RE}.}
  \item{gamma.star}{a numeric vector that contains the simulated true positive rate random effects for each given level of the random effects included in the true positive rate portion of the observation model. Only relevant when true positive rate random effects are included in the model.}
  \item{delta.star}{a numeric vector that contains the simulated false positive rate random effects for each given level of the random effects included in the false positive rate portion of the observation model. Only relevant when false positive rate random effects are included in the model.}
  \item{beta.star}{a numeric vector that contains the simulated occurrence random effects for each given level of the random effects included in the occurrence model. Only relevant when occurrence random effects are included in the model.}
}

\examples{
set.seed(9928)
J.x <- 10
J.y <- 10
n.rep <- rep(4, J.x * J.y)
beta <- c(0.5, -0.15)
gamma <- c(2, 0.4)
delta <- c(1, -0.8)
phi <- 3 / .6
sigma.sq <- 2
psi.RE <- list(levels = 10,
               sigma.sq.psi = 1.2)
mu.RE <- list(levels = 15,
              sigma.sq.mu = 0.8)
omega.RE <- list(levels = 11,
                 sigma.sq.omega = 1.2)
prop.val <- 0.05
dat <- simFPOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, gamma = gamma,
                delta = delta, psi.RE = psi.RE, mu.RE = mu.RE, omega.RE = omega.RE,
                sp = TRUE, sigma.sq = sigma.sq, phi = phi, cov.model = "exponential",
                prop.val = prop.val, family = "Poisson")
}
