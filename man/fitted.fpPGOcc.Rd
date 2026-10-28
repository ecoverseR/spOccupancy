\name{fitted.fpPGOcc}

\alias{fitted.fpPGOcc}

\title{Extract Model Fitted Values for fpPGOcc Object}

\description{
  Method for extracting model fitted values and detection probabilities from a fitted single-species false positive occupancy (\code{fpPGOcc}) model. 
}

\usage{
\method{fitted}{fpPGOcc}(object, \dots)
}

\arguments{
  \item{object}{object of class \code{fpPGOcc}.}
  \item{\dots}{currently no additional arguments}
}

\details{
  A method to the generic \code{\link[stats]{fitted}} function to extract fitted values and true and false positive rates for fitted model objects of class \code{fpPGOcc}. 
}

\value{
  A list comprised of: 

  \item{y.rep.samples}{A three-dimensional numeric array of fitted values of the observation data for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, sites, and replicates.}
  
  \item{d.rep.samples}{A three-dimensional numeric array of fitted values of the human validation data for use in Goodness of Fit assessments. Array dimensions correspond to MCMC samples, sites, and replicates.}

  \item{mu.samples}{A three-dimensional numeric array of true positive rate values. Array dimensions correspond to MCMC samples, sites, and replicates.}
  
  \item{omega.samples}{A three-dimensional numeric array of false positive rate values. Array dimensions correspond to MCMC samples, sites, and replicates.}
}

\keyword{model}

