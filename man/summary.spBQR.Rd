% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.spBQR.R
\name{summary.spBQR}
\alias{summary.spBQR}
\title{Summary methods for Bayesian quantile regression models.}
\usage{
\method{summary}{spBQR}(object, burnin = 1000, ci = 0.95, ...)
}
\arguments{
\item{object}{This is an object of class "spBQR", produced by a call to the
spBQR function.}

\item{burnin}{Initial part of the chain, which is to be discarded. Default
value is 1000.}

\item{ci}{Credible interval coefficient. Default value is 0.95.}

\item{...}{other summary parameters (currently not used)}
}
\value{
A data frame with summary information about the quantile regression
 parameters.
}
\description{
Returns a summary data.frame for a Bayesian quantile regression fit for more
 than one quantile.
}
\examples{
set.seed(1)
}
