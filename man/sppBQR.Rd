% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sppBQR.R
\name{sppBQR}
\alias{sppBQR}
\title{Spatial Bayesian quantile regression models using predictive processes}
\usage{
sppBQR(formula, tau = 0.5, data, itNum, thin = 1, betaValue = NULL,
  sigmaValue = 1, spCoord1, spCoord2, lambdaVec, lambda = sample(lambdaVec,
  1), shapeL = 1, rateL = 50, tuneP = 1, m,
  indexes = sample((1:dim(data)[1] - 1), size = m), alpha = 0.5,
  tuneA = 1000, priorVar = 100, refresh = 100, quiet = T,
  jitter = 1e-10, includeAlpha = TRUE, tuneV = 0.5, kMT = 5,
  discLambda = FALSE)
}
\arguments{
\item{formula}{a formula object, with the response on the left of a ~
operator, and the terms, separated by + operators, on the right.}

\item{tau}{Quantile of interest.}

\item{data}{A data.frame from which to find the variables defined in the
formula}

\item{itNum}{Number of iterations.}

\item{thin}{Thinning parameter. Default value is 1.}

\item{betaValue}{Initial values for the parameter beta for the continuous
part.}

\item{sigmaValue}{Initial value for the scale parameter.}

\item{spCoord1}{Name of the first spatial coordinate, as character.}

\item{spCoord2}{Name of the second spatial coordinate, as character.}

\item{lambdaVec}{Vector of lambdas to be used in the estimation process.}

\item{lambda}{Initial value for the parameter in the covariance matrix.}

\item{shapeL}{Shape hyperparameter value for Gamma prior for the lambda
parameter.}

\item{rateL}{Rate hyperparameter value for Gamma prior for the lambda
parameter.}

\item{tuneP}{Tuning parameter for the Metropolis-Hastings algorithm to draw
samples from the posterior distribution of kappa.}

\item{m}{Number of knots.}

\item{indexes}{Vector with indexes from the knots. These are randomly
selected given the integer m.}

\item{alpha}{Value between 0 and 1 of the pure error variance in the
covariance matrix. Default is 0.5.}

\item{tuneA}{Tuning parameter for the Metropolis_Hastings algorithm to draw
samples from the posterior distribution of alpha.}

\item{priorVar}{Value that multiplies an identity matrix in the elicition
process of the prior variance of the regression parameters.}

\item{refresh}{Interval between printing a message during the iteration
process. Default is set to 100.}

\item{quiet}{If TRUE, the default, it does not print messages to check if
the MCMC is actually updating. If FALSE, it will use the value of refresh
to print messages to control the iteration process.}

\item{jitter}{Ammount of jitter added to help in the process for inverting
the covariance matrix. Default is 1e-10.}

\item{includeAlpha}{If TRUE, the default, the model will include the alpha
parameter. If FALSE, alpha is set to zero for all draws of the chain.}

\item{tuneV}{Tuning parameter to the multiple-try Metropolis to sample for
the posterior distribution of the latent variables. Default value is 0.5.}

\item{kMT}{Integer, number of Metropolis samples in the multiple-try
Metropolis. Default value is 5.}

\item{discLambda}{use discretization for the lambda parameter}
}
\value{
A list with the chains of all parameters of interest.
}
\description{
This function estimates a spatial Bayesian quantile regression model
The Asymmetric Laplace Predictive Process (ALPP) is considered to fit this
 spatial model.
}
\examples{
set.seed(1)
}
\references{
Lum and Gelfand (2012) - Spatial Quantile Multiple Regression
 Using the Asymmetric Laplace process. Bayesian Analysis.
}
