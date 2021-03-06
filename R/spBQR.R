#' Spatial Bayesian quantile regression models
#'
#' This function estimates a spatial Bayesian quantile regression model
#' The Asymmetric Laplace Process (ALP) is considered to fit this
#'  spatial model.
#'
#' @param formula a formula object, with the response on the left of a ~
#'  operator, and the terms, separated by + operators, on the right.
#' @param tau Quantile of interest.
#' @param data A data.frame from which to find the variables defined in the
#' formula
#' @param itNum Number of iterations.
#' @param thin Thinning parameter. Default value is 1.
#' @param betaValue Initial values for the parameter beta for the continuous
#' part.
#' @param sigmaValue Initial value for the scale parameter.
#' @param spCoord1 Name of the first spatial coordinate, as character.
#' @param spCoord2 Name of the second spatial coordinate, as character.
#' @param lambda Initial value for the parameter in the covariance matrix.
#' @param tuneP Tuning parameter for the Metropolis-Hastings algorithm to draw
#'  samples from the posterior distribution of kappa.
#' @param alpha Value between 0 and 1 of the pure error variance in the
#'  covariance matrix. Default is 0.5.
#' @param tuneA Tuning parameter for the Metropolis_Hastings algorithm to draw
#'  samples from the posterior distribution of alpha.
#' @param priorVar Value that multiplies an identity matrix in the elicition
#'  process of the prior variance of the regression parameters.
#' @param refresh Interval between printing a message during the iteration
#'  process. Default is set to 100.
#' @param quiet If TRUE, the default, it does not print messages to check if
#'  the MCMC is actually updating. If FALSE, it will use the value of refresh
#'  to print messages to control the iteration process.
#' @param jitter add a small value to invert spatial correlation matrix
#' @param includeAlpha If TRUE, the default, the model will include the alpha
#'  parameter. If FALSE, alpha is set to zero for all draws of the chain.
#' @param tuneV Tuning parameter to the multiple-try Metropolis to sample for
#'  the posterior distribution of the latent variables. Default value is 0.5.
#' @param kMT Integer, number of Metropolis samples in the multiple-try
#'  Metropolis. Default value is 5.
#' @return A list with the chains of all parameters of interest.
#' @references Lum and Gelfand (2012) - Spatial Quantile Multiple Regression
#'  Using the Asymmetric Laplace process. Bayesian Analysis.
#' @useDynLib baquantreg
#' @examples
#' set.seed(1)

spBQR <- function(formula, tau = 0.5, data, itNum, thin=1,
                      betaValue = NULL, sigmaValue=1, spCoord1, spCoord2,
                      lambda = 0.5, tuneP = 1, alpha = 0.5, tuneA = 1e3,
                      priorVar = 100,
                      refresh = 100, quiet = T, jitter = 0,
                      includeAlpha = TRUE, tuneV = 0.5, kMT = 5){

  y <- as.numeric(stats::model.extract(stats::model.frame(formula, data), 'response'))
  X <- stats::model.matrix(formula, data)

  if (is.null(betaValue)) betaValue <- rep(0, dim(X)[2])

  spatial1 <- data[,spCoord1]
  spatial2 <- data[,spCoord2]

  matDist <- outer(spatial1, spatial1, "-")^2 +
    outer(spatial2, spatial2, "-")^2

  output <- list()

  output$chains <- lapply(tau, function(a){
    spBayesQR(tau = a, y = y, X = X, itNum = itNum, thin = thin,
            betaValue = betaValue, sigmaValue = sigmaValue,
            matDist = matDist, lambda = lambda,
            tuneP = tuneP, alphaValue = alpha, tuneA = tuneA,
            priorVar = priorVar, refresh = refresh, quiet = quiet,
            jitter = jitter, includeAlpha = includeAlpha,
            tuneV = tuneV, kMT = kMT)
  })

  output$acceptRateKappa <- lapply(output$chains, function(b){
    1-sum(diff(b$kappa1)==0)/(length(b$kappa1)-1)
  })

  output$acceptRateAlpha <- lapply(output$chains, function(b){
    1-sum(diff(b$alphaSample)==0)/(length(b$alphaSample)-1)
  })

  output$acceptRateV <- lapply(output$chains, function(b){
    1-apply(b$vSample, 2, function(bb){
      sum(diff(bb)==0)/length(diff(bb))
    })
  })

  output$tau <- tau
  output$formula <- formula
  output$data <- data

  class(output) <- "spBQR"
  return(output)
}
