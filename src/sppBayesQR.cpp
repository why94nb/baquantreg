#include <RcppArmadillo.h>
#include <RcppGSL.h>

#include "helperLV.h"
#include "helperRD.h"
#include "helperAlpha.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List sppBayesQR(double tau, arma::colvec y, arma::mat X, int itNum,
                   int thin, arma::colvec betaValue, double sigmaValue,
                   arma::mat matDist,
                   NumericVector lambdaVec, double lambda,
                   double shapeL, double rateL,
                   double tuneP, arma::uvec indices, int m,
                   double alphaValue, double tuneA, double priorVar,
                   bool quiet, int refresh, double jitter, bool includeAlpha,
                   double tuneV, int kMT, bool discLambda){


   int n = X.n_rows;
   int p = X.n_cols;

   double theta, psi2, s0, n0, nTilde, sTilde, delta2;

   NumericVector sigmaSample(itNum), termsSum(n);

   arma::colvec b0(p), zSample(n), mu(p), resVec(n);

   arma::mat B0(p, p), SigmaMinusOne(p, p), diagU, Sigma,
	   betaSample(itNum, p), vSample(itNum, n), sigmaDot(n, n);

   B0 = priorVar * B0.eye(p, p);
   b0.fill(0);

   theta = (1 - 2 * tau) / (tau*(1 - tau));
   psi2 = 2 / (tau*(1 - tau));

   zSample.fill(1);

   n0 = 3.0;
   s0 = 0.1;

   arma::mat covMatAux(n, m, arma::fill::zeros), covMat(n, n, arma::fill::zeros),
	   covMat2(m, m, arma::fill::zeros), covMatInv(m, n, arma::fill::zeros),
	   auxCov(m, m, arma::fill::zeros), cholCov(m, m), cholCov2(m, m),
	   matAux(n, n, arma::fill::zeros), matM(m, m, arma::fill::zeros),
	   matM2(m, n, arma::fill::zeros), matM3(n, n, arma::fill::zeros),
	   CovCov(n, n, arma::fill::zeros);

   arma::colvec lambdaSample(itNum), alphaSample(itNum, arma::fill::zeros);

   lambdaSample[0] = lambda;
   if (includeAlpha) alphaSample[0] = alphaValue;
   else alphaValue = 0.0;

   NumericVector lambdaPrior = logPriorKappa2(lambdaVec, shapeL, rateL);

   IntegerVector seqRefresh = seq(1, itNum / refresh)*(refresh);

   
   return Rcpp::List::create(
	   Rcpp::Named("BetaSample") = n,
   Rcpp::Named("SigmaSample") = p);
   
}
