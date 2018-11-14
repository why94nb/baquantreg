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

   RNGScope scope;

   int n = X.n_rows;
   int p = X.n_cols;

   
   return Rcpp::List::create(
        Rcpp::Named("BetaSample") = n,
        Rcpp::Named("SigmaSample") = p
 
}
