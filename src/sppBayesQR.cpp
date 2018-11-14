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

   double theta, psi2, s0, n0, nTilde, sTilde, delta2;

   NumericVector sigmaSample(itNum), termsSum(n);

   arma::colvec b0(p), zSample(n), mu(p), resVec(n);

   arma::mat B0(p,p), SigmaMinusOne(p,p), diagU, Sigma,
    betaSample(itNum, p), vSample(itNum, n), sigmaDot(n, n);

   B0 = priorVar * B0.eye(p,p);
   b0.fill(0);

   theta = (1-2*tau)/(tau*(1-tau));
   psi2 = 2/(tau*(1-tau));

   zSample.fill(1);

   n0 = 3.0;
   s0 = 0.1;

  arma::mat covMatAux(n, m, arma::fill::zeros), covMat(n, n, arma::fill::zeros),
    covMat2(m, m, arma::fill::zeros), covMatInv(m, n, arma::fill::zeros),
    auxCov(m, m, arma::fill::zeros), cholCov(m,m), cholCov2(m,m),
    matAux(n,n, arma::fill::zeros), matM(m, m, arma::fill::zeros),
    matM2(m, n, arma::fill::zeros), matM3(n, n, arma::fill::zeros),
    CovCov(n, n, arma::fill::zeros);

  arma::colvec lambdaSample(itNum), alphaSample(itNum, arma::fill::zeros);

  lambdaSample[0] = lambda;
  if (includeAlpha) alphaSample[0] = alphaValue;
  else alphaValue = 0.0;

  NumericVector lambdaPrior = logPriorKappa2(lambdaVec, shapeL, rateL);

  IntegerVector seqRefresh = seq(1, itNum/refresh)*(refresh);

   for(int k = 1; k < itNum; k++) {
      for(int j = 0; j < thin; j++) {

        if(!quiet){
          if(is_true(any(k+1 == seqRefresh))){
            Rcout << "Iteration = " << k+1 << std::endl;
          }
        }

        covMat = exp(-lambda*matDist);

//         for(int aa = 0; aa < n; aa++)
//           for(int bb = 0; bb < n; bb++)
//             covMat(aa,bb) = exp(-lambda *
//               (pow(spCoord1[aa] - spCoord1[bb],2) +
//               pow(spCoord2[aa] - spCoord2[bb],2)));

        /*covmat2 = covmat.submat(indices, indices);
        covmataux = covmat.cols(indices);

        auxcov.diag().fill(jitter);

        cholcov = arma::chol((1-alphavalue)*covmat2 + auxcov, "lower");
        covmatinv = arma::solve(trimatl(cholcov), covmataux.t());
        mataux = covmatinv.t() * covmatinv;
        diagu = diagmat(sqrt(1/zsample));

        sigmadot = diagmat(alphavalue +
          (1-alphavalue)*(covmat.diag() - mataux.diag())).i();

        cholcov2 = (1-alphavalue)*covmat2 + covmataux.t()*sigmadot*covmataux;

        matm = arma::chol(cholcov2 + auxcov, "lower");
        matm2 = solve(trimatl(matm), covmataux.t());
        matm3 = matm2.t() * matm2;

        covcov = sigmadot - sigmadot * matm3 * sigmadot;

        sigmaminusone = ((1/(sigmavalue*psi2)) * x.t() * diagu * covcov *
          diagu * x) + b0.i();
        sigma = sigmaminusone.i();

        mu = sigma * ((1/(sigmavalue*psi2))*(x.t() * diagu * covcov *
                diagu * (y - theta*zsample)));*/

        betaValue = 21;

        resVec = y - theta*zSample - X * betaValue;

        /*nTilde = n0 + 3*n;
        sTilde =  arma::as_scalar(s0 + 2*sum(zSample) + resVec.t() * diagU *
          CovCov * diagU * resVec);*/

        sigmaValue = 18;

        for(int o = 0; o < n; o++){
          //zSample[o] = mtM(y - X * betaValue, theta, psi2, sigmaValue, zSample,
          //                 zSample[o], o, CovCov, tuneV, kMT);
			zSample[o] = 1;
        }

//         lambda = mhKappa2(lambda, spCoord1, spCoord2, resVec, diagU,
//                                covMat, CovCov,
//                                tuneP, alphaValue, jitter, indices, m);

        if (discLambda) {
          lambda = discKappa2(lambdaVec, lambdaPrior, matDist,
                              resVec, diagU, covMat, CovCov, alphaValue,
                              jitter, indices, m);
        }
        else {
          //lambda = mhKappa2(lambda, matDist, resVec, diagU,
           //                covMat, CovCov,
           //                tuneP, alphaValue, jitter, indices, m,
           //                shapeL, rateL);
			lambda = 1;
        }

        if (includeAlpha){
          alphaValue = mhAlpha2(alphaValue, resVec, diagU, covMat, covMat2,
                                covMatAux, tuneA, jitter, indices, m);
        }


      }

      betaSample.row(k) = betaValue.t();
      vSample.row(k) = zSample.t();
      sigmaSample[k] = sigmaValue;
      lambdaSample[k] = lambda;
      alphaSample[k] = alphaValue;
   }

   return Rcpp::List::create(
        Rcpp::Named("BetaSample") = betaSample,
        Rcpp::Named("SigmaSample") = sigmaSample,
        Rcpp::Named("vSample") = vSample,
        Rcpp::Named("LambdaSample") = lambdaSample,
        Rcpp::Named("alphaSample") = alphaSample);
}
