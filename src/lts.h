#ifndef LTS_H
#define LTS_H

#ifndef ARMA_DONT_PRINT_ERRORS
#define ARMA_DONT_PRINT_ERRORS
#endif

#ifndef  ARMA_USE_CXX11
#define ARMA_USE_CXX11
#endif

#include "RcppArmadillo.h" 
#include "LocScaleEstimators.h"

struct CstepOut {
  double alpha;
  arma::vec beta;
  arma::uvec idx;
  double objective;
};


double getRandomStart(const arma::mat& X, const arma::vec& y,
                      arma::uvec& idx, 
                      const double lambda, arma::uword h);

CstepOut executeCsteps(const arma::mat& X, const arma::vec& y,
                       arma::uvec& idx,
              const arma::uword maxIts,
              const double precScale, const double lambda);

std::tuple<double, double> updateIdx(arma::vec residuals, 
                                     arma::uvec& idx, arma::uword h);

Rcpp::List lts_cpp(const arma::mat& X, const arma::vec& y, const double alpha,
                   const arma::uword maxIts,
                   const arma::uword nFinal, 
                   const double precScale);


arma::vec getCoef(const arma::mat& X, const arma::vec& y,
                  const arma::uvec& idx,
                  double lambda);


#endif
