#ifndef cellMCD_H
#define cellMCD_H

#ifndef ARMA_DONT_PRINT_ERRORS
#define ARMA_DONT_PRINT_ERRORS
#endif

#ifndef  ARMA_USE_CXX11
#define ARMA_USE_CXX11
#endif

#include "RcppArmadillo.h" 

arma::mat subinverse_cpp(const arma::mat& Sigma,
                         const arma::mat& Sigmai,
                         arma::uvec indx);

void uniqueRows(const arma::umat& W,
                std::unordered_map<std::string,std::vector<arma::uword>>& Wmap);

arma::vec Deltacalc_cpp(const arma::mat & X,
                        const arma::umat & W,
                        const arma::mat & Sigma,
                        const arma::mat & Sigmai,
                        const arma::vec & mu,
                        arma::uword j);


void Cstep(const arma::mat &X,
              arma::umat &W,
              arma::vec &mu,
              arma::mat &Sigma,
              const arma::mat &Sigmai,
              const arma::vec &lambda,
              const arma::uword &h,
              const bool fixedCenter);

void truncEig(arma::mat& Sigma,
              arma::mat& Sigmai,
              double lmin);
#endif
