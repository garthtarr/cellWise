
#include "cellMCD.h"


arma::mat subinverse_cpp(const arma::mat& Sigma,
                         const arma::mat& Sigmai,
                         arma::uvec indx) {
  // # Returns inverse of a submatrix of the symmetric matrix Sigma.
  // # The submatrix is determined by the indices in indx.
  // # Sigma= original matrix, Sigmai = inverse of Sigma
  
  arma::uword p = Sigma.n_cols;
  arma::uword nindx = indx.n_elem;
  arma::uvec indx_neg = arma::regspace<arma::uvec>(0, p - 1);
  indx_neg.shed_rows(indx);
  arma::mat result(nindx, nindx); 
  if ((nindx < (arma::uword)(p)) &&
      (nindx > (arma::uword)(p / 2.0))){
    result = Sigmai(indx, indx) -
      Sigmai(indx, indx_neg) *
      arma::inv(Sigmai(indx_neg, indx_neg)) *
      Sigmai(indx_neg, indx);
  } else if((nindx < (arma::uword)(p)) &&
    (nindx <= (arma::uword)(p / 2.0))){ // don't use the trick above
    result = arma::inv(Sigma(indx,indx));
  } else { // here nindx=d, so submatrix is entire matrix
    result = Sigmai;
  }
  return(result);
}

void uniqueRows(arma::umat W,
                std::unordered_map<std::string, arma::uvec>& Wmap){
  //
  // constructs an unordered map of pairs (string, arma::vec)
  // string is the binary string resulting from concatenating the row in W
  // the second element is a vector containing the row indices with 
  // these patterns
  
  arma::uword n = W.n_rows;
  arma::uword p = W.n_cols;
  for(arma::uword i = 0; i < n; i++)
  {
    std::string s = "";
    
    for(arma::uword j = 0; j < p; j++)
      s += std::to_string(W(i,j));
    
    arma::uvec idx(1);
    idx(0) = i;
    if(Wmap.count(s) == 0)
    {
      Wmap.insert({s, idx});
    } else {
      arma::uvec idxold = Wmap.at(s); 
      arma::uvec idxnew = arma::join_cols(idxold, idx);
      Wmap[s] = idxnew;
    }
  }
}

arma::vec Deltacalc_cpp(const arma::mat & X,
                        const arma::umat & W,
                        const arma::mat & Sigma,
                        const arma::mat & Sigmai,
                        const arma::vec & mu,
                        arma::uword j){
  // convert to 0-based indexing, for calling directly from R
  j = j - 1;
  
  arma::uword n = X.n_rows;
  arma::uword p = X.n_cols;
  arma::vec deltas(n);
  deltas.fill(arma::datum::inf);
  
  std::unordered_map<std::string, arma::uvec> Wmap;
  uniqueRows(W, Wmap);
  
  arma::vec  x(n);
  arma::uvec inds;
  arma::uvec finiteInds;
  arma::urowvec w0(p);
  arma::urowvec w1(p);
  for (std::pair<std::string, arma::uvec> w : Wmap) {
    x = X.col(j);
    inds = w.second; // all indices with rowpattern w
    finiteInds = arma::find_finite(x(inds));
    
    if (finiteInds.n_elem > 0) {
      inds = inds(finiteInds);
      x = x(inds);
      
      w0 = W.row(inds(0));
      w1 = W.row(inds(0));
      w0(j) = 0;
      
      if (arma::any(w0)) {
        w1(j) = 1;
        arma::uvec index0 = arma::find(w0);
        arma::uvec index1 = arma::find(w1);
        
        arma::rowvec mu0  = mu(index0).t();
        arma::rowvec mu1  = mu(index1).t();
        
        arma::mat Sigma1  = Sigma(index1, index1);
        arma::mat Sigma0i = subinverse_cpp(Sigma, Sigmai, index0);
        
        arma::uvec jtemp   = arma::find(index1 == j);
        arma::uvec jtemp_n = arma::find(index1 != j);
        
        arma::mat Xtemp = X(inds, index0);
        Xtemp.each_row() -= mu0;
        arma::vec x1h = arma::as_scalar(mu1(jtemp)) +
          Xtemp * Sigma0i * Sigma1(jtemp_n, jtemp);
        
        double C1 = arma::as_scalar(Sigma1(jtemp, jtemp) -
                                    Sigma1(jtemp, jtemp_n) *
                                    Sigma0i * Sigma1(jtemp_n, jtemp));
        
        deltas(inds) = arma::pow(x - x1h, 2) /  C1 + 
          std::log(C1) + std::log(2 *  arma::datum::pi);
      } else {
        deltas(inds) = arma::pow(x - mu(j), 2) / (double) Sigma(j, j) + 
          std::log(Sigma(j, j)) + std::log(2 *  arma::datum::pi);
      }
    }
  }
  return(deltas);
}







arma::vec Deltacalc_cpp2(const arma::mat & X,
                         const arma::umat & W,
                         const arma::mat & Sigma,
                         const arma::mat & Sigmai,
                         const arma::vec & mu,
                         arma::uword j){
  // convert to 0-based indexing, for calling directly from R
  j = j - 1;
  
  arma::uword n = X.n_rows;
  arma::uword p = X.n_cols;
  arma::vec deltas(n);
  deltas.fill(arma::datum::inf);
  
  std::unordered_map<std::string, arma::uvec> Wmap;
  uniqueRows(W, Wmap);
  
  arma::vec  x(n);
  arma::uvec inds;
  arma::uvec finiteInds;
  arma::urowvec w0(p);
  arma::urowvec w1(p);
  for (std::pair<std::string, arma::uvec> w : Wmap) {
    x = X.col(j);
    inds = w.second; // all indices with rowpattern w
    finiteInds = arma::find_finite(x(inds));
    
    if (finiteInds.n_elem > 0) {
      inds = inds(finiteInds);
      x = x(inds);
      
      w0 = W.row(inds(0));
      w1 = W.row(inds(0));
      w0(j) = 0;
      
      if (arma::any(w0)) {
        w1(j) = 1;
        arma::uvec index0 = arma::find(w0);
        arma::uvec index1 = arma::find(w1);
        
        arma::uvec jtemp   = arma::find(index1 == j);
        arma::uvec jtemp_n = arma::find(index1 != j);
        
        double C1 = 1;
        arma::vec x1h(inds.n_elem, arma::fill::zeros);
        
        if (index0.n_elem <= (arma::uword)((double)p / 2.0)) {
          // compute v  as a solution to \Sigma_{index0,  index0} v = Sigma1_{jtemp_n, jtemp}
          arma::mat A = Sigma(index0, index0);
          arma::vec b = Sigma(index0, arma::uvec {j});
          arma::mat L = arma::chol(A, "lower");
          arma::vec v = arma::solve(arma::trimatu(L.t()),
                                    arma::solve(arma::trimatl(L), b));
          arma::mat Xtemp = X(inds, index0);
          Xtemp.each_row() -= mu(index0).t();
          x1h = mu(j) + Xtemp * v;
          C1 = Sigma(j,j) - arma::as_scalar(Sigma(arma::uvec {j}, index0) *  v);
        } else {
          // in this case we can use mu[m] - ((\sigma^{-1})_{m,m})^{-1} \sigma_{m,o} (x_o - mu_o)
          // this requires some attention though, since we need only 
          // the prediction for j, we need to extract the jth row of
          // ((\sigma^{-1})_{m,m})^{-1} \sigma_{m,o}
          
          arma::uvec index0_neg = arma::regspace<arma::uvec>(0, p - 1);
          index0_neg.shed_rows(index0);
          
          arma::mat A = Sigmai(index0_neg, index0_neg);
          arma::mat L = arma::chol(A, "lower");
          
          arma::vec ej = arma::zeros<arma::vec>(index0_neg.n_elem);
          arma::uvec tempidx = arma::find(index0_neg == j);
          ej(tempidx(0)) = 1.0;
          
          arma::vec Y = arma::solve(trimatu(L.t()),
                                    arma::solve(arma::trimatl(L), ej));
          arma::vec v = Sigmai(index0, index0_neg) * Y;
          
          arma::mat Xtemp = X(inds, index0);
          Xtemp.each_row() -= mu(index0).t();
          x1h = mu(j) - Xtemp * v;
          
          // now conditional variance
          
          C1 =  Y(tempidx(0)); 
        }
        
        deltas(inds) = arma::pow(x - x1h, 2) /  C1 + 
          std::log(C1) + std::log(2 *  arma::datum::pi);
        
      } else {
        deltas(inds) = arma::pow(x - mu(j), 2) / (double) Sigma(j, j) + 
          std::log(Sigma(j, j)) + std::log(2 *  arma::datum::pi);
      }
    }
  }
  return(deltas);
}
