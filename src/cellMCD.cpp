
#include "cellMCD.h"

void uniqueRows(const arma::umat& W,
                 std::unordered_map<std::string,std::vector<arma::uword>>& Wmap){
  //
  // constructs an unordered map of pairs (string, arma::vec)
  // string is the binary string resulting from concatenating the row in W
  // the second element is a vector containing the row indices with 
  // these patterns
  arma::uword n = W.n_rows;
  arma::uword p = W.n_cols;
  
  std::string s;
  s.reserve(p);
  
  for(arma::uword i = 0; i < n; i++) {
    s.clear();
    for(arma::uword j = 0; j < p; j++) {
      s += (char)('0' + W(i,j));
    }
    
    Wmap[s].push_back(i);  
  }
}

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
  
  
  
  std::unordered_map<std::string, std::vector<arma::uword>> Wmap;
  uniqueRows(W, Wmap);
  
  arma::vec  x(n);
  arma::uvec inds;
  arma::uvec finiteInds;
  arma::urowvec w0(p);
  arma::urowvec w1(p);
  for (std::pair<std::string, arma::uvec> w : Wmap) {
    x = X.col(j);
    inds = arma::conv_to<arma::uvec>::from(w.second); // all indices with rowpattern w
    finiteInds = arma::find_finite(x(inds));
    
    if (finiteInds.n_elem > 0) {
      inds = inds(finiteInds);
      x = x(inds);
      
      w0 = W.row(inds(0));
      w0(j) = 0;
      
      if (arma::any(w0)) {
        w1 = W.row(inds(0));
        w1(j) = 1;
        arma::uvec index0 = arma::find(w0);
        arma::uvec index1 = arma::find(w1);
        
        double C1 = 1;
        arma::vec x1h(inds.n_elem, arma::fill::zeros);
        
        if (index0.n_elem <= (arma::uword)((double)p / 2.0)) {
          arma::mat A = Sigma(index0, index0);
          arma::vec b = Sigma(index0, arma::uvec {j});
          
          arma::vec v(index0.n_elem);
          arma::solve(v, A, b,  arma::solve_opts::fast + arma::solve_opts::likely_sympd);
          x1h = mu(j) - arma::dot(mu(index0), v) + X(inds, index0) * v;
          
          C1 = Sigma(j,j) - arma::as_scalar(Sigma(arma::uvec {j}, index0) *  v);
        } else {
          // in this case we can use mu[m] - ((\sigma^{-1})_{m,m})^{-1} \sigma_{m,o} (x_o - mu_o)
          // this requires some attention though, since we need only 
          // the prediction for j, we need to extract the jth row of
          // ((\sigma^{-1})_{m,m})^{-1} \sigma_{m,o}
          
          arma::uvec index0_neg = arma::regspace<arma::uvec>(0, p - 1);
          index0_neg.shed_rows(index0);
          
          arma::mat A = Sigmai(index0_neg, index0_neg);
          
          arma::vec ej = arma::zeros<arma::vec>(index0_neg.n_elem);
          arma::uvec tempidx = arma::find(index0_neg == j);
          ej(tempidx(0)) = 1.0;
          arma::vec Y(index0_neg.n_elem);
          arma::solve(Y, A, ej,  arma::solve_opts::fast +arma::solve_opts::likely_sympd);
          
          arma::vec v = Sigmai(index0, index0_neg) * Y;
          
          
          // 
          x1h = mu(j) + arma::dot(mu(index0), v) - X(inds, index0) * v;
          
          
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




// [[Rcpp::export]]
Rcpp::List allpreds_cpp(arma::mat& X,
                         arma::mat& S,
                         arma::vec& mu,
                         arma::umat& W) {
  try
  {
    
    arma::uword n = W.n_rows;
    arma::uword d = W.n_cols;
    //arma::mat Sinv = arma::inv_sympd(arma::symmatu(S));
    arma::mat Sinv = arma::inv(S);
    arma::mat preds(n, d);
    arma::mat cvars(n, d);
    
    
    std::unordered_map<std::string,std::vector<arma::uword>> Wmap;
    uniqueRows(W, Wmap);
    
    
    for (std::pair<std::string, std::vector<arma::uword>> s : Wmap) {
      // for (std::pair<std::string, arma::uvec> s : Wmap) {
      
      arma::uvec myrows = arma::conv_to<arma::uvec>::from(s.second); // all indices with rowpattern(string) s
      arma::urowvec w = W.row(myrows(0)); // the actual rowpattern
      
      arma::uvec obs = arma::find(w);
      arma::uvec mis = arma::find(1 - w);
      
      if (mis.n_elem < d) {
        if (mis.n_elem > 0) {
          
          arma::mat Sjinv = subinverse_cpp(S, Sinv, obs);
          
          arma::vec cvar = arma::diagvec(S(mis, mis) -
            S(mis, obs) * Sjinv * S(obs, mis));
          for (arma::uword i = 0; i < myrows.n_elem; i++) {
            arma::uvec rowi(1);
            rowi(0) = myrows(i);
            cvars(rowi, mis) = cvar.t();
            preds(rowi, mis) = mu(mis).t() + (X(rowi, obs) - mu(obs).t()) * Sjinv * S(obs, mis);
          }
          
          // now compute preds and cvars for all available entries
          //
          
          if (obs.n_elem == 1) {
            for (arma::uword i = 0; i < myrows.n_elem; i++) {
              arma::uvec rowi(1);
              rowi(0) = myrows(i);
              preds(rowi, obs) = mu(obs);
              cvars(rowi, obs) = S(obs, obs);
            }
          } else {
            for (arma::uword j = 0; j < obs.n_elem; j++) {
              arma::uvec obsj(1);
              obsj(0) = obs(j);
              arma::uvec helps = obs(arma::find(obs != obs(j)));
              arma::mat M2inv = subinverse_cpp(S, Sinv, helps);
              
              arma::mat cvar = S(obsj, obsj) - S(obsj, helps) *
                M2inv * S(helps, obsj);
              
              for (arma::uword i = 0; i < myrows.n_elem; i++) {
                arma::uvec rowi(1);
                rowi(0) = myrows(i);
                cvars(rowi, obsj) = cvar;
                preds(rowi, obsj) = mu(obsj).t() +
                  (X(rowi, helps) - mu(helps).t()) * M2inv * S(helps, obsj);
              }
            }
          }
          
        } else { //no missings in this pattern
          
          arma:: vec cvar = (1 / arma::diagvec(Sinv));
          
          for (arma::uword i = 0; i < myrows.n_elem; i++) {
            cvars.row(myrows(i)) =  cvar.t();
          }
          
          for (arma::uword j = 0; j < d; j++) {
            arma::uvec jvec(1);
            jvec(0) = j;
            arma::uvec j_neg = arma::regspace<arma::uvec>(0, d - 1);
            j_neg.shed_rows(jvec);
            
            arma::mat Sjinv = subinverse_cpp(S, Sinv, j_neg);
            arma::mat Xtemp = X(myrows, j_neg);
            Xtemp.each_row() -= mu(j_neg).t();
            
            preds(myrows, jvec) = (arma::as_scalar(mu(j)) +
              Xtemp * Sjinv * S(j_neg, jvec));
          }
        }
      } else { // all elements of this pattern are missing
        for (arma::uword i = 0; i < myrows.n_elem; i++) {
          preds.row(myrows(i)) = mu.t();
          cvars.row(myrows(i)) =  arma::diagvec(Sinv).t();
        }
      }
      
    }
    
    
    return Rcpp::List::create(Rcpp::Named("preds") = preds,
                              Rcpp::Named("cvars") = cvars);
  } catch( std::exception& __ex__ )
  {
    forward_exception_to_r( __ex__ );
  } catch(...)
  {
    ::Rf_error( "c++ exception " "(unknown reason)" );
  }
  
  return Rcpp::wrap(NA_REAL);
}





// [[Rcpp::export]]
double Objective_cpp(const arma::mat &X,
                      const arma::umat &W,
                      const arma::vec &mu,
                      const arma::mat &Sigma,
                      const arma::mat &Sigmai) {
  try
  {
    
    double objective = 0;
    double objective_add = 0;
    
    
    
    std::unordered_map<std::string, std::vector<arma::uword>> Wmap;
    uniqueRows(W, Wmap);
    
    
    for (std::pair<std::string, std::vector<arma::uword>> s : Wmap) {
      
      arma::uvec myrows = arma::conv_to<arma::uvec>::from(s.second); // all indices with rowpattern(string) s
      arma::urowvec w = W.row(myrows(0)); // the actual rowpattern
      
      arma::uvec obs = arma::find(w);
      arma::uvec mis = arma::find(1 - w);
      
      arma::vec subMu = mu(obs);
      arma::mat subSigma = Sigma(obs, obs);
      arma::mat invsubSigma = subinverse_cpp(Sigma, Sigmai, obs);
      
      
      arma::mat Xtemp = X(myrows, obs);
      Xtemp.each_row() -= subMu.t();
      arma::mat MD = arma::sum((Xtemp * invsubSigma) % Xtemp, 1);
      objective_add = arma::as_scalar(arma::sum(MD + 
        arma::log_det_sympd(arma::symmatu(subSigma)) + 
        std::log(2 * arma::datum::pi) * (obs.n_elem)));
      objective    = objective + objective_add;
    } 
    
    return(objective);
  } catch( std::exception& __ex__ )
  {
    forward_exception_to_r( __ex__ );
  } catch(...)
  {
    ::Rf_error( "c++ exception " "(unknown reason)" );
  }
  
  return NA_REAL;
}



// [[Rcpp::export]]
arma::umat updateW_cpp(const arma::mat &X,
                       arma::umat W,
                       const arma::vec &mu,
                       const arma::mat &Sigma,
                       const arma::mat &Sigmai,
                       const arma::vec &lambda,
                       const arma::uword &h) {
  try
  {
    
    arma::uword n = W.n_rows;
    arma::uword d = W.n_cols;
    arma::uvec ord = arma::stable_sort_index(arma::sum(W, 0));
    arma::vec Delta(n);
    arma::uvec goodCells;
    arma::uvec wnew(n);
    arma::uvec Delta_rank(n);
    for (arma::uword j = 0; j < d; j++) {
      arma::uword jtemp = ord(j);
      Delta = Deltacalc_cpp(X, W, Sigma, Sigmai, mu, jtemp + 1);
      goodCells = arma::find(Delta <= lambda(jtemp));
      if (goodCells.n_elem < h) {
        Delta_rank = arma::stable_sort_index(Delta - lambda(jtemp));
        goodCells = Delta_rank.head(h);
      }
      wnew.zeros();
      wnew(goodCells).fill(1);
      W.col(jtemp) = wnew;
    }
    return(W);
  } catch( std::exception& __ex__ )
  {
    forward_exception_to_r( __ex__ );
  } catch(...)
  {
    ::Rf_error( "c++ exception " "(unknown reason)" );
  }
  arma::umat error_out(1,1);
  error_out.fill(0);
  return error_out;
}





// [[Rcpp::export]]
Rcpp::List iterMCD_cpp(const arma::mat& X,
                   const arma::vec& initmu,
                   const arma::mat& initSigma,
                   const arma::uword h, // ceil(alpha*n)
                   const arma::vec& lambdas,
                   const double crit,
                   const arma::uword noCits,
                   const double lmin,
                   const double precScale,
                   const bool fixedCenter,
                   const bool silent) {
  
  try {
    arma::uword d = X.n_cols;
    arma::uword n = X.n_rows;
    
    arma::umat W(n, d, arma::fill::ones);
    W(arma::find_nonfinite(X)).zeros();
    
    arma::vec objvals;
    
    
    arma::vec mu = initmu;
    arma::mat Sigma = initSigma;
    arma::mat Sigmai;
    truncEig(Sigma, Sigmai, lmin);
    
    
    double convcrit = crit + 1.0;
    arma::uword nosteps = 0;
    arma::mat oldSigma = Sigma;
    
    if (!silent) {
      objvals.resize(noCits + 1);
      objvals.fill(arma::datum::nan);
      double penalty = arma::sum(arma::sum(1 - W).t() % lambdas);
      objvals(0) = Objective_cpp(X, W, initmu, initSigma, Sigmai) + penalty;
      Rcpp::Rcout << "Objective at step " << nosteps << " = " << std::round(objvals(0) * 10000) / 10000 << std::endl;
    } 
    
    
    
    while (convcrit > crit && nosteps < noCits){
      Cstep(X, W, mu, Sigma, Sigmai, lambdas, 
             h, fixedCenter);
      
      // now truncate the eigenvalues and compute Sigmai
      truncEig(Sigma, Sigmai, lmin);
      
      if (!silent) {
        double penalty = arma::sum(arma::sum(1 - W).t() % lambdas);
        objvals(nosteps + 1) = Objective_cpp(X, W, initmu, Sigma, Sigmai) + penalty;
        Rcpp::Rcout << "Objective at step " << nosteps << " = " << std::round(objvals(nosteps + 1) * 10000) / 10000 << std::endl;
      }
      
      
      convcrit = (arma::abs(Sigma - oldSigma)).max();
      oldSigma = Sigma;
      nosteps++;
    }
    
    return Rcpp::List::create(Rcpp::Named("W") = W,
                              Rcpp::Named("mu") = mu,
                              Rcpp::Named("Sigma") = Sigma,
                              Rcpp::Named("nosteps") = nosteps);
    
  } catch( std::exception& __ex__ )
  {
    forward_exception_to_r( __ex__ );
  } catch(...)
  {
    ::Rf_error( "c++ exception " "(unknown reason)" );
  }
  return Rcpp::List::create(Rcpp::Named("error") = true);
}







void Cstep(const arma::mat &X,
            arma::umat &W,
            arma::vec &mu,
            arma::mat &Sigma,
            const arma::mat &Sigmai,
            const arma::vec &lambda,
            const arma::uword &h,
            const bool fixedCenter) {
  // W, mu and Sigma are updated in-place.
  // returns the "convcrit"= max(abs((oldSigma - newSigmaSigma))) 
  // First update W
  
  arma::uword n = W.n_rows;
  arma::uword d = W.n_cols;
  arma::uvec ord = arma::stable_sort_index(arma::sum(W, 0));
  arma::vec Delta(n);
  arma::uvec goodCells;
  arma::uvec wnew(n);
  arma::uvec Delta_rank(n);
  for (arma::uword j = 0; j < d; j++) {
    arma::uword jtemp = ord(j);
    Delta = Deltacalc_cpp(X, W, Sigma, Sigmai, mu, jtemp + 1);
    goodCells = arma::find(Delta <= lambda(jtemp));
    if (goodCells.n_elem < h) {
      Delta_rank = arma::stable_sort_index(Delta - lambda(jtemp));
      goodCells = Delta_rank.head(h);
    }
    wnew.zeros();
    wnew(goodCells).fill(1);
    W.col(jtemp) = wnew;
  }
  
  
  
  // Now execute one EM step:
  // first compute Ximp
  
  arma::mat Ximp = X;
  
  std::unordered_map<std::string, std::vector<arma::uword>> Wmap;
  uniqueRows(W, Wmap);
  
  arma::mat bias(d, d, arma::fill::zeros);
  
  for (std::pair<std::string, std::vector<arma::uword>> s : Wmap) {
    
    
    arma::uvec myrows = arma::conv_to<arma::uvec>::from(s.second); // all indices with rowpattern(string) s
    arma::urowvec w = W.row(myrows(0)); // the actual rowpattern
    
    arma::uvec obs = arma::find(w);
    arma::uvec mis = arma::find(1 - w);
    
    if (mis.n_elem < d) {
      if (mis.n_elem > 0) {
        
        arma::mat Sigmai_temp = subinverse_cpp(Sigmai, Sigma, mis); //solve(Sigmai[mis, mis])  
        bias(mis,mis) += (Sigmai_temp * (double)myrows.n_elem);
        
        for (arma::uword i = 0; i < myrows.n_elem; i++) {
          arma::uvec rowi(1);
          rowi(0) = myrows(i);
          Ximp(rowi, mis) = mu(mis).t() - (Sigmai_temp * Sigmai(mis, obs) *(X(rowi, obs).t() - mu(obs))).t();
        }
      }
    } else { // all elements of this pattern are missing
      bias += (Sigma * (double)myrows.n_elem);
      for (arma::uword i = 0; i < myrows.n_elem; i++) {
        Ximp.row(myrows(i)) = mu.t();
      }
    }
  }
  
  arma::vec newmu = mu;
  if (!fixedCenter) {
    newmu = mean(Ximp, 0).t();
  }
  
  bias /= (double)n;
  
  arma::mat newSigma = bias;
  if (!fixedCenter) {
    newSigma += arma::cov(Ximp, 1);
  } else {
    newSigma += (Ximp.t() * Ximp) / (double)n ;
  }
  
  mu = newmu;
  Sigma = newSigma;
}


void truncEig(arma::mat& Sigma,
              arma::mat& Sigmai,
              const double lmin) {
  // Truncates the eigenvalues of S from below by lmin
  
  // first ensure perfect symmetry
  Sigma += Sigma.t();
  Sigma *= 0.5;
  
  arma::mat U;
  arma::vec s;
  arma::mat V;
  arma::svd_econ(U, s, V, Sigma); // could also use eig_sym, but this might be slightly more stable.
  
  if (!std::isnan(lmin)) {
    s.transform([lmin](double v) { return std::max(v, lmin); });
  }
  
  U.each_row() %= arma::sqrt(s).t();
  Sigma = U * U.t();
  
  // Compute inverse: U * diag(1/s) * U'
  // Reuse U but scale by 1/s instead
  U.each_row() %= (1.0 / (s)).t();
  Sigmai = U * U.t();
  
}
