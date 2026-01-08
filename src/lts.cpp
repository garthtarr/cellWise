
#include "lts.h"






arma::mat getRandomStart2(const arma::mat& X, const arma::vec& y,
                          const double lambda, arma::uword h,
                          arma::uword nrep) {
  // 
  
  arma::mat betas(X.n_cols, 10, arma::fill::zeros);
  arma::uword n = X.n_rows;
  arma::uvec subsetIdx(h);
  arma::vec beta;
  arma::vec residuals(n);
  arma::uword nsuccess = 0;
  
  std::tuple<double, double> alfaObjective;
  
  
  double newObjective = arma::datum::inf;
  double maxObjective = arma::datum::inf;
  arma::vec ObjVals(10); 
  ObjVals.fill(arma::datum::inf);
  
  for (arma::uword i = 0; i < nrep; i++) { // take a random subset of size h and execute 2 Csteps
    subsetIdx =  arma::randperm(n, h);  
    beta = getCoef(X, y, subsetIdx,  lambda);
    
    if (std::isnan(beta(0))){
      continue;
    }
    
    residuals = y - X * beta;
    alfaObjective = updateIdx(residuals, subsetIdx, h);
    
    beta = getCoef(X, y, subsetIdx,  lambda);
    
    if (std::isnan(beta(0))){
      continue;
    }
    
    residuals = y - X * beta;
    alfaObjective = updateIdx(residuals, subsetIdx, h);
    beta = getCoef(X, y, subsetIdx,  lambda);
    
    if (std::isnan(beta(0))){
      continue;
    }
    alfaObjective = updateIdx(residuals, subsetIdx, h);
    
    // replace the worst beta by the new beta in case of a better objective:
    newObjective = std::get<1>(alfaObjective);
    if (newObjective < maxObjective) {
      arma::uword maxIdx = arma::index_max(ObjVals);
      ObjVals(maxIdx)    = newObjective;
      betas.col(maxIdx)  = beta;
      maxObjective       = arma::max(ObjVals);
      nsuccess++;
    }
  }
  if (nsuccess < 10) {
    Rcpp::stop("Could not find enough non-singular subsets to generate starting values.");
  }
  
  return(betas);
}



// [[Rcpp::export]]
Rcpp::List lts_cpp(const arma::mat& X,
                   const arma::vec& y,
                   const double alpha,
                   const double lambda,
                   const arma::uword nrep,
                   const arma::uword maxIts,
                   const arma::uword nFinal, // number of final candidate solutions to collect before iterating until convergence
                   const double precScale,
                   const arma::uword h) {
  // Lts with a ride penalty lambda * ||\beta||_2^2
  try{
    
    // align the RNG with R
    Rcpp::RNGScope scope;
    
    arma::uword n = X.n_rows;
    arma::uword d = X.n_cols;
    arma::uvec idx(h);
    // retain 10 best candidate solutions
    arma::umat bestIdx(h, nFinal);
    arma::vec ObjVals(nFinal); 
    ObjVals.fill(arma::datum::inf);
    
    
    
    if (X.n_rows <= 1500) { // n <= 1500
      
      double newObjective = arma::datum::inf;
      double maxObjective = arma::datum::inf;
      
      for (arma::uword i = 0; i < nrep; i++) {
        newObjective = getRandomStart(X, y, idx, lambda, h);
        if (newObjective < maxObjective) {
          arma::uword maxIdx  = arma::index_max(ObjVals);
          ObjVals(maxIdx)     = newObjective;
          bestIdx.col(maxIdx) = idx;
          maxObjective        = arma::max(ObjVals);
        }
      }
      
    } else { // n > 1500
      
      // the first 300 are the first subset, the next the second, etc.
      arma::uvec subsetIdx = arma::randperm(n, 1500);
      
      
      // within each subset of size 300, we need 100 replications
      // of a draw of a random subset of size hsub, do 2 Csteps
      // and keep the best 10 results
      arma::uword hsub = (arma::uword) (300 * (double) h / (double) n);
      arma::uword nrep_local = (arma::uword)((double)nrep / 5);
      arma::mat bestBeta(d, 50, arma::fill::zeros);
      
      arma::mat localBestBeta(d, 10, arma::fill::zeros);
      
      
      for (arma::uword i = 0; i < 5; i ++) {
        arma::uvec localsubset = subsetIdx.subvec(i * 300, (i + 1) * 300 - 1);
        arma::mat X_local = X.rows(localsubset);
        arma::vec y_local = y.elem(localsubset);
        
        arma::mat localBestBeta = getRandomStart2(X_local, y_local, lambda, hsub, nrep_local);
        
        bestBeta.cols(i * 10, (i + 1) * 10 - 1) = localBestBeta;
      }
      
      
      // now within the merged set of 1500, do 2 C steps for the 50 current solutions
      hsub = (arma::uword) (1500 * (double) h / (double) n);
      arma::mat X_local = X.rows(subsetIdx);
      arma::vec y_local = y.elem(subsetIdx);
      
      arma::vec beta(d);
      arma::vec residuals(1500);
      std::tuple<double, double> alfaObjective;
      localBestBeta.reshape(d, nFinal);
      localBestBeta.zeros();
      double newObjective = arma::datum::inf;
      double maxObjective = arma::datum::inf;
      
      for (arma::uword i = 0; i < bestBeta.n_cols; i++) { // take a random subset of size h and execute 2 Csteps
        
        beta = bestBeta.col(i);
        
        if (std::isnan(beta(0))){
          continue;
        }
        
        residuals = y_local - X_local * beta;
        alfaObjective = updateIdx(residuals, subsetIdx, hsub);
        beta = getCoef(X, y, subsetIdx,  lambda);
        
        if (std::isnan(beta(0))){
          continue;
        }
        
        residuals = y_local - X_local * beta;
        alfaObjective = updateIdx(residuals, subsetIdx, hsub);
        beta = getCoef(X, y, subsetIdx,  lambda);
        
        if (std::isnan(beta(0))){
          continue;
        }
        alfaObjective = updateIdx(residuals, subsetIdx, hsub);
        
        // replace the worst beta by the new beta in case of a better objective:
        newObjective = std::get<1>(alfaObjective);
        if (newObjective < maxObjective) {
          arma::uword maxIdx = arma::index_max(ObjVals);
          ObjVals(maxIdx)    = newObjective;
          localBestBeta.col(maxIdx)  = beta;
          maxObjective       = arma::max(ObjVals);
        }
      }
      
      // finally compute the h-subset for the nFinal best betas so far.
      for (arma::uword i = 0; i < nFinal; i++){
        residuals = y - X * beta;
        alfaObjective = updateIdx(residuals, idx, h);
        bestIdx.col(i) = idx;
      }
      
    }
    
    // final step: iterate through the 10 candidate subsets
    // and refine them until convergence.
    
    double bestObjective = arma::datum::inf;
    CstepOut bestCstepOut;
    
    for (arma::uword i = 0; i < nFinal; i++) {
      idx = bestIdx.col(i);
      CstepOut result = executeCsteps(X, y, idx, maxIts, precScale, lambda);
      if (result.objective < bestObjective) {
        bestCstepOut = result;
      }
    }
    
    return Rcpp::List::create(Rcpp::Named("idx") = bestCstepOut.idx,
                              Rcpp::Named("intercept") = bestCstepOut.alpha,
                              Rcpp::Named("beta") = bestCstepOut.beta);
    
  } catch( std::exception& __ex__ )
  {
    forward_exception_to_r( __ex__ );
  } catch(...)
  {
    ::Rf_error( "c++ exception " "(unknown reason)" );
  }
  return Rcpp::wrap(NA_REAL);
}




double getRandomStart(const arma::mat& X, const arma::vec& y,
                      arma::uvec& idx, 
                      const double lambda, arma::uword h) {
  // fill in random start indices in idx
  // they are generated from p+1 observations on which an ols is fit,
  // then the smallest h residuals are put in &idx.
  // if the subset of p+1 observations leads to a rank-deficient system,
  // observations are added one by one until it doesn't.
  // returns the rss of this subset
  
  arma::uvec subsetIdx = arma::randperm(X.n_rows);
  arma::vec beta(X.n_cols, arma::fill::zeros);
  
  arma::uword i = X.n_cols + 1;
  
  while (i < X.n_rows) {
    
    beta = getCoef(X, y, subsetIdx.head(i),  lambda);
    
    if (!std::isnan(beta(0))) {
      break;
    }
  }
  
  std::tuple<double, double> alfaObjective;
  
  if (!std::isnan(beta(0))){
    arma::vec residuals = y - X * beta;
    alfaObjective = updateIdx(residuals, idx, h);
  } else {
    Rcpp::Rcout << "Could not find a suitable starting value" << std::endl;
  }
  
  return(std::get<1>(alfaObjective));
}





CstepOut executeCsteps(const arma::mat& X, const arma::vec& y,
                       arma::uvec& idx,
                       const arma::uword maxIts,
                       const double precScale, const double lambda) {
  // idx: initial subset
  
  arma::vec beta(X.n_cols, arma::fill::zeros);
  double objective = arma::datum::inf;
  
  arma::vec residuals(X.n_rows, arma::fill::zeros);
  std::tuple<double, double> alfaObjective{0,0}; 
  
  bool converged = false;
  const arma::uword h = idx.n_elem;
  
  for (arma::uword i = 0; i < maxIts; i++) {
    
    beta = getCoef(X, y, idx, lambda);
    
    if (std::isnan(beta(0))) {
      break;
    }
    
    residuals = y - X * beta;
    
    // update idx, get intercept and objective function
    alfaObjective = updateIdx(residuals, idx, h);
    
    converged = (objective - std::get<1>(alfaObjective)) < precScale;
    if (converged) {
      break;
    }
    objective = std::get<1>(alfaObjective);
  }
  
  double alpha = std::get<0>(alfaObjective);
  
  
  CstepOut output;
  output.alpha = alpha;
  output.beta = beta;
  output.objective = objective;
  output.idx = idx;
  
  return(output);
}


std::tuple<double, double> updateIdx(arma::vec residuals, 
                                     arma::uvec& idx,
                                     arma::uword h) {
  // updates the idx of the smallest h subset.
  // resturns a tuple with (intercept, RSS)
  // wts is not used yet.
  //
  
  //determine quan
  arma::uword quan = h;
  arma::uword nbHSubsets = residuals.n_elem - quan + 1; // number of h-subsets
  
  arma::uvec I = arma::sort_index(residuals);
  residuals = residuals(I);
  arma::uword bestSubset = 0;
  
  double sh_0 = arma::sum(residuals.head(quan));
  double sh2_0 = std::pow(sh_0, 2) / (double)quan;
  double sq_0 = arma::as_scalar(arma::sum(arma::pow(residuals.head(quan), 2)) - sh2_0);
  double sh_1 = 0;
  double sh2_1 = 0;
  double sq_1 = 0;
  double sh_b = sh_0;
  double sq_b = sq_0;
  
  for (arma::uword i = 1; i < nbHSubsets ; i++) {
    sh_1 = sh_0 - residuals(i - 1) + residuals(i + quan - 1);
    sh2_1 = std::pow(sh_1, 2) / (double)quan;
    sq_1 = sq_0 - std::pow(residuals(i - 1), 2) + std::pow(residuals(i + quan - 1), 2) - sh2_1 + sh2_0;
    
    if (sq_1 < sq_b) {
      sh_b = sh_1;
      bestSubset = i;
      sq_b = sq_1;
    }
    sh_0 = sh_1;
    sh2_0 = sh2_1;
    sq_0 = sq_1;
  }
  
  double mean = sh_b / (double)quan;
  double rss  = sq_b;
  
  idx = I.subvec(bestSubset, static_cast<arma::uword>(bestSubset + h - 1));
  
  std::tuple<double, double> output(mean, rss);
  
  return(output);
}



arma::vec getCoef(const arma::mat& X,
                  const arma::vec& y, 
                  const arma::uvec& idx,
                  double lambda) {
  // Get coefficients of weighted least squares 
  // on the observations with indices in idx.
  // lambda is the size of the ridge penalty.
  // if the system is rank-defficient, the first
  // element of the returned coefficient vector is nan
  //
  
  
  // Subset the data and weights based on the provided indices
  arma::mat X_sub = X.rows(idx);
  arma::vec y_sub = y.rows(idx);
  
  arma::rowvec mu_X = arma::mean(X_sub, 0);
  double mu_y = arma::mean(y_sub);
  
  // center
  X_sub.each_row() -= mu_X;
  y_sub -= mu_y;
  
  
  // Compute Normal Equations using the weighted data (more stable and efficient)
  arma::mat H = X_sub.t() * X_sub;
  
  if (lambda > 0) {
    H.diag() += lambda;
  }
  
  const arma::vec G = X_sub.t() * y_sub;
  
  // Solve the linear system H * beta = G
  // Armadillo's solve uses numerically stable decomposition methods
  // (e.g., Cholesky, QR) appropriate for the matrix structure.
  
  arma::vec beta(X.n_cols, arma::fill::zeros);
  bool success = arma::solve(beta, H, G, arma::solve_opts::likely_sympd);
  
  if (!success) {
    beta(0) = arma::datum::nan;
  }
  
  return beta;
}
