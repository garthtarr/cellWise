

cellMCD <- function(X, alpha = 0.75, quant = 0.99, crit = 1e-4, 
                    noCits = 100, lmin = 1e-4, fixedCenter = FALSE,
                    checkPars = list()) {
  #
  # checkPars is as in cellWise::transfo(). If not coreOnly,
  # we run checkDataSet() like we did in other functions.
  #
  #
  #
  # Arguments:
  # X         : A n by d data matrix or data frame to be 
  #             analyzed. Its columns are the variables.
  # alpha     : In each column, at least n*alpha cells must
  #             remain unflagged. Defaults to 75%, should not
  #             be set (much) lower.
  # quant     : Determines the cutoff value to flag cells. 
  #             Defaults to 0.99
  # crit      : The iteration stops when successive covariance
  #             matrices (of the standardized data) differ by
  #             less than crit. Defaults to 1e-4.
  # noCits    : The maximal number of C-steps used.
  # lmin      : a lower bound on the eigenvalues of the 
  #             estimated covariance matrix on the 
  #             standardized data. Defaults to 1e-04.
  #             Should not be smaller than 1e-6.
  # lmax      : if not NULL, an upper bound on the eigenvalues 
  #             of the estimated covariance matrix on the 
  #             standardized data. Could e.g. be set to
  #             2*ncol(X) since the largest eigenvalue is at
  #             most d max_{ij} |C_{ij}| = d max_j |C_{jj}|.
  # fixedCenter: if TRUE, cellMCD is fit with a fixed center at O.
  # checkPars : Optional list of parameters used in the call 
  #            to cellMCD. The options are:
  #     - coreOnly: If TRUE, skip the execution of checkDataset. 
  #         Defaults to FALSE
  #     - numDiscrete:  A column that takes on numDiscrete or 
  #         fewer values will be considered discrete and not
  #         retained in the cleaned data. Defaults to 5.
  #     - precScale: Only consider columns whose scale is larger 
  #         than precScale. Here scale is measured by the median 
  #         absolute deviation. Defaults to 1e-12.
  #     - silent: Whether or not the function progress messages 
  #         should be printed. Defaults to FALSE.
  
  # First some auxiliary functions that are only used here, so 
  # they don't show up in the list of functions:
  
  
  # Here the main function starts:
  #
  X <- as.matrix(X) 
  if (!"coreOnly" %in% names(checkPars)) {
    checkPars$coreOnly <- FALSE
  }
  if (!"numDiscrete" %in% names(checkPars)) {
    checkPars$numDiscrete <- 5
  }
  if (!"precScale" %in% names(checkPars)) {
    checkPars$precScale <- 1e-12
  }
  if (!"silent" %in% names(checkPars)) {
    checkPars$silent <- FALSE
  }
  
  if (!"fracNA" %in% names(checkPars)) {
    checkPars$fracNA <- 0.5
  }
  
  if (!checkPars$coreOnly) {
    X <- checkDataSet(X, fracNA = 0.5, 
                      numDiscrete = checkPars$numDiscrete, 
                      precScale = checkPars$precScale, 
                      silent = checkPars$silent)$remX
  }
  # check dimension of the data w.r.t. the sample size
  if (nrow(X) < 5 * ncol(X)) {
    warning(paste0("There are fewer than 5 cases per dimension",
                   " in this data set.\n",
                   "It is not recommended to run cellMCD",
                   " on these data.\n",
                   "Consider reducing the number of variables."))
  }
  if (lmin < 1e-6) stop("lmin should be at least 1e-6.")
  
  #
  # Robustly standardize the data
  #
  if (fixedCenter) {
    locsca  <- cellWise::estLocScale(X, center = FALSE)
  } else {
    locsca  <- cellWise::estLocScale(X)
  }
  
  rscales <- locsca$scale
  Xs      <- scale(X, center = locsca$loc, scale = rscales)
  #
  # Check whether there are not too many bad cells:
  #
  margfrac <- colSums(abs(Xs) > sqrt(qchisq(0.99,1)), na.rm = TRUE) / nrow(Xs)
  if (max(margfrac) > 1 - alpha) {
    cat(paste0("\nAt least one variable of X has more than ",
               "100*(1-alpha)% = ", 100 * (1 - alpha), "%",
               "\nof marginal outliers.",
               "\nThe percentages per variable are:\n"))
    print(round(100 * margfrac, 2))
    stop("Too many marginal outliers.")
  }
  badfrac <- colMeans((abs(Xs) > sqrt(qchisq(0.99,1))) | (is.na(Xs))) # also includes real NAs in X
  if (max(badfrac) > 1 - alpha) {
    cat(paste0("\nAt least one variable of X has more than ",
               "100*(1-alpha)% = ", 100 * (1 - alpha), "%",
               "\nof marginal outliers plus NA's.",
               "\nThe percentages per variable are:\n"))
    print(round(100 * badfrac, 2))
    stop("Too many marginal outliers plus NA's.")
  }
  #
  DDCWout <- DDCWcov(Xs, maxCol = 1 - alpha, lmin = lmin, lmax = NULL,
                     fixedCenter = fixedCenter)
  initEst <- list(mu = DDCWout$center, Sigma = DDCWout$cov)
  #
  # We set the marginal outliers to NA so they stay flagged:
  Xs[abs(Xs) > 3] = NA
  #
  logC    <- -log(diag(solve(initEst$Sigma)))
  cutoff  <- qchisq(quant, df = 1)
  lambdas <- cutoff + logC + log(2 * pi)
  h       <- ceiling(alpha * nrow(Xs))
  temp    <- iterMCD_cpp(Xs,
                         initEst$mu,
                         initEst$Sigma,
                         h,
                         lambdas,
                         crit,
                         noCits,
                         lmin,
                         checkPars$precScale,
                         fixedCenter,
                         checkPars$silent)
  W <- temp$W
  rownames(W) <- rownames(X)
  colnames(W) <- colnames(X)
  out   <- allpreds_cpp(Xs, temp$Sigma, temp$mu, W)
  preds <- out$preds
  if (sum(is.na(as.vector(preds))) > 0 ) stop("There are missing preds!")
  cvars <- out$cvars
  if (sum(is.na(as.vector(cvars))) > 0 ) stop("There are missing cvars!") 
  if (min(as.vector(cvars)) <= 0) stop("There are cvars <= 0 !")
  rownames(preds) <- rownames(cvars) <- rownames(X)
  colnames(preds) <- colnames(cvars) <- colnames(X)
  
  if (!checkPars$silent) {
    percflag <- 100 * colMeans(1 - W, na.rm = TRUE)
    cat("Percentage of flagged cells per variable:\n")
    print(round(percflag,2))
  }
  mu    <- locsca$loc + temp$mu * rscales
  raw.S <- diag(rscales) %*% temp$Sigma %*% diag(rscales)
  S <- diag(rscales) %*% cov2cor(temp$Sigma) %*% diag(rscales)
  colnames(S) = rownames(S) = colnames(raw.S) = rownames(raw.S) = colnames(X)
  
  preds <- scale(preds, center = FALSE, scale = 1 / rscales)
  preds <- scale(preds, center = -locsca$loc, scale = FALSE)
  csds  <- scale(sqrt(cvars), center = FALSE, scale = 1 / rscales)
  Ximp  <- X
  Ximp[which(W == 0)] <- preds[which(W == 0)]
  Zres  <- (X - preds) / csds
  return(list(mu = mu, S = S,
              W = W, preds = preds,
              csds = csds, Ximp = Ximp,
              Zres = Zres, raw.S = raw.S,
              locsca = locsca,
              nosteps = temp$nosteps,
              X = X, quant = quant))
} 





plot_cellMCD = function(cellout, type = "Zres/X", whichvar = NULL,
                        horizvar = NULL, vertivar = NULL,  
                        hband = NULL, vband = NULL, drawellipse = T,
                        opacity = 0.5, identify = FALSE, 
                        ids=NULL, labelpoints = T, vlines = FALSE,
                        clines = TRUE, main = NULL,
                        xlab = NULL, ylab = NULL, xlim = NULL,
                        ylim = NULL, cex = 1, cex.main = 1.2, 
                        cex.txt = 0.8, cex.lab = 1, line=2.0){
  #
  # Function for making plots based on cellMCD output.
  #
  # Arguments:
  # cellout      output of function cellMCD()
  # type         "index", "Zres/X", "Zres/pred", "X/pred", or
  #              "bivariate".
  # whichvar     number or name of the variable to be plotted.  
  #              Not applicable when type == "bivariate".
  # horizvar     number or name of the variable to be plotted on the 
  #              horizontal axis. Only when type == "bivariate".
  # vertivar     number or name of the variable to be plotted on the 
  #              vertical axis. Only when type == "bivariate".
  # hband        draw a horizontal tolerance band? TRUE or FALSE. 
  #              NULL yields TRUE for types "index", "Zres/X", 
  #              and "Zres/pred".
  # vband        draw a vertical tolerance band? TRUE or FALSE.
  #              NULL yields TRUE for types "Zres/X", "Zres/pred", 
  #              and "X/pred".
  # drawellipse  whether to draw a 99% tolerance ellipse. Only
  #              for type == "bivariate".
  # opacity      opacity of the plotted points: 1 is fully opaque,
  #              less is more transparent.
  # identify     if TRUE, identify cases by mouseclick, then Esc.
  # ids          vector of case numbers to be emphasized in the plot.
  #              If NULL or of length zero, none are emphasized.
  # labelpoints  if TRUE, labels the points in ids by their
  #              row name in X.
  # vlines       for the points in ids, draw dashed vertical lines 
  #              from their standardized residual to 0 when type is
  #              "index", "Zres/X", or "Zres/pred". Draws dashed
  #              vertical ines to the diagonal for type "X/pred". 
  #              Can be TRUE or FALSE, default is FALSE.
  # clines       only for type == "bivariate". If TRUE, draws
  #              a red connecting line from each point in ids to 
  #              its imputed point, shown in blue.
  #
  # The following arguments are plot options, to finalize plots
  # for presentation:
  #
  # main         main title of the plot. If NULL, it is constructed
  #              automatically from the arguments.  
  # xlab         overriding label for x-axis, unless NULL.
  # ylab         overriding label for y-axis, unless NULL.
  # xlim         overriding limits of horizontal axis.
  # ylim         overriding limits of vertical axis.
  # cex          size of plotted points.
  # cex.main     size of the main title.
  # cex.lab      size of the axis labels.
  # cex.txt      size of the point labels.
  # line         distance of axis labels to their axis.
  #
  # Invisible output:
  # out      NULL, except when identify == TRUE. Then a list with:
  #          $ids    : the case number(s) that were identified
  #          $coords : coordinates of all points in the plot.
  
  # First some auxiliary functions:
  #
  identfy = function(xcoord,ycoord){
    # identify points in a plot(x,y)
    coordinates <- cbind(xcoord,ycoord)
    message("Press the escape key to stop identifying.")
    iout <- identify(coordinates, order = TRUE)
    ids <- iout$ind[order(iout$order)]
    if (length(ids) > 0) {
      cat("Identified point(s): ")
      print(ids)
    }   
    return(list(ids=ids, coords=coordinates[ids,,drop=F]))
  }
  
  addlabels = function(xcoord, ycoord, ids, labs, cex.txt = 0.8,
                       labtype = NULL){
    # Adds labels to plot. When labs = "i" it is the case number.
    # Also, labs can be the set of row names of the dataset.
    # For labs = "letters" we plot a, b, c, ...  
    #
    len = length(ids)
    if(len == 0) stop("ids has no elements")
    if(is.null(labtype)) { mylabs = labs[ids] 
    } else {
      if(labtype == "i") mylabs = ids
      if(labtype == "letters") mylabs = letters[seq_len(len)]
    }
    text(x = xcoord[ids], y = ycoord[ids], labels = mylabs, 
         cex = cex.txt)
  }
  
  vlines2diag = function(xcoord, ycoord, ids, lty=2, col="red"){
    # For indices in isd, plot vertical residual line 
    for(i in seq_len(length(ids))){ # i=2
      xys = c(xcoord[ids[i]],xcoord[ids[i]],
              ycoord[ids[i]],xcoord[ids[i]])
      lines(matrix(xys, ncol = 2, byrow=F), lty=lty, col=col)
    }
  }
  
  vlines2zero = function(xcoord, ycoord, ids, lty=2, col="red"){
    # For indices in isd, plot vertical residual line 
    for(i in seq_len(length(ids))){ # i=2
      xys = c(xcoord[ids[i]],xcoord[ids[i]],
              ycoord[ids[i]],0)
      lines(matrix(xys, ncol = 2, byrow=F), lty=lty, col=col)
    }
  }
  
  # Replaced ellipse::ellipse() by the simple function below,
  # so we do not have a dependency on library(ellipse).
  ellipsepoints = function(covmat, mu, quant=0.99, npoints = 100)
  { # computes points of the ellipse t(x-mu)%*%covmat%*%(x-mu) = c
    # with c = qchisq(quant,df=2)
    if (!all(dim(covmat) == c(2, 2))) stop("covmat is not 2 by 2")
    eig = eigen(covmat)
    U = eig$vectors
    R = U %*% diag(sqrt(eig$values)) %*% t(U) # square root of covmat
    angles = seq(0, 2*pi, length = npoints+1)
    xy = cbind(cos(angles),sin(angles)) # points on the unit circle
    fac = sqrt(qchisq(quant, df=2))
    scale(fac*xy%*%R, center = -mu, scale=FALSE)
  }  
  
  # Here the main function starts:
  #
  cutf = sqrt(qchisq(cellout$quant, df=1))
  mycol <- adjustcolor("black", alpha.f = opacity)
  redcol <- adjustcolor("red", alpha.f = 1) # = opacity)
  X = as.matrix(cellout$X)
  if(length(dim(X)) != 2) stop("cellout$X is not a matrix")
  ncol = ncol(X)
  nrow = nrow(X)
  cn = colnames(X)
  if(is.null(cn)) cn = seq_len(ncol)
  rn = rownames(X)
  if(is.null(rn)) rn = seq_len(nrow)
  if(type %in% c("index", "Zres/X", "Zres/pred", "X/pred")){
    if(is.null(whichvar)){
      stop("You must specify the variable whichvar to plot.")
    }
    if(whichvar %in% seq_len(ncol)){
      j = whichvar
      varlab = cn[j] 
    } else { 
      if(whichvar %in% cn){
        j = which(cn == whichvar)
        varlab = whichvar
      } else { stop(paste0("whichvar = ",whichvar," is not valid")) }
    }
    Xj     = X[,j]
    preds  = cellout$preds[,j]
    
    Zres = cellout$Zres[, j]
    flagged = which(abs(Zres) > cutf)
    #
    if(type == "index"){
      ycoord = Zres 
      xcoord = seq_len(length(ycoord))
      hlim = c(-cutf, cutf)
      obsp = which(!is.na(ycoord)) # points to plot
      if(is.null(ylim)) ylim = range(c(ycoord[obsp],hlim), na.rm=T)    
      plot(xcoord[obsp], ycoord[obsp], xlim=xlim, ylim=ylim, pch=16, 
           xlab="", ylab="", col=mycol, cex=cex) 
      if(is.null(xlab)) xlab = "index"
      title(xlab = xlab, line=line, cex.lab = cex.lab)
      if(is.null(ylab)) ylab = paste0("standardized residual of ",varlab)
      title(ylab = ylab, line=line, cex.lab = cex.lab)
      if(is.null(main)) main = 
        paste0("index plot: standardized residual of ",varlab)
      title(main = main, line = 1, cex.main = cex.main)
      if(is.null(hband) || hband==T){
        abline(h = hlim, lwd = 3, col="darkgray") }
      # in an index plot we cannot draw a meaningful vband.
      points(xcoord[flagged], ycoord[flagged], pch=16, col=redcol)
      if(length(ids) > 0){
        if(labelpoints) {
          addlabels(xcoord, ycoord, ids, labs = rn, cex.txt = cex.txt)
        }
        if(vlines == TRUE){
          vlines2zero(xcoord, ycoord, ids, lty=2, col="red")
        }
      }
    }
    #
    if(type == "Zres/X"){
      xcoord = Xj
      ycoord = Zres
      lsx  = cellWise::estLocScale(xcoord) 
      vlim = c(lsx$loc - cutf*lsx$scale, lsx$loc + cutf*lsx$scale)
      hlim = c(-cutf, cutf)
      obsp = which(!is.na(xcoord) & !is.na(ycoord)) # points to plot
      if(is.null(xlim)) xlim = range(c(xcoord[obsp],vlim), na.rm=T)
      if(is.null(ylim)) ylim = range(c(ycoord[obsp],hlim), na.rm=T)    
      plot(xcoord[obsp], ycoord[obsp], xlim=xlim, ylim=ylim, pch=16, 
           xlab="", ylab="", col=mycol, cex=cex) 
      if(is.null(xlab)) xlab=varlab
      title(xlab = xlab, line=line, cex.lab = cex.lab)
      if(is.null(ylab)) ylab = paste0("standardized residual of ",varlab)
      title(ylab = ylab, line=line, cex.lab = cex.lab)
      if(is.null(main)) main = paste0("standardized residual versus X for ",varlab)
      title(main = main,line = 1, cex.main = cex.main) 
      if(is.null(hband) || hband==T){
        abline(h = hlim, lwd = 3, col="darkgray") }  
      if(is.null(vband) || vband==T){
        abline(v = vlim, lwd = 3, col="darkgray") }
      points(xcoord[flagged], ycoord[flagged], pch=16, col=redcol)
      if(length(ids) > 0){
        if(labelpoints) {      
          addlabels(xcoord, ycoord, ids, labs = rn, cex.txt = cex.txt)
        }
        if(vlines == TRUE){
          abline(h=0) 
          vlines2zero(xcoord, ycoord, ids, lty=2, col="red") 
        }
      }
    }
    #
    if(type == "Zres/pred"){  
      xcoord = preds
      ycoord = Zres
      lsx  = cellWise::estLocScale(xcoord) 
      vlim = c(lsx$loc - cutf*lsx$scale, lsx$loc + cutf*lsx$scale)
      hlim = c(-cutf, cutf)
      obsp = which(!is.na(xcoord) & !is.na(ycoord)) # points to plot
      if(is.null(xlim)) xlim = range(c(xcoord[obsp],vlim), na.rm=T)
      if(is.null(ylim)) ylim = range(c(ycoord[obsp],hlim), na.rm=T)    
      plot(xcoord[obsp], ycoord[obsp], xlim=xlim, ylim=ylim, pch=16, 
           xlab="", ylab="", col=mycol, cex=cex) 
      if(is.null(xlab)) xlab = paste0("predicted ",varlab)
      title(xlab=xlab, line=line, cex.lab = cex.lab)
      if(is.null(ylab)) ylab = paste0("standardized residual of ",varlab)
      title(ylab = ylab, line=line, cex.lab = cex.lab)
      if(is.null(main)) main = paste0(
        "standardized residual versus prediction for ",varlab)
      title(main = main,line = 1, cex.main = cex.main)
      if(is.null(hband) || hband==T){
        abline(h = hlim, lwd = 3, col="darkgray") }  
      if(is.null(vband) || vband==T){
        abline(v = vlim, lwd = 3, col="darkgray") }
      points(xcoord[flagged], ycoord[flagged], pch=16, col=redcol)
      if(length(ids) > 0){
        if(labelpoints) {
          addlabels(xcoord, ycoord, ids, labs = rn, cex.txt = cex.txt)
        }
        if(vlines == TRUE){
          abline(h=0) 
          vlines2zero(xcoord, ycoord, ids, lty=2, col="red") 
        }
      }
    }
    #
    if(type == "X/pred"){
      xcoord = preds
      ycoord = Xj
      lsx  = cellWise::estLocScale(xcoord) 
      lsy  = cellWise::estLocScale(ycoord)
      vlim = c(lsx$loc - cutf*lsx$scale, lsx$loc + cutf*lsx$scale)
      hlim = c(lsy$loc - cutf*lsy$scale, lsy$loc + cutf*lsy$scale)
      obsp = which(!is.na(xcoord) & !is.na(ycoord)) # points to plot
      if(is.null(xlim)) xlim = range(c(xcoord[obsp],vlim), na.rm=T)
      if(is.null(ylim)) ylim = range(c(ycoord[obsp],hlim), na.rm=T)    
      plot(xcoord[obsp], ycoord[obsp], xlim=xlim, ylim=ylim, pch=16, 
           xlab="", ylab="", col=mycol, cex=cex) 
      if(is.null(xlab)) xlab = paste0("predicted ",varlab)
      title(xlab = xlab, line=line, cex.lab = cex.lab)
      if(is.null(ylab)) ylab = paste0("observed ",varlab)
      title(ylab = ylab, line=line, cex.lab = cex.lab)
      if(is.null(main)) main = paste0(varlab," versus its prediction")
      title(main = main,line = 1, cex.main = cex.main)
      abline(0,1)
      if(!is.null(hband) && hband == TRUE){
        abline(h = hlim, lwd = 3, col="darkgray") }
      if(!is.null(vband) && vband == TRUE){
        abline(v = vlim, lwd = 3, col="darkgray") }
      points(xcoord[flagged], ycoord[flagged], pch=16, col=redcol)
      if(length(ids) > 0){
        if(labelpoints) {
          addlabels(xcoord, ycoord, ids, labs = rn, cex.txt = cex.txt)
        }
        if(vlines == TRUE){
          vlines2diag(xcoord, ycoord, ids, lty=2, col="red")
        }   
      }
    }
  } else {
    if(type == "bivariate"){
      if(is.null(horizvar)) {
        stop(paste0("You must specify the variable horizvar\n",
                    "to plot on the horizontal axis.")) }
      if(is.null(vertivar)) {
        stop(paste0("You must specify the variable verticar\n",
                    "to plot on the vertical axis.")) }
      if(horizvar %in% seq_len(ncol)){
        jj = horizvar
        hlab = cn[jj] 
      } else { 
        if(horizvar %in% cn){
          jj = which(cn == horizvar)
          hlab = horizvar
        } else { stop(paste0("horizvar = ",horizvar," is not valid")) }
      }
      if(vertivar %in% seq_len(ncol)){
        kk = vertivar
        vlab = cn[kk] 
      } else { 
        if(vertivar %in% cn){
          kk = which(cn == vertivar)
          vlab = vertivar
        } else { stop(paste0("vertivar = ",vertivar," is not valid")) }
      }
      if(jj==kk) stop("horivar and vertivar should differ.")
      xcoord = X[,jj]
      ycoord = X[,kk]
      imp = cellout$Ximp
      ell = ellipsepoints(covmat = cellout$S[c(jj,kk),c(jj,kk)],
                          mu = cellout$mu[c(jj,kk)], quant=0.99,
                          npoints=500)
      lsx = cellWise::estLocScale(xcoord) 
      lsy = cellWise::estLocScale(ycoord)
      vlim = c(lsx$loc - cutf*lsx$scale, lsx$loc + cutf*lsx$scale)
      hlim = c(lsy$loc - cutf*lsy$scale, lsy$loc + cutf*lsy$scale)
      obsp = which(!is.na(xcoord) & !is.na(ycoord)) # points to plot
      if(is.null(xlim)) {
        xlim = range(c(xcoord[obsp],imp[obsp,jj],ell[,1],vlim), na.rm=T)
      }
      if(is.null(ylim)) {
        ylim = range(c(ycoord[obsp],imp[obsp,kk],ell[,2],hlim),na.rm=T)
      }
      flagged = which(cellout$W[,jj]*cellout$W[,kk] == 0) 
      plot(xcoord[obsp], ycoord[obsp], xlim=xlim, ylim=ylim, pch=16, 
           xlab="", ylab="", col=mycol, cex=cex) 
      if(is.null(xlab)) xlab = hlab
      title(xlab = xlab, line=line, cex.lab = cex.lab)
      if(is.null(ylab)) ylab = vlab
      title(ylab = ylab, line=line, cex.lab = cex.lab)
      if(is.null(main)) main = paste0(vlab," versus ",hlab)
      title(main = main,line = 1, cex.main = cex.main)
      if(is.null(hband)) hband = F
      if(hband == T) abline(h = hlim, lwd=3, col="darkgray")
      if(is.null(vband)) vband = F
      if(vband == T) abline(v = vlim, lwd=3, col="darkgray")
      if(drawellipse) lines(ell, lwd=3, col="darkgray")
      points(xcoord[flagged], ycoord[flagged], pch=16, col=redcol)
      if(length(ids) > 0){
        if(labelpoints) {
          addlabels(xcoord, ycoord, ids, labs = rn, cex.txt = cex.txt)
        }
        if(is.null(clines) || clines==T){
          imp = cellout$Ximp
          for(i in ids) {
            if(i %in% obsp){
              lines(x=c(X[i,jj],imp[i,jj]),y=c(X[i,kk],imp[i,kk]), 
                    lty=1, col="red")
              points(x=imp[i,jj],y=imp[i,kk],pch=16,col="blue")
            }
          }
        }   
      }
    } else  stop(paste0("type = \"",type,"\" is invalid"))
  }
  if(identify) invisible(identfy(xcoord,ycoord))
}





truncEig <- function(S, lmin = NULL, lmax = NULL) {
  # Truncates the eigenvalues of S to between lmin and lmax.
  nrS <- nrow(S)
  if (ncol(S) != nrS) stop(" S is not square")
  if (mean(as.vector(abs(S - t(S)))) > 1e-8) {
    stop(" S is not symmetric")
  }
  Sout <- S
  if (!(is.null(lmin) && is.null(lmax))) {
    if (!is.null(lmin)) {
      if (lmin < 0) stop(paste0("lmin = ",lmin," must be >= 0"))
      if (!(lmin < Inf)) stop(paste0("lmin = ",lmin," must be finite"))
    }
    if (!is.null(lmax)) {
      if (lmax < 0) stop(paste0("lmax = ",lmax," must be >= 0"))
      if (!(lmax < Inf)) stop(paste0("lmax = ",lmax," must be finite"))
    }
    eig <- eigen(0.5 * (S + t(S)))
    vals <- eig$values
    if (!is.null(lmin)) { # there is a lower bound
      vals <- pmax(vals, lmin) 
    }
    if (!is.null(lmax)) { # there is an upper bound
      vals <- pmin(vals, lmax) 
    }    
    Sout <- eig$vectors %*% diag(vals) %*% t(eig$vectors) 
  } 
  return(Sout)
}
