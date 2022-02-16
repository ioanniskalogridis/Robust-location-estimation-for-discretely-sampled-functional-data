# Huber-type smoothing spline for location estimation from discretely sampled functional data.
# This function implements the estimator of Kalogridis and Van Aelst (2021), currently without any scale.
# The implementation is with a penalized version of IRLS described in (Maronna 2011).
# Selection of the penalty parameter is with weighted GCV.

library(fda)
library(SparseM)
library(MASS)

huber.smsp <- function(Y, k = 0.70, maxit = 1000, toler = 1e-12, r = 2, interval = NULL){
  # Y is the matrix containing the data values, each row corresponds to an observation.
  # Missing values, i.e., irregularly sampled data, are allowed.
  # k is the tuning parameter that regulates the behaviour of the Huber loss.
  # maxit refers to the maximum number of iterations allowed for the irls algorithm.
  # r is the order of the penalty; the order of the spline is 2*r.
  # interval is the interval in which the algorithm will look for the optimal value of lambda,
  # i.e., the value of lambda that minimizes the weighted GCV.
  # if interval = NULL, a rough grid search is performed followed by optimization in the most promising interval.

  rho.huber <- function(x, k = 1.345) ifelse(abs(x)<= k, x^2/2, k*abs(x)-k^2/2)
  psi.huber <- function(x, k = 1.345) ifelse(abs(x)<= k, x, k*sign(x))
  weights.huber <- function(x, k = 1.345) ifelse( x==0, 1, psi.huber(x, k = k)/x)
  
  Y=Y[rowSums(is.na(Y)) !=ncol(Y), ]
  Y = as.matrix(Y)
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  grid <- 1:p/p
  b.basis <- create.bspline.basis(rangeval = c(grid[1], grid[p]), breaks = grid, norder = 2*r)
  b.basis.e <- eval.basis(b.basis, grid)
  T.m <- t(apply(Y, 1, FUN = function(x) x <- 1:dim(Y)[2] ))
  T.m[is.na(Y)] <- NA
  B <- matrix(0, nrow = (p+2*r-2), ncol = (n*p)-sum(is.na(T.m))  )
  # B2 <- matrix(0, nrow = (p+2), ncol = (n*p)-sum(is.na(T.m))  )
  
  for(j in 1:((n*p)-sum(is.na(T.m)))){
    B[, j] <- b.basis.e[ na.omit(as.vector(t(T.m)))[j], ]
  }
  ms <- as.vector(apply(T.m, 1, FUN = function(x) p-sum(is.na(x))))
  h.m <- 1/mean(1/ms)
  ms <- rep(ms, times = c(ms))
  B.s <- scale(B, center = FALSE, scale =  ms)
  B.p <- B.s%*%t(B)
  
  # Penalty matrix, see, e.g., de Boor (2001)
  grid.a <- c(0, grid)
  if(r==1){
    Pen.matrix =  (diag(1/ (diff(grid.a)))%*%diff(diag(length(grid.a)), differences = 1) )%*%t(diag(1/ (diff(grid.a)))%*%diff(diag(length(grid.a)), differences = 1) )
  } else{
    Pen.matrix =  bsplinepen(b.basis, Lfdobj= r)
  }

  # Initial values, least-squares estimator with "optimal" penalty parameter
  par.in = 100*(h.m*n)^{-2*r/(2*r+1)}
  fit.in.c <- ginv(B.p+par.in*Pen.matrix)%*%(B.s%*%na.omit(as.vector(as.matrix(t(Y)))))
  resids.in <- na.omit(as.vector(t(Y)))-t(B)%*%fit.in.c
  
  # Iteratively reweighted least-squares
  huber.irls <- function(X, X.s, y, k = 1.345, tol = toler, lambda, Pen.matrix, resids.in, maxit = 20){
    
    ic = 0
    istop = 0
    
    while(istop == 0 & ic <= maxit){
      ic = ic + 1
      weights.prelim <- as.vector(weights.huber(resids.in, k = k))
      M1 <-  t(t(X.s)*weights.prelim)%*%t(X) + 2*lambda*Pen.matrix
      M2 <-  t(t(X.s)*weights.prelim)%*%na.omit(as.vector(as.matrix(t(y))))
      v1 = SparseM::solve(M1, M2)
      resids1 <- as.vector(na.omit(as.vector(t(y)))-t(X)%*%v1)
      check = max( abs(resids1-resids.in) ) 
      if(check < tol){istop =1}
      resids.in <- resids1
    }
    weights1 = as.vector(weights.huber(resids1,  k = k) )
    # fast trace of the hat-matrix
    hat.tr <- sum( diag(SparseM::solve(M1,  t(t(X.s)*weights1))%*%t(X)))/length(na.omit(as.vector(as.matrix(t(y)))))
    return(list(resids = resids1, beta.hat = v1, hat.tr = hat.tr, ic = ic,
                weights = weights1 ) ) 
  }
  
  GCV <- function(lambda){
    fit.r <- huber.irls(X = B, X.s = B.s, y = Y, resids.in = resids.in, 
                           lambda = lambda, Pen.matrix = Pen.matrix, k = k, maxit = maxit )
    GCV.scores <- mean( fit.r$weights*1/ms*(fit.r$resids)^2/((1-fit.r$hat.tr)^2)  )
    return(GCV.scores)
  }
  if(is.null(interval)){
    lambda.cand <- c(1e-09, 1e-08, 3e-08, 6e-08, 9e-08, 1e-07, 3e-07, 6e-07, 9e-07, 1e-06, 3e-06, 6e-06, 9e-06, 
                     1e-05, 3e-05, 6e-05, 9e-05, 1e-04, 3e-04, 6e-04, 9e-04,  1e-03, 4e-03, 7e-03,
                     1e-02, 4e-02, 7e-02, 1e-01, 6e-01, 2)
    lambda.e <- sapply(lambda.cand, FUN  = GCV)
    wm <- which.min(lambda.e)
    if(wm == 1){wm <- 2}
    if(wm == length(lambda.cand)){wm <- (length(lambda.cand)-1)  }
    lambda1 <- optimize(f = GCV, lower = lambda.cand[wm-1], upper = lambda.cand[wm+1])$minimum
  } else {
    lambda1 <- optimize(f = GCV, interval = interval)$minimum}
  
  fit.f <- huber.irls(X = B, X.s = B.s, y=  Y, resids.in = resids.in, 
                         lambda = lambda1, Pen.matrix = Pen.matrix, k = k, maxit = maxit)
  mu = b.basis.e%*%fit.f$beta.hat
  
  # plot(t, mu.t, type = "l", col = "black", lwd = 3)
  # lines(t, mu, col = "blue", lwd = 3)
  
  return(list(mu = mu, weights = fit.f$weights, Pen.matrix = Pen.matrix,
              lambda = lambda1))
}

