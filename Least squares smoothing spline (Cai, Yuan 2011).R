require(fda)
require(SparseM)

ls.smsp <-  function(Y, interval = NULL){
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  grid <- 1:p/p
  
  b.basis <- create.bspline.basis(rangeval = c(0,1), breaks = grid)
  b.basis.e <- eval.basis(b.basis, grid)
  
  Pred.big <- b.basis.e
  Pen.matrix <- bsplinepen(b.basis)
  
  Y.v = as.vector(t(Y))
  X.v = matrix(rep(t(Pred.big), dim(Y)[1]), ncol = ncol(Pred.big), byrow = TRUE)
  
  GCV <- function(lambda){
    M1 <-  t(X.v)%*%X.v + lambda * Pen.matrix
    fit.ls <- as.vector(SparseM::solve( M1, t(X.v)%*%Y.v ))
    resids.ls <- as.vector(Y.v - X.v%*%fit.ls)
    hat.values <- mean(diag(X.v%*%SparseM::solve( M1, t(X.v))))
    # hat.tr <- mean(t(X.v)*SparseM::solve(M1, t(X.v) ) )*(p+2)
    GCV.scores <- mean( (resids.ls)^2  )/( (1-hat.values)^2 )
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
  
  fit.lsf <- SparseM::solve( t(X.v)%*%X.v + lambda1 * Pen.matrix, t(X.v)%*%Y.v )
  
  mu <- as.vector(Pred.big%*%fit.lsf)
  
  return(list(mu = mu, Pen.matrix = Pen.matrix,lambda = lambda1))
}
