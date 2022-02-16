require(robustbase)
require(pcaPP)

M.est <- function(x, tol = 1e-10, max.it = 200, k = 1.345){
  # Functional Huber M-estimator of location (Sinova et al. 2018)
  # Starting curve is the L1 median
  # k  is the tuning parameter for the Huber loss function
  
  x = as.matrix(x)
  n <- dim(x)[1]
  istop = 0
  ic = 0 
  gkm.in <- l1median(x)
  
  cc = k
  
  while(istop == 0 && (ic <= max.it) ){
    ic = ic + 1
    Distances.o <- t(apply(x, 1, FUN = function(y) ( y-gkm.in )^2 /dim(x)[2] )  )
    Distances.o <- apply(Distances.o, 1, FUN = function(x) sqrt(sum(x))  )
    
    J.o = mean( Mpsi(Distances.o, cc = cc, psi = "huber", deriv = -1)  )
    u <- Mwgt(Distances.o, cc = cc, psi = psi)
    u <- u /sum(u) 
    
    gkm <- diag(u) %*%x
    gkm <- apply(gkm, 2, sum)
    Distances.n <- t(apply(x, 1, FUN = function(y) ( y-gkm )^2 /dim(x)[2] )  )
    Distances.n <- apply(Distances.n, 1, FUN = function(x) sqrt(sum(x))  )
    
    J.n = mean( Mpsi(Distances.n, cc = cc, psi = "huber", deriv = -1)  )
    
    
    check <- abs(J.n- J.o)/J.o
    
    if(check < tol){istop = 1}
    gkm.in = gkm
  }
  list(gkm = gkm, check = check, ic = ic)
}
