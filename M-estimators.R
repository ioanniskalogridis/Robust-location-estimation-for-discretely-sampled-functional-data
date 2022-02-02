library(robustbase)
library(pcaPP)

# rho.huber <- function(x, k = 0.193){
#   rho.huber <- ifelse(abs(x) <= k, x^2/2, k*(abs(x)-k/2) )
#   return(rho.huber)
# }

# psi.huberc <- function(x, k = 0.193){
#   d <- ifelse( abs(x) <= k, x, k*sign(x))
#   return(d)
# }

#psi.huberc(2)
#psi.huber(2, k = 0.733)
#curve(psi.huberc, from = -5, to = 5 )
#curve(psi.huber, from = -5, to = 5)

# Eff.h <- function(x){
#   Eff.h <- (2*pnorm(x)-1)^2/(2*x^2*(1-pnorm(x)) + 2*pnorm(x) - 1 - 2*x*dnorm(x))
#   return(Eff.h)
# }
#Eff.h(x = 1.345)
# Eff.h(x = 0.733) #85% Efficiency
#Eff.h(x = 0.193) #70% Efficiency
#Eff.h(x = 0.04) #65% Efficiency

#Eff.h(1.55)
#BIC.rob <- function(x){   NOT USED
#  n <- length(x$fitted.values)
#  k <- length(x$coefficients)-1
#  BIC.rob <- -2 * sum(  rho.huber(residuals(x)/x$scale) + log(x$scale) ) - k*log(n)
#  return(-BIC.rob)
#}

M.est <- function(x, tol = 1e-10, max.it = 50, psi = "huber"){
  #Functiona M-estimator of location (Sinova et al. 2017)
  #Starting curve is the L1 median
  psi = psi
  x = as.matrix(x)
  n <- dim(x)[1]
  istop = 0
  ic = 0 
  gkm.in <- l1median(x)
  Distances.prelim <- t(apply(x, 1, FUN = function(y) ( y-gkm.in )^2 /dim(x)[2] )  )
  Distances.prelim <- apply(Distances.prelim, 1, FUN = function(x) sqrt(sum(x))  )
  alpha <- quantile(Distances.prelim, probs = 0.5, names = F)
  beta <- quantile(Distances.prelim, probs = 0.75, names = F)
  ce <- quantile(Distances.prelim, probs = 0.85, names = F)
  
  if(psi == "bisquare") {
    cc = beta
  } else if (psi == "huber") {
      cc = alpha
      } else 
        cc = c(alpha, beta, ce)
  
  while(istop == 0 && (ic <= max.it) ){
   ic = ic + 1
   Distances.o <- t(apply(x, 1, FUN = function(y) ( y-gkm.in )^2 /dim(x)[2] )  )
   Distances.o <- apply(Distances.o, 1, FUN = function(x) sqrt(sum(x))  )
   
   if(psi == "huber") {
     J.o = mean( Mpsi(Distances.o, cc = cc, psi = "huber", deriv = -1)  )
   } else {J.o = mean( Mpsi(Distances.o, cc = cc, psi = psi, deriv = -1)  )
   }
   u <- Mwgt(Distances.o, cc = cc, psi = psi)
   u <- u /sum(u) 
   
   gkm <- diag(u) %*%x
   gkm <- apply(gkm, 2, sum)
   Distances.n <- t(apply(x, 1, FUN = function(y) ( y-gkm )^2 /dim(x)[2] )  )
   Distances.n <- apply(Distances.n, 1, FUN = function(x) sqrt(sum(x))  )
   
   if(psi == "huber") {
     J.n = mean( rho.huber(Distances.n, k = cc)  )
   } else {J.n = mean( Mpsi(Distances.n, cc = cc, psi = psi, deriv = -1)  )
   }
   
   check <- abs(J.n- J.o)/J.o
   
   if(check < tol){istop = 1}
   gkm.in = gkm
  }
  list(gkm = gkm, check = check, ic = ic)
}

M.est.nd <- function(x, tol = 1e-10, max.it = 50, psi = "huber", k = 1.345){
  # Functiona M-estimator of location (Sinova et al. 2017)
  # Starting curve is the L1 median
  # k  is the tuning
  psi = psi
  x = as.matrix(x)
  n <- dim(x)[1]
  istop = 0
  ic = 0 
  gkm.in <- l1median(x)
  # Distances.prelim <- t(apply(x, 1, FUN = function(y) ( y-gkm.in )^2 /dim(x)[2] )  )
  # Distances.prelim <- apply(Distances.prelim, 1, FUN = function(x) sqrt(sum(x))  )
  # alpha <- quantile(Distances.prelim, probs = 0.5, names = F)
  # beta <- quantile(Distances.prelim, probs = 0.75, names = F)
  # ce <- quantile(Distances.prelim, probs = 0.85, names = F)
  
  cc = k
  
  while(istop == 0 && (ic <= max.it) ){
    ic = ic + 1
    Distances.o <- t(apply(x, 1, FUN = function(y) ( y-gkm.in )^2 /dim(x)[2] )  )
    Distances.o <- apply(Distances.o, 1, FUN = function(x) sqrt(sum(x))  )
    
    if(psi == "huber") {
      J.o = mean( Mpsi(Distances.o, cc = cc, psi = "huber", deriv = -1)  )
    } else {J.o = mean( Mpsi(Distances.o, cc = cc, psi = psi, deriv = -1)  )
    }
    u <- Mwgt(Distances.o, cc = cc, psi = psi)
    u <- u /sum(u) 
    
    gkm <- diag(u) %*%x
    gkm <- apply(gkm, 2, sum)
    Distances.n <- t(apply(x, 1, FUN = function(y) ( y-gkm )^2 /dim(x)[2] )  )
    Distances.n <- apply(Distances.n, 1, FUN = function(x) sqrt(sum(x))  )
    
    if(psi == "huber") {
      J.n = mean( Mpsi(Distances.n, cc = cc, psi = "huber", deriv = -1)  )
    } else {J.n = mean( Mpsi(Distances.n, cc = k, psi = psi, deriv = -1)  )
    }
    
    check <- abs(J.n- J.o)/J.o
    
    if(check < tol){istop = 1}
    gkm.in = gkm
  }
  list(gkm = gkm, check = check, ic = ic)
}


#x <- matrix(rnorm(5000, 200, 100), 500, 100)
#M.est(x)
# 
# radar <- read.table(file = file.choose(), header = F)
# x = as.matrix(radar)
# dim(radar)
# M.est(radar)
# x = radar
# dim(x)
# M.est(radar)
# par(mfrow = c(1, 1))
# matplot(t(radar), col = "gray", type = "l")
# plot(1:70, M.est(radar, psi = "bisquare")$gkm, type = "l")
# lines(1:70, M.est(radar, psi = "hampel")$gkm, type = "l", lwd = 2)
# lines(1:70, M.est(radar, psi = "huber")$gkm, type = "l", lwd = 2)
# 
# fit.v <- M.pen.sp.r(grid = 1:70, Y=radar, n.basis = 40, par.in = 1)
# lines(fit.v$mu, col = "blue", lwd = 3)



# s <- create.bspline.basis(c(1, 70), nbasis = 100)
# M.dat.h <- Data2fd(1:70, y = as.vector(M.est(radar, psi = "hampel")$gkm), basisobj = s)
# plot(M.dat.h)
# M.dat.b <- Data2fd(1:70, y = as.vector(M.est(radar, psi = "bisquare")$gkm), basisobj = s)
# plot(M.dat.b, add = TRUE)
# M.dat.b <- Data2fd(1:70, y = as.vector(M.est(radar, psi = "huber")$gkm), basisobj = s)
# plot(M.dat.b, add = TRUE)

#Eff.bisq <- function(c){
#  Num <- integrate(function(x) (Mpsi(x, cc = c, psi = "bisquare"))^2 *dnorm(x), -Inf, Inf )$value
#  Denom <- (integrate( function(x) Mpsi(x, cc = c, psi = "bisquare", deriv = 1)*dnorm(x), -Inf, Inf )$value)^2
#  return(Denom/Num)
#}

#Eff.bisq(2.6973)#70% efficiency
#Eff.bisq(2.524)#65% efficiency
#Eff.bisq(2.8978) #75% efficiency
#Eff.bisq(3.446)
#Eff.bisq(2.367) #60% efficiency
