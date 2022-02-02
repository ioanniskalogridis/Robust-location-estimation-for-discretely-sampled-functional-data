require(fda)
require(sn)
require(EnvStats)

m <- 50
n <- 60
sigma <- 0.2
nrep <- 50

t <- seq(1/m, (m-1)/m, len = m)
# mu <- function(x) sin(6*pi*x)*(x+1)
# mu.t <- sapply(t, FUN = mu)
mu <- function(x) 3*exp(-(x-0.25)^2/0.3)
mu.t <- sapply(t, FUN = mu)

mse.hsp <- rep(NA, nrep)
mse.hf <- rep(NA, nrep)
mse.rs <- rep(NA, nrep)
mse.lp <- rep(NA, nrep)
mse.lssp <- rep(NA, nrep)
mse.l1sp <- rep(NA, nrep)

mse.hsp.m <- matrix(NA, nrow = m, ncol = nrep)
mse.hf.m <- matrix(NA, nrow = m, ncol = nrep)
mse.cao.m <- matrix(NA, nrow = m, ncol = nrep)
mse.lp.m <- matrix(NA, nrow = m, ncol = nrep)
mse.lssp.m <- matrix(NA, nrow = m, ncol = nrep)
mse.l1sp.m <- matrix(NA, nrow = m, ncol = nrep)

for(k in 1:nrep){
  print(k)
  X <- matrix(NA, nrow = n, ncol = m)
  for(i in 1:n){
    X[i,] <- mu.t 
    for(j in 1:50){
      X[i,] <- X[i, ] +  sqrt(2)*rt(1, df = 5)*sapply(t, FUN = function(x) sin((j-1/2)*pi*x)/((j-1/2)*pi) )
    }
  }
  # Y <- X + matrix( sigma*rnorm(m*n), nrow = n, ncol = m )
  # Y <- X + matrix( sigma*rt(m*n, df = 3), nrow = n, ncol = m )
  # matplot(t(Y), lwd = 3, lty = 1, col = "gray", type = "p", pch = 20)
  # Y <- X + matrix( sigma*rt(m*n, df = 3), nrow = n, ncol = m )
  # Y <- X + matrix( sigma*rst(m*n, nu = 3, alpha = 0.5), nrow = n, ncol = m )
  # Y <- X + matrix( sigma*rnormMix(m*n, mean1 = 0, sd1 = 1, mean2 = 0, sd2 = 9, p.mix = 0.15), nrow = n, ncol = m )
  Y <- X + matrix( sigma*rnorm(m*n)/runif(m*n), nrow = n, ncol = m )
  
  library(reshape2)
  library(ggplot2)
  data <- data.frame(t(Y))
  data$id <- 1:nrow(data)/nrow(data)
  #reshape to long format
  plot_data <- melt(data, id.var="id" )
  #plot
  p <- ggplot(plot_data, aes(x=id, y=value, group=variable, colour=variable)) + geom_point(col = "gray", size = 3.2)
  p <- p + theme_bw(base_size = 35) + labs(x = "t", y = "")
  p <- p + stat_function(aes(x = id, y = value), fun = mu, size = 1.4, colour = "black") + ylim(-7, 11)
  p
  
  fit.l1 <- quan.smsp(Y, interval = c(1e-06, 1e-05), alpha = 0.5)
  fit.h <- huber.smsp(Y, interval = c(9e-07, 9e-06))
  
  fit.sin <- M.est.nd(Y, k = 0.70)
  fit.cao <- Cao(Y=Y, k = 0.70)
  fit.deg <- Degr(Y = Y)
  fit.ls <- ls.smsp(Y= Y, interval = c(7e-05, 7e-04))
  
  mse.hsmsp[k] <- mean((fit.h$mu-mu.t)^2)
  mse.sin[k] <-  mean((fit.sin$gkm-mu.t)^2)
  mse.cao[k] <- mean((fit.cao$mu-mu.t)^2)
  mse.deg[k] <- mean((fit.deg$mu-mu.t)^2)
  mse.ls[k] <-  mean((fit.ls$mu-mu.t)^2)
  mse.l1sp[k] <- mean((fit.l1$mu-mu.t)^2)
  
  mse.hsmsp.m[, k] <- fit.h$mu
  mse.sin.m[, k] <- fit.sin$gkm
  mse.cao.m[, k] <- fit.cao$mu
  mse.deg.m[, k] <- fit.deg$mu
  mse.ls.m[, k] <- fit.ls$mu
  mse.l1sp.m[, k] <- fit.l1$mu
  
}
mean(mse.ls) ;  sd(mse.ls)/sqrt(nrep)
mean(mse.hsmsp) ; sd(mse.hsmsp)/sqrt(nrep)
mean(mse.deg) ; sd(mse.deg)/sqrt(nrep)
mean(mse.sin) ; sd(mse.sin)/sqrt(nrep)
mean(mse.cao) ; sd(mse.cao)/sqrt(nrep)
mean(mse.l1sp) ; sd(mse.l1sp)/sqrt(nrep)

matplot(mse.hsmsp.m, lty = 1, lwd = 3, col = "gray", type = "l", cex.lab = 1.3, cex.axis = 1.3) ; grid()
matplot(mse.sin.m, lty = 1, lwd = 3, col = "gray", type = "l", cex.lab = 1.3, cex.axis = 1.3) ; grid()
matplot(mse.ls.m, lty = 1, lwd = 3, col = "gray", type = "l", cex.lab = 1.3, cex.axis = 1.3) ; grid()
matplot(mse.l1sp.m, lty = 1, lwd = 3, col = "gray", type = "l", cex.lab = 1.3, cex.axis = 1.3) ; grid()

m <- 50
n <- 60
sigma <- 1
nrep <- 1000

t <- seq(1/m, (m-1)/m, len = m)
# mu <- function(x) sin(6*pi*x)*(x+1)
# mu.t <- sapply(t, FUN = mu)
mu <- function(x) 3*exp(-(x-0.25)^2/0.1)
mu.t <- sapply(t, FUN = mu)

mse.hsmsp <- rep(NA, nrep)
mse.sin <- rep(NA, nrep)
mse.cao <- rep(NA, nrep)
mse.deg <- rep(NA, nrep)
mse.ls <- rep(NA, nrep)

mse.hsmsp.m <- matrix(NA, nrow = m, ncol = nrep)
mse.sin.m <- matrix(NA, nrow = m, ncol = nrep)
mse.cao.m <- matrix(NA, nrow = m, ncol = nrep)
mse.deg.m <- matrix(NA, nrow = m, ncol = nrep)
mse.ls.m <- matrix(NA, nrow = m, ncol = nrep)


for(k in 1:nrep){
  print(k)
  X <- matrix(NA, nrow = n, ncol = m)
  for(i in 1:n){
    X[i,] <- mu.t 
    for(j in 1:50){
      X[i,] <- X[i, ] +  sqrt(2)*rt(1, df = 5)*sapply(t, FUN = function(x) sin((j-1/2)*p*x)/((j-1/2)*pi) )
    }
  }
  Y <- X + matrix( sigma*rnorm(m*n), nrow = n, ncol = m )
  # Y <- X + matrix( sigma*rt(m*n, df = 3), nrow = n, ncol = m )
  # Y <- X + matrix( sigma*rst(m*n, nu = 3, alpha = 0.5), nrow = n, ncol = m )
  # Y <- X + matrix( sigma*rnormMix(m*n, mean1 = 0, sd1 = 1, mean2 = 0, sd2 = 9, p.mix = 0.15), nrow = n, ncol = m )
  # Y <- X + matrix( sigma*rnorm(m*n)/runif(m*n), nrow = n, ncol = m )
  
  fit.h <- huber.smsp.cd(Y, interval = c(1.5e-03, 1.5e-02))
  fit.sin <- M.est.nd(Y, k = 0.70)
  fit.cao <- Cao(Y=Y, k = 0.70)
  fit.deg <- Degr(Y = Y)
  fit.ls <- ls.smsp(Y= Y, interval = c(5e-03, 5e-02))
  
  mse.hsmsp[k] <- mean((fit.h$mu-mu.t)^2)
  mse.sin[k] <-  mean((fit.sin$gkm-mu.t)^2)
  mse.cao[k] <- mean((fit.cao$mu-mu.t)^2)
  mse.deg[k] <- mean((fit.deg$mu-mu.t)^2)
  mse.ls[k] <-  mean((fit.ls$mu-mu.t)^2)
  
  mse.hsmsp.m[, k] <- fit.h$mu
  mse.sin.m[, k] <- fit.sin$gkm
  mse.cao.m[, k] <- fit.cao$mu
  mse.deg.m[, k] <- fit.deg$mu
  mse.ls.m[, k] <- fit.ls$mu
  
}
mean(mse.ls) ;  sd(mse.ls)/sqrt(nrep)
mean(mse.hsmsp) ; sd(mse.hsmsp)/sqrt(nrep)
mean(mse.deg) ; sd(mse.deg)/sqrt(nrep)
mean(mse.sin) ; sd(mse.sin)/sqrt(nrep)
mean(mse.cao) ; sd(mse.cao)/sqrt(nrep)

matplot(mse.hsmsp.m, lty = 1, lwd = 3, col = "gray", type = "l", cex.lab = 1.3, cex.axis = 1.3) ; grid()
matplot(mse.sin.m, lty = 1, lwd = 3, col = "gray", type = "l", cex.lab = 1.3, cex.axis = 1.3) ; grid()
matplot(mse.ls.m, lty = 1, lwd = 3, col = "gray", type = "l", cex.lab = 1.3, cex.axis = 1.3) ; grid()
matplot(mse.cao.m, lty = 1, lwd = 3, col = "gray", type = "l", cex.lab = 1.3, cex.axis = 1.3) ; grid()


############################################## Peak and partial contamination ###################################################


require(fda)
require(sn)
require(EnvStats)

nrep <- 1000
m <- 50
n <- 60
sigma <- 0.5

t <- seq(1/m, (m-1)/m, len = m)
mu <- function(x) sin(6*pi*x)*(x+1)
# mu.t <- sapply(t, FUN = mu)
# mu <- function(x) 3*exp(-(x-0.25)^2/0.3)
mu.t <- sapply(t, FUN = mu)

mse.hsp <- rep(NA, nrep)
mse.hf <- rep(NA, nrep)
mse.cao <- rep(NA, nrep)
mse.lp <- rep(NA, nrep)
mse.ls <- rep(NA, nrep)
mse.l1sp <- rep(NA, nrep)

mse.hsp.m <- matrix(NA, nrow = m, ncol = nrep)
mse.hf.m <- matrix(NA, nrow = m, ncol = nrep)
mse.cao.m <- matrix(NA, nrow = m, ncol = nrep)
mse.lp.m <- matrix(NA, nrow = m, ncol = nrep)
mse.ls.m <- matrix(NA, nrow = m, ncol = nrep)
mse.l1sp.m <- matrix(NA, nrow = m, ncol = nrep)

for(k in 1:nrep){
  print(k)
  X <- matrix(NA, nrow = n, ncol = m)
  for(i in 1:n){
    X[i,] <- mu.t 
    for(j in 1:50){
      X[i,] <- X[i, ] +  sqrt(2)*rt(1, df = 5)*sapply(t, FUN = function(x) sin((j-1/2)*pi*x)/((j-1/2)*pi) )
    }
    T.r <- runif(1, 1, m)
    sigma.r <- sample(c(-1, 1), size = 1, replace=  TRUE, prob = c(1/2, 1/2))
    epsilon.r <- rbinom(1, 1, prob = 0.1)
    
    X[i, seq_along(1:m)>T.r] <- X[i, seq_along(1:m)>T.r] + epsilon.r*sigma.r*5
  }
  # Y <- X + matrix( sigma*rnorm(m*n), nrow = n, ncol = m )
  # Y <- X + matrix( sigma*rt(m*n, df = 5), nrow = n, ncol = m )
  matplot(t(Y), lwd = 3, lty = 1, col = "gray", type = "l")
  Y <- X + matrix( sigma*rt(m*n, df = 3), nrow = n, ncol = m )
  # Y <- X + matrix( sigma*rst(m*n, nu = 3, alpha = 0.5), nrow = n, ncol = m )
  # Y <- X + matrix( sigma*rnormMix(m*n, mean1 = 0, sd1 = 1, mean2 = 0, sd2 = 9, p.mix = 0.15), nrow = n, ncol = m )
  # Y <- X + matrix( sigma*rnorm(m*n)/runif(m*n), nrow = n, ncol = m )
  
  # library(reshape2)
  # library(ggplot2)
  # data <- data.frame(t(Y))
  # data$id <- 1:nrow(data)/nrow(data)
  # #reshape to long format
  # plot_data <- melt(data, id.var="id" )
  # #plot
  # p <- ggplot(plot_data, aes(x=id, y=value, group=variable, colour=variable)) + geom_point(col = "gray", size = 2.5)
  # p <- p + theme_bw(base_size = 35) + labs(x = "t", y = "")
  # p <- p + stat_function(aes(x = id, y = value), fun = mu, size = 1.4, colour = "black") 
  # p
  

  fit.h <- huber.smsp(Y)
  fit.sin <- M.est.nd(Y, k = 0.70)
  # plot(t, mu.t, lwd = 3, type = "l")
  # lines(t,fit.h$mu, lwd = 3, col = "blue")
  # lines(t, fit.sin$gkm, lwd = 3, col = "red")
  # 
  # mean((fit.sin$gkm-mu.t)^2)*100
  # mean((fit.h$mu-mu.t)^2)*100
  # 
  fit.l1 <- quan.smsp(Y, alpha = 0.5)
  fit.cao <- Cao(Y=Y, k = 0.70)
  fit.deg <- Degr(Y = Y)
  fit.ls <- ls.smsp(Y= Y)
  
  mse.hsp[k] <- mean((fit.h$mu-mu.t)^2)
  mse.hf[k] <-  mean((fit.sin$gkm-mu.t)^2)
  mse.cao[k] <- mean((fit.cao$mu-mu.t)^2)
  mse.lp[k] <- mean((fit.deg$mu-mu.t)^2)
  mse.ls[k] <-  mean((fit.ls$mu-mu.t)^2)
  mse.l1sp[k] <- mean((fit.l1$mu-mu.t)^2)
  
  mse.hsp.m[, k] <- fit.h$mu
  mse.hf.m[, k] <- fit.sin$gkm
  mse.cao.m[, k] <- fit.cao$mu
  mse.lp.m[, k] <- fit.deg$mu
  mse.ls.m[, k] <- fit.ls$mu
  mse.l1sp.m[, k] <- fit.l1$mu
  
}
  mean(mse.ls, na.rm = TRUE)*100 ;  sd(mse.ls)*10/sqrt(nrep)
mean(mse.hsp, na.rm = TRUE)*100 ; sd(mse.hsp)*10/sqrt(nrep)
mean(mse.lp, na.rm = TRUE)*100 ; sd(mse.lp)*10/sqrt(nrep)
mean(mse.hf, na.rm = TRUE)*100 ; sd(mse.hf)*10/sqrt(nrep)
mean(mse.cao, na.rm = TRUE)*100 ; sd(mse.cao)*10/sqrt(nrep)
mean(mse.l1sp, na.rm = TRUE)*100 ; sd(mse.l1sp)*10/sqrt(nrep)

matplot(mse.hsp.m, lty = 1, lwd = 3, col = "gray", type = "l", cex.lab = 1.3, cex.axis = 1.3) ; grid()
matplot(mse.hf.m, lty = 1, lwd = 3, col = "gray", type = "l", cex.lab = 1.3, cex.axis = 1.3) ; grid()
matplot(mse.ls.m, lty = 1, lwd = 3, col = "gray", type = "l", cex.lab = 1.3, cex.axis = 1.3) ; grid()
matplot(mse.l1sp.m, lty = 1, lwd = 3, col = "gray", type = "l", cex.lab = 1.3, cex.axis = 1.3) ; grid()

