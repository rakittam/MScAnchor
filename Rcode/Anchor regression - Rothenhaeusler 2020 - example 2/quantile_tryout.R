# Quantile tryout with example 2 in Anchor Regression - Rothenhaeuser 2020
#
# Author: Maic Rakitta
# Date: 24.09.20
##########################################################################
set.seed(1)

n <- 1000 # number of samples from unpertubed and pertubed distribution

# initialize training data
library(extraDistr) # for rademacher distribution
A <- rsign(n)
epsH.train <- rnorm(n)
epsX.train <- rnorm(n)
epsY.train <- rnorm(n)

H.train <- epsH.train
X.train <- A+H.train+epsX.train
Y.train <- X.train+2*H.train+epsY.train

# initialize test data
epsH.test <- rnorm(n)
epsX.test <- rnorm(n)
epsY.test <- rnorm(n)

H.test <- epsH.test
X.test <- 1.8+H.test+epsX.test
Y.test <- X.test+2*H.test+epsY.test

##########################################################################
# Anchor regression
anchor.regression <- function(X, Y, A, gamma, n){
  
  P.A <- A%*%solve(t(A)%*%A)%*%t(A)
  W <- diag(n)-(1-sqrt(gamma))*P.A
  
  Y.tilde <- W%*%Y
  X.tilde <- W%*%X
  
  fit <- lm(Y.tilde~X.tilde-1)
}


##########################################################################
# Quantile usage for gamma calculation

alpha <- (1:100)/101
gamma.quantile <- qchisq(alpha, df=1) 

MSE.train.vec.quantile <- numeric(length(gamma.quantile))
MSE.test.vec.quantile <- numeric(length(gamma.quantile))
b.vec.quantile <- numeric(length(gamma.quantile))
for (i in 1:length(gamma.quantile)){
  
  gamma <- gamma.quantile[i]
  
  fit <- anchor.regression(X.train, Y.train, A, gamma, n)
  
  MSE.train.vec.quantile[i] <- mean((Y.train - t(X.train)*coef(fit)) ^ 2)
  MSE.test.vec.quantile[i] <- mean((Y.test - t(X.test)*coef(fit)) ^ 2)
  
  b.vec.quantile[i] <- coef(fit)
}

plot(gamma.quantile,MSE.train.vec.quantile)
plot(gamma.quantile,MSE.test.vec.quantile)

plot(alpha,MSE.train.vec.quantile, type = "l", ylim = c(min(min(MSE.train.vec.quantile),min(MSE.test.vec.quantile)),max(max(MSE.train.vec.quantile),max(MSE.test.vec.quantile))))
lines(alpha,MSE.test.vec.quantile, col=2)
legend(0.5, 6.5, legend=c("train", "test"),
       col=c(1,2), cex=0.8,lty=1)

quantile.MSE.train <- as.vector(quantile(MSE.train.vec.quantile, alpha))
quantile.MSE.test <- as.vector(quantile(MSE.test.vec.quantile, alpha))

plot(alpha,quantile.MSE.train, type = "l", ylim = c(min(min(quantile.MSE.train),min(quantile.MSE.train)),max(max(quantile.MSE.test),max(quantile.MSE.test))))
lines(alpha,quantile.MSE.test, col=2)
legend(0, 6.5, legend=c("train", "test"),
       col=c(1,2), cex=0.8,lty=1)


gamma.quantile[which(MSE.train.vec.quantile==min(MSE.train.vec.quantile))]
alpha.optimal <- alpha[which(MSE.train.vec.quantile==min(MSE.train.vec.quantile))]
alpha.optimal

qchisq(alpha.optimal, df=1) 
