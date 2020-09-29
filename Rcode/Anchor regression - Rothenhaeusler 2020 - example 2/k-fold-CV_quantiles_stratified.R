# Optimal gamma calculations using k-fold cross validation
# in example 2 in Anchor Regression - Rothenhaeuser 2020
#
# Author: Maic Rakitta
# Date: 23.09.20
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


##########################################################################
# Anchor regression
anchor.regression <- function(X, Y, A, gamma, n){
  
  P.A <- A%*%solve(t(A)%*%A)%*%t(A)
  W <- diag(n)-(1-sqrt(gamma))*P.A
  
  Y.tilde <- W%*%Y
  X.tilde <- W%*%X
  
  fit <- lm(Y.tilde~X.tilde-1)
}

# PA = AR for gamma = 0
fit.PA <- anchor.regression(X.train, Y.train, A, 0, n)
b.PA <- coef(fit.PA)

# OLS = AR for gamma = 1
fit.OLS <- anchor.regression(X.train, Y.train, A, 1, n)
b.OLS <- coef(fit.OLS)

# IV - 2SLS
fit.2SLS.step <- lm(X.train~A-1)
X.train.hat <- fitted.values(fit.2SLS.step)
fit.IV <- lm(Y.train~X.train.hat-1)
b.IV <- coef(fit.IV)

##########################################################################
# k-fold Cross Validation
data <- data.frame(Y=Y.train,X=X.train, A=A) # create data frame for CV

# initialize
k <- 2
alpha.vec <- (1:100)/101

gamma.optimal <- numeric(length(alpha.vec))
quantile.CV.matrix <- matrix(nrow = length(alpha.vec), ncol = k)
for (a in 1:length(alpha.vec)) {
  
  alpha <- alpha.vec[a] # step 1: choose alpha
  gamma.CV <- qchisq(alpha, df=1)
  
  folds <- c(-1,1) # step 2: create folds
  
  # step 3: for varying gamma train and test
  for (out in 1:k) { # iterating over folds
    
    # split the data into CV training and test sets
    train <- data[A!=folds[out],]
    test <- data[A==folds[out],]
    
    # build the models
    model <- anchor.regression(train$X, train$Y, train$A, gamma.CV, nrow(train))
    
    # make predictions on CV test set and compute test MSE
    predictions <- test$X * model$coefficients
    quantile.CV.matrix[a,out] <- quantile((test$Y - predictions) ^ 2, alpha)
  }
}
quantile.CV <- apply(quantile.CV.matrix, 1, mean)
quantile.CV

# step 5: choose optimal gamma
alpha.vec[which(quantile.CV==min(quantile.CV))]
gamma.optimal <- qchisq(alpha.vec[which(quantile.CV==min(quantile.CV))], df=1)
gamma.optimal

plot(alpha.vec,quantile.CV)
plot(qchisq(alpha.vec, df=1),quantile.CV)

##########################################################################
# Fit AR with optimal gamma
fit<- anchor.regression(X.train, Y.train, A, gamma.optimal, n)
b <- coef(fit)

MSE.train <- mean((Y.train - t(X.train)*b) ^ 2)

##########################################################################
# Test CV output on data with different pertubations
v.vec <- seq(-5,5, by=0.1)

nit <- 100
MSE.test <- matrix(nrow=length(v.vec), ncol=nit)
MSE.test.OLS <- matrix(nrow=length(v.vec), ncol=nit)
MSE.test.PA <- matrix(nrow=length(v.vec), ncol=nit)
MSE.test.IV <- matrix(nrow=length(v.vec), ncol=nit)
for (vi in 1:length(v.vec)) {
  v <- v.vec[vi]
  
  for (i in 1:nit) {
    # initialize test data
    epsH.test <- rnorm(n)
    epsX.test <- rnorm(n)
    epsY.test <- rnorm(n)
    
    H.test <- epsH.test
    X.test <- v+H.test+epsX.test
    Y.test <- X.test+2*H.test+epsY.test
    
    MSE.test[vi,i] <- mean((Y.test - t(X.test)*b) ^ 2)
    MSE.test.OLS[vi,i] <- mean((Y.test - t(X.test)*b.OLS) ^ 2)
    MSE.test.PA[vi,i] <- mean((Y.test - t(X.test)*b.PA) ^ 2)
    MSE.test.IV[vi,i] <- mean((Y.test - t(X.test)*b.IV) ^ 2)
  }
}

MSE.test <- apply(MSE.test, 1, mean)
MSE.test.OLS <- apply(MSE.test.OLS, 1, mean)
MSE.test.PA <- apply(MSE.test.PA, 1, mean)
MSE.test.IV <- apply(MSE.test.IV, 1, mean)

plot(v.vec,MSE.test, type = "l", col=2, ylim = c(3,8), ylab = "MSE", xlab = "v", main = "MSE for varying shifts v")
abline(h=MSE.train)
lines(v.vec, MSE.test.OLS, col = 3)
lines(v.vec, MSE.test.PA, col = 4)
lines(v.vec, MSE.test.IV, col = 5)
legend(2.8, 8, legend=c("CV train MSE", "CV", "OLS", "PA", "IV"),
       col=c(1, 2, 3, 4, 5), cex=0.8,pch=16)

plot(v.vec,MSE.test, type = "l", col=2, ylim = c(3,8), xlim= c(-3,3), ylab = "MSE", xlab = "v", main = "MSE for varying shifts v")
abline(h=MSE.train)
lines(v.vec, MSE.test.OLS, col = 3)
lines(v.vec, MSE.test.PA, col = 4)
lines(v.vec, MSE.test.IV, col = 5)
legend(1.6, 8, legend=c("CV train MSE", "CV", "OLS", "PA", "IV"),
       col=c(1, 2, 3, 4, 5), cex=0.8,pch=16)

##########################################################################
# Comparison to specific shifts
v.spec <- c(1.8,-1.8,3.5)

nit <- 100
b.test <- matrix(nrow=length(v.spec), ncol=nit)
b.OLS <- matrix(nrow=length(v.spec), ncol=nit)
b.PA <- matrix(nrow=length(v.spec), ncol=nit)
b.IV <- matrix(nrow=length(v.spec), ncol=nit)
for (vi in 1:length(v.spec)) {
  v <- v.vec[vi]
  
  for (i in 1:nit) {
    # initialize test data
    epsH.test <- rnorm(n)
    epsX.test <- rnorm(n)
    epsY.test <- rnorm(n)
    
    H.test <- epsH.test
    X.test <- v+H.test+epsX.test
    Y.test <- X.test+2*H.test+epsY.test
    
    fit <- anchor.regression(X.test, Y.test, rep(v,length(Y.test)), gamma.optimal, length(Y.test))
    b.test[vi,i] <- coef(fit)
    
    fit.OLS <- anchor.regression(X.test, Y.test, rep(v,length(Y.test)), 1, length(Y.test))
    b.OLS[vi,i] <- coef(fit.OLS)
    
    fit.PA <- anchor.regression(X.test, Y.test, rep(v,length(Y.test)), 0, length(Y.test))
    b.PA[vi,i] <- coef(fit.PA)
    
    fit.2SLS.step <- lm(X.test~rep(v,length(Y.test))-1)
    X.test.hat <- fitted.values(fit.2SLS.step)
    fit.IV <- lm(Y.test~X.test.hat-1)
    b.IV[vi,i] <- coef(fit.IV)
  }
}

apply(b.test, 1, mean)
apply(b.OLS, 1, mean)
apply(b.PA, 1, mean)
apply(b.IV, 1, mean)