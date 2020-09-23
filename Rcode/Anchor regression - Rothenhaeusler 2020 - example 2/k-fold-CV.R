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

# test MSE for specific cases
MSE.PA <- mean((Y.test - t(X.test)*b.PA) ^ 2)
MSE.OLS <- mean((Y.test - t(X.test)*b.OLS) ^ 2)
MSE.IV <- mean((Y.test - t(X.test)*b.IV) ^ 2)

##########################################################################
# k-fold Cross Validation 

data <- data.frame(Y=Y.train,X=X.train, A=A) # create data frame for CV

# initialize
k <- 10
gamma.CV <- seq(0,5,by=0.01)
MSE.CV.matrix <- matrix(nrow = length(gamma.CV), ncol = k)

for (g in 1:length(gamma.CV)) { # iterating over different gammas

  folds <- sample(1:k, nrow(data), replace=T)
  
  for (out in 1:k) { # iterating over folds
    
    # split the data into CV training and test sets
    train <- data[folds!=out,]
    test <- data[folds==out,]
    
    # build the models
    model <- anchor.regression(train$X, train$Y, train$A, gamma.CV[g], nrow(train))
      
    # make predictions on CV test set and compute test MSE
    predictions <- test$X * model$coefficients
    MSE.CV.matrix[g,out] <- mean((test$Y - predictions) ^ 2)
  }
}

MSE.CV <- apply(MSE.CV.matrix,1,mean) # average over the k folds

# chosen optimal gamma with CV
gamma.opt.CV <- gamma.CV[which(MSE.CV==min(MSE.CV))]
gamma.opt.CV

# plot CV test MSE
plot(gamma.CV,MSE.CV)
abline(v=gamma.opt.CV, col="forestgreen")

##########################################################################
# Compare CV MSE to test MSE of pertubed data

MSE.pert <- numeric(length(gamma.CV))

for (i in 1:length(gamma.CV)){
  
  gamma <- gamma.CV[i]
  
  fit <- anchor.regression(X.train, Y.train, A, gamma, n)
  
  MSE.pert[i] <- mean((Y.test - t(X.test)*coef(fit)) ^ 2)
}

# true optimal gamma for given shift
gamma.opt.pert <- gamma.CV[which(MSE.pert==min(MSE.pert))]
gamma.opt.pert

# performance plots
plot(gamma.CV,MSE.CV)
plot(gamma.CV,MSE.pert)

plot(gamma.CV,MSE.CV, type = "l", ylim = c(min(min(MSE.CV),min(MSE.pert)),max(max(MSE.CV),max(MSE.pert))))
lines(gamma.CV,MSE.pert, col = 2, lty = "dashed")
abline(v=gamma.opt.CV, col="forestgreen")
abline(v=gamma.opt.pert, col="blue")
legend(3.5, 6.5, legend=c("CV MSE", "true MSE", "CV optimal g", "true optimal g"),
       col=c(1, 2, "forestgreen", "blue"), cex=0.8,pch=16)
