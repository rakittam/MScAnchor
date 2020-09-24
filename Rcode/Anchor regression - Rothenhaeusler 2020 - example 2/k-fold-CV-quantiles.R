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
alpha.vec <- (1:100)/101

MSE.CV <- numeric(length(alpha.vec))
gamma.optimal <- numeric(length(alpha.vec))
for (a in 1:length(alpha.vec)) {
  
  alpha <- alpha.vec[a] # step 1: choose alpha
  gamma.CV <- qchisq(alpha, df=1)
  
  folds <- sample(1:k, nrow(data), replace=T) # step 2: create folds
  
  # step 3: for varying gamma train and test
  MSE.CV.matrix <- numeric(k)
  for (out in 1:k) { # iterating over folds
    
    # split the data into CV training and test sets
    train <- data[folds!=out,]
    test <- data[folds==out,]
    
    # build the models
    model <- anchor.regression(train$X, train$Y, train$A, gamma.CV, nrow(train))
    
    # make predictions on CV test set and compute test MSE
    predictions <- test$X * model$coefficients
    MSE.CV.matrix[out] <- mean((test$Y - predictions) ^ 2)
  }
  # step 4: average over the folds
  MSE.CV[a] <- mean(MSE.CV.matrix)
}
MSE.CV

# step 5: choose optimal gamma
alpha.vec[which(MSE.CV==min(MSE.CV))]
gamma.optimal <- qchisq(alpha.vec[which(MSE.CV==min(MSE.CV))], df=1)
gamma.optimal

plot(alpha.vec,MSE.CV)
plot(qchisq(alpha.vec, df=1),MSE.CV)
