# Optimal gamma calculations using k-fold cross validation
# in example 2 in Anchor Regression - Rothenhaeuser 2020
#
# Author: Maic Rakitta
# Date: 14.10.20
##########################################################################
set.seed(1)

n <- 1000 # Number of samples from unpertubed and pertubed distribution

# Initialize training data
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

##########################################################################
# Calculating parameter beta for specific cases

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
# k-fold Cross Validation for anchor regression
data <- data.frame(Y=Y.train,X=X.train, A=A) # create data frame for CV

# Initialize
k <- 2
alpha.vec <- seq(0,1,by = 0.01)
gamma.vec <- seq(0,3,by = 0.1)

result.array <- array(dim=c(length(gamma.vec),length(alpha.vec),k))

for (g in 1:length(gamma.vec)) {
  
  gamma.temp <- gamma.vec[g]
  for (a in 1:length(alpha.vec)) {
    
    alpha <- alpha.vec[a]
    folds <- c(-1,1)
    for (out in 1:k) {
      
      # Split the data into CV training and test sets
      train <- data[A!=folds[out],]
      test <- data[A==folds[out],]
      
      # Build the models
      model <- anchor.regression(train$X, train$Y, train$A, gamma.temp, nrow(train))
      
      # Make predictions on CV test set and compute test quantile
      predictions <- test$X * model$coefficients
      result.array[g,a,out] <- quantile((test$Y - predictions) ^ 2, alpha)
      
    }
  }
}

averaged.folds <- apply(result.array,1:2,mean) # average over folds for each gamma alpha combo
min.index <- apply(averaged.folds,2,which.min) # find minimum for each alpha
gamma.optimal.vec <- gamma.vec[min.index]
gamma.optimal.vec

plot(alpha.vec, gamma.optimal.vec) # plot optimal gamma for each alpha
plot(apply(averaged.folds,2,min)) # plot min for each alpha

# Choose specific alpha and give corresponding optimal gamma
alpha <- 0.9
gamma.optimal <- gamma.optimal.vec[which(alpha.vec==alpha)]
gamma.optimal

##########################################################################
# Fit AR with optimal gamma on hole training data
fit <- anchor.regression(X.train, Y.train, A, gamma.optimal, n)
b <- coef(fit)
b

# Train MSE with optimal gamma chosen by CV
quantile.train <- quantile((Y.train - t(X.train)*b) ^ 2, alpha)
quantile.train

# Compute corresponding optimal beta for CV chosen optimal gamma
beta.optimal.vec <- numeric(length(gamma.optimal.vec))
for (g in 1:length(gamma.optimal.vec)) {
  gamma <- gamma.optimal.vec[g]
  beta.optimal.vec[g] <- coef(anchor.regression(X.train, Y.train, A, gamma, n))
}

##########################################################################
# Test CV output on data with different pertubations
v.vec <- seq(-10,10, length.out = 100)
nit <- 100

n <- 10000

# Initialize result arrays
quantile.test.CV.array <- array(dim=c(length(alpha.vec), length(v.vec), nit))
quantile.test.noCV.array <- array(dim=c(3, length(v.vec), nit))
dimnames(quantile.test.noCV.array)[[1]] <- c("OLS", "PA", "IV")

for (vi in 1:length(v.vec)) {
  v <- v.vec[vi]
  
  for (i in 1:nit) {
    # Initialize by CV unseen test data
    epsH.test <- rnorm(n)
    epsX.test <- rnorm(n)
    epsY.test <- rnorm(n)
    
    H.test <- epsH.test
    X.test <- v+H.test+epsX.test
    Y.test <- X.test+2*H.test+epsY.test
    
    for (a in 1:length(alpha.vec)) {
      alpha <- alpha.vec[a]
      beta <- beta.optimal.vec[a]
      quantile.test.CV.array[a, vi, i] <- quantile((Y.test - X.test*beta) ^ 2, alpha)
    }
    
    quantile.test.noCV.array["OLS", vi, i] <- quantile((Y.test - X.test*b.OLS) ^ 2, alpha)
    quantile.test.noCV.array["PA", vi, i] <- quantile((Y.test - X.test*b.PA) ^ 2, alpha)
    quantile.test.noCV.array["IV", vi, i] <- quantile((Y.test - X.test*b.IV) ^ 2, alpha)
  }
}

quantile.test.CV.averaged.iterations <- apply(quantile.test.CV.array,1:2,mean)
quantile.test.noCV.averaged.iterations <- apply(quantile.test.noCV.array,1:2,mean)

index0.9 <- which(alpha.vec==0.9)
plot(v.vec, quantile.test.CV.averaged.iterations[index0.9, ],
     type="l", ylim = c(0, 120),
     xlab = "v",
     ylab= "Quantile of test squared residual",
     main = "Benchmarking of different interventions v for alpha = 0.9")
lines(v.vec, quantile.test.noCV.averaged.iterations["OLS",], col=2)
lines(v.vec, quantile.test.noCV.averaged.iterations["PA",], col=3)
lines(v.vec, quantile.test.noCV.averaged.iterations["IV",], col=4)
legend(-1, 120, legend=c("optimal gamma", "OLS", "PA", "IV"),
       col=c(1, 2, 3, 4), cex=0.8, pch=16)

# Playground for different alphas
index.temp <- which(alpha.vec==0.8)
lines(v.vec, quantile.test.CV.averaged.iterations[index.temp, ], col = 5)
