# Investigating anchor approach using likelihood
#
# Author: Maic Rakitta
# Date: 28.09.20
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

# Objective data
X <- X.train
Y <- Y.train
H <- H.train
A <- A

# Assuming linear model for Anchor dependency
P.A <- A%*%solve(t(A)%*%A)%*%t(A)

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
# Fitting AR with objective fct
objec.fct <- function(b){
  
  return(sum((Y-b*X)^2)+(gamma-1)*sum(((P.A%*%Y-b*P.A%*%X))^2))
}

gamma <- 0.5
optim <- optimize(objec.fct, interval = c(-20,20))
optim$minimum

fit <- anchor.regression(X, Y, A, gamma, n)
b <- coef(fit)
b

# matrix notation
objec.fct <- function(b){
  
  return(t(Y-b*X)%*%(Y-b*X)+(gamma-1)*t(Y-b*X)%*%P.A%*%(Y-b*X))
}



##########################################################################
# Fitting AGLM for binomial data with likelihood objective

# Binomial Data
psi <- 0.5
##---Example 2: logistic model and no covariates---
#m <- plogis(1+0.8*X-0.39*H)
#Y <- rbinom(n, 1, plogis(psi*X+log(m/(1-m))))
Y <- rbinom(n, 1, plogis(psi*X))

data <- data.frame(A, X, Y)

objec.binom <- function(b){
  
  return(-sum(Y*b*X-log(1+exp(b*X)))-(gamma-1)*sum(P.A%*%Y*b*P.A%*%X-log(1+exp(b*P.A%*%X))))
}

# classic MLE
gamma <- 1
optim <- optimize(objec.binom, interval = c(-20,20))
optim$minimum

fit <- glm(Y~X-1, family="binomial")
coef(fit)

# Instrumental variables with GLM
gamma <- 1000000
optim <- optimize(objec.binom, interval = c(-20,20))
optim$minimum

library(ivtools)
fitX.HA <- glm(formula=X~A-1, family="gaussian", data=data) #two-stage estimation
fitY.HX <- glm(formula=Y~X-1, family="binomial", data=data)
fitIV <- ivglm(estmethod="ts", fitX.LZ=fitX.HA, fitY.LX=fitY.HX, data=data,
               ctrl=F)
summary(fitIV)












### IN CONSTRUCTION ###

##########################################################################
# Fitting AGLM for poisson data with likelihood objective

# Poisson Data

### TO INSERT HERE ####

objec.poisson <- function(b){

  return(-sum(Y*b*X-exp(b*X))-(gamma-1)*sum(P.A%*%Y*b*P.A%*%X-exp(b*P.A%*%X)))
}

# classic MLE
gamma <- 1
optim <- optimize(objec.poisson, interval = c(-20,20))
optim$minimum

fit <- glm(Y~X-1, family="poisson")
coef(fit)
