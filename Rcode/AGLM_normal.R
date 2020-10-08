# Anchor regression using GLM approach and CCXR package
#
# Author: Maic Rakitta
# Date: 07.10.20
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

# Orthogonal projection on column space of anchor A
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
# Fitting AGLM with CVXR
library(CVXR)

AGLM_normal <- function(xi){
  
  # Step 1. Define the variable to be estimated
  b.hat <- Variable(1) 
  
  # Step 2. Define the objective to be optimized
    loss <- sum((Y-X%*%b.hat)^2)
  
  anchor_penalty <- function(b.hat){
    p.hat <- X%*%b.hat # inverse of identity link
    r.D <- sqrt(2)*(Y-p.hat)  # deviance residuals
    return(quad_form(r.D, P.A))
  }
  
  objective <- 1/n*(loss + xi * anchor_penalty(b.hat))
  
  # Step 3. Create a problem to solve
  problem <- Problem(Minimize(objective))
  
  # Step 4. Solve it!
  result <- solve(problem)
  
  # Step 5. Extract solution and objective value
  b.AGLM <- result$getValue(b.hat)
  
  return(b.AGLM)
}

# Set gamma / xi
gamma <- 2
xi <- gamma-1

# Run AGLM for normal relation
AGLM_normal(xi)

# AR
fit <- anchor.regression(X, Y, A, gamma, n)
b.AR <- coef(fit)
b.AR
