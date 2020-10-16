# Anchor regression using GLM approach and CCXR, optim, alabama
# Normal linear regression
#
# Author: Maic Rakitta
# Date: 07.10.20
##########################################################################
set.seed(1)

n <- 300 # number of samples from unpertubed and pertubed distribution

# Anchor coefficients
g1 <- rnorm(n=1)
g2 <- rnorm(n=1)
g3 <- -2

# Initialize training data
A.train <- matrix(nrow = n, ncol = 2)
H.train <- matrix(nrow = n, ncol = 1)
X.train <- matrix(nrow = n, ncol = 10)
Y.train <- matrix(nrow = n, ncol = 1)

for (i in 1:n) {
  
  A.train[i,] <- rnorm(n=2, mean=0, sd=1)
  
  epsH.train <- rnorm(n=1, mean=0, sd=1)
  H.train[i] <- epsH.train
  
  epsX.train <- rnorm(n=10, mean=0, sd=1)
  X.train[i,] <- g1*A.train[i,1]+g2*A.train[i,2]+H.train[i]+epsX.train
  
  epsY.train <- rnorm(n=1, mean=0, sd=0.25^2)
  Y.train[i] <- 3*X.train[i,2]+3*X.train[i,3]+H.train[i]+g3*A.train[i,1]+epsY.train
}

# Objective data
X <- X.train[,2:3]
Y <- Y.train
H <- H.train
A <- A.train

# Orthogonal projection on column space of anchor A
P.A <- A%*%solve(t(A)%*%A)%*%t(A)

##########################################################################
# Linear Anchor Regression
anchor.regression <- function(X, Y, A, gamma, n){
  
  P.A <- A%*%solve(t(A)%*%A)%*%t(A)
  W <- diag(n)-(1-sqrt(gamma))*P.A
  
  Y.tilde <- W%*%Y
  X.tilde <- W%*%X
  
  fit <- lm(Y.tilde~X.tilde-1)
}

##########################################################################
# AGLM with CVXR
library(CVXR)

AGLM_CVXR <- function(xi){
  
  # Step 1. Define the variable to be estimated
  b.hat <- Variable(ncol(X)) 
  
  # Step 2. Define the objective to be optimized
  loss <- -sum((Y-X%*%b.hat)^2)
  
  anchor_penalty <- function(b.hat){
    p.hat <- X%*%b.hat # inverse of identity link
    r.D <- Y-p.hat  # deviance residuals
    return(quad_form(r.D, P.A))
  }
  
  objective <- 1/n*(-loss + xi * anchor_penalty(b.hat))
  
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

# Classic linear AR Method
fit <- anchor.regression(X, Y, A, gamma, n)
b.AR <- coef(fit)
b.AR

# CVXR AR
b.AGLM_CVXR <- AGLM_CVXR(xi)
b.AGLM_CVXR

##########################################################################
# AGLM using alabama, optim, optimize

AGLM <- function(xi){
  # Step 1. Define the objective loss
  loss <- function(b.hat){
    return(-sum((Y-X%*%b.hat)^2))
  }
  
  # Step 2. Define anchor penalty
  anchor_penalty <- function(b.hat){
    p.hat <- X%*%b.hat # inverse of identity link
    r.D <- Y-p.hat  # deviance residuals
    return(t(r.D)%*%P.A%*%r.D)
  }
  
  # Step 3. Contruct objective by 1. and 2.
  objective <- function(b.hat){
    return(1/n*(-loss(b.hat) + xi * anchor_penalty(b.hat)))
  }
  
  # For one dimensional unconstrained optimization
  #ans1 <- optimize(f = objective, interval = c(-20,20))
  #ans1
  
  # For multidimensional unconstrained optimization
  ans2 <- optim(par=c(1,1), fn=objective, hessian = TRUE)
  b.AGLM <- ans2$par
  b.AGLM
  hess.mat <- ans2$hessian
  
  # Is b.AGLM local minimum?
  det(hess.mat)>0 & hess.mat[1,1]>0
  
  # For constrained optimization
  #ans3 <- auglag(par=c(1,1), fn=objective)
  #ans3
  return(b.AGLM)
}

b.AGLM <- AGLM(1)
b.AGLM

##########################################################################
# Iterating over different hyper parameters for Rothenhaeusler e2 plot

g1.test <- 1
g2.test <- 1
g3.test <- -2
  
# Initialize test data
A.test <- matrix(nrow = n, ncol = 2)
H.test <- matrix(nrow = n, ncol = 1)
X.test <- matrix(nrow = n, ncol = 10)
Y.test <- matrix(nrow = n, ncol = 1)

for (i in 1:n) {
  A.test[i,] <- rnorm(n=2, mean=0, sd=1)
  epsH.test <- rnorm(n=1, mean=0, sd=1)
  H.test[i] <- epsH.test
  epsX.test <- rnorm(n=10, mean=0, sd=1)
  X.test[i,] <- g1.test*A.test[i,1]+g2.test*A.test[i,2]+H.test[i]+epsX.test
  epsY.test <- rnorm(n=1, mean=0, sd=0.25^2)
  Y.test[i] <- 3*X.test[i,2]+3*X.test[i,3]+H.test[i]+g3.test*A.test[i,1]+epsY.test
}

# Objective data
X.test <- X.test[,2:3]
Y.test <- Y.test
H.test <- H.test
A.test <- A.test

# Iterating over xi
xi.vec <- seq(-1,10,by=0.1)
b.AGLM.matrix <- matrix(nrow=length(xi.vec), ncol = 2)
deviance.vec <- numeric(length(xi.vec))

for (i in 1:length(xi.vec)) {
  xi <- xi.vec[i]
  b.AGLM.matrix[i,] <- AGLM(xi)
  
  p.hat.test <- X.test%*%b.AGLM.matrix[i,]# inverse of logit link
  
  r.D.test <- Y-p.hat.test # deviance residuals
  deviance.vec[i] <- 1/n*t(r.D.test)%*%r.D.test
}

# Plot like in ex2 of Rothenhaeusler
plot(xi.vec, deviance.vec, type = "l")
xi.vec[which.min(deviance.vec)]

plot(b.AGLM.matrix[,1], deviance.vec)
plot(b.AGLM.matrix[,2], deviance.vec)

plot(b.AGLM.matrix[,1], deviance.vec, type="l")
lines(b.AGLM.matrix[,2], deviance.vec)
