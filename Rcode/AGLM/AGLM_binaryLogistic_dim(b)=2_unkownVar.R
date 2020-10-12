# Anchor regression using GLM approach and CCXR, optim, alabama
# Normal linear regression
#
# Author: Maic Rakitta
# Date: 07.10.20
##########################################################################
set.seed(1)

n <- 300 # number of samples from unpertubed and pertubed distribution

# sample anchor coefficients
g1 <- rnorm(n=1)
g2 <- rnorm(n=1)

# initialize training data
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
  Y.train[i] <- 3*X.train[i,2]+3*X.train[i,3]+H.train[i]-2*A.train[i,1]+epsY.train
}

# Objective data
X <- X.train[,2:3]
Y <- Y.train
H <- H.train
A <- A.train

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
# Using alabama, optim, optimize

AGLM <- function(xi){
  # hat.vec is vector that combines parameter to be estimated b.hat and sigma^2
  loss <- function(hat.vec){
    b.hat <- hat.vec[1:2]
    s2 <- hat.vec[3]
    return(-n/2*log(2*pi*s2)*-1/(2*s2)*sum((Y-X%*%b.hat)^2))
  }
  
  anchor_penalty <- function(hat.vec){
    b.hat <- hat.vec[1:2]
    s2 <- hat.vec[3]
    p.hat <- X%*%b.hat # inverse of identity link
    r.D <- 1/sqrt(s2)*(Y-p.hat)  # deviance residuals
    return(t(r.D)%*%P.A%*%r.D)
  }
  objective <- function(hat.vec){
    return(1/n*(-loss(hat.vec) + xi * anchor_penalty(hat.vec)))
  }
  
  # For one dimensional unconstrained optimization
  #ans1 <- optimize(f = objective, interval = c(-20,20))
  #ans1
  
  # For multidimensional unconstrained optimization
  ans2 <- optim(par=c(1,1,1), fn=objective, hessian = TRUE)
  b.AGLM <- ans2$par
  return(b.AGLM)
  hess.mat <- ans2$hessian
  
  # Is b.AGLM local minimum?
  det(hess.mat)>0 & hess.mat[1,1]>0
  
  # For constrained optimization
  #ans3 <- auglag(par=c(1,1), fn=objective)
  #ans3
}

xi=1
AGLM(xi)
