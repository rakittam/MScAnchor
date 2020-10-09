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
b.CVXR <- AGLM_normal(xi)
b.CVXR

# AR
fit <- anchor.regression(X, Y, A, gamma, n)
b.AR <- coef(fit)
b.AR


##########################################################################
# Using alabama, optim, optimize
xi=1
# Step 2. Define the objective to be optimized
loss <- function(b.hat){
  return(sum((Y-X*b.hat)^2))
}

anchor_penalty <- function(b.hat){
  p.hat <- X*b.hat # inverse of identity link
  r.D <- sqrt(2)*(Y-p.hat)  # deviance residuals
  return(t(r.D)%*%P.A%*%r.D)
}
objective <- function(b.hat){
  return(1/n*(loss(b.hat) + xi * anchor_penalty(b.hat)))
}

# For one dimensional unconstrained optimization
ans1 <- optimize(f = objective, interval = c(-20,20))
ans1

# For multidimensional unconstrained optimization
ans2 <- optim(par=1, fn=objective, method = "Brent", lower = -20, upper = 20, hessian = TRUE)
b.AGLM <- ans2$par
b.AGLM
hess.mat <- ans2$hessian

# Is b.AGLM local minimum?
det(hess.mat)>0 & hess.mat[1,1]>0

# For constrained optimization
#ans3 <- auglag(par=p0, fn=objective, heq=heq, hin=hin, gr=gr)
#ans3






# 1 dimensional b - gradient and hessian
gr <- function(b.hat) {
  g <- 1/n * (2*b.hat*t(X)%*%X - t(Y)%*%X + 2*xi * (-t(X)%*%P.A%*%Y + b.hat* t(X)%*%P.A%*%X))
  g
}

hessian <- function(b.hat) {
  hess <- 2/n * (t(X)%*%X + t(X)%*%P.A%*%X)
  hess
}