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
# Fitting AGLM with CVXR
library(CVXR)

AGLM_normal <- function(xi){
  
  # Step 1. Define the variable to be estimated
  b.hat <- Variable(ncol(X)) 
  
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


##########################################################################
# Using alabama, optim, optimize
xi=1
# Step 2. Define the objective to be optimized
loss <- function(b.hat){
  return(sum((Y-X%*%b.hat)^2))
}

anchor_penalty <- function(b.hat){
  p.hat <- X%*%b.hat # inverse of identity link
  r.D <- sqrt(2)*(Y-p.hat)  # deviance residuals
  return(t(r.D)%*%P.A%*%r.D)
}
objective <- function(b.hat){
  return(1/n*(loss(b.hat) + xi * anchor_penalty(b.hat)))
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
