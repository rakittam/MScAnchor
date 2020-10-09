# Anchor regression using GLM approach and CCXR, optim, alabama
# Binary logistic regression
#
# Author: Maic Rakitta
# Date: 07.10.20
##########################################################################
set.seed(1)

n <- 1000 # number of samples from unpertubed and pertubed distribution

# initialize data
library(extraDistr) # for rademacher distribution

m <- 5

# sample anchor coefficients
g1 <- rnorm(n=1)
g2 <- rnorm(n=1)

# initialize training data
A.train <- matrix(nrow = n, ncol = 2)
H.train <- matrix(nrow = n, ncol = 1)
X.train <- matrix(nrow = n, ncol = 10)
Y.train <- matrix(nrow = n, ncol = 1)

for (i in 1:n) {
  
  A.train[i,] <- rsign(n=2)
  
  epsH.train <- rnorm(n=1, mean=0, sd=1)
  H.train[i] <- epsH.train
  
  epsX.train <- rnorm(n=10, mean=0, sd=1)
  X.train[i,] <- g1*A.train[i,1]+g2*A.train[i,2]+H.train[i]+epsX.train
  
  Y.train[i] <- rbinom(n=1, size=m, plogis(3*X.train[i,2]+3*X.train[i,3]+H.train[i]-2*A.train[i,1]))
}

# Objective data
X <- X.train[,2:3]
Y <- Y.train
H <- H.train
A <- A.train

# Orthogonal projection on column space of anchor A
P.A <- A%*%solve(t(A)%*%A)%*%t(A)

##########################################################################
# Anchor GLM for binary logistic regression
library(CVXR)

AGLM_normal <- function(xi){
  
  # Step 1. Define the variable to be estimated
  b.hat <- Variable(1) 
  # Step 2. Define the objective to be optimized
  
  #loss <- -1/n*sum(Y*X%*%b.hat-m*log(1+exp(X%*%b.hat)))
  loss <- -sum(X[Y==1]%*%b.hat)+sum(m*logistic(X%*%b.hat))
  
  anchor_penalty <- function(b.hat){
    p.hat <- exp(X%*%b.hat)/(1+exp(X%*%b.hat)) # inverse of logit link
    r.D <- (Y/m-p.hat)/abs(Y/m-p.hat)*sqrt(2*(Y*log(Y/(m*p.hat))+(m-Y)*log((m-Y)/(m-m*p.hat)))) # deviance residuals
    return(quad_form(r.D, P.A))
  }
  
  #objective <- 1/n*loss
  objective <- 1/n*(loss + xi * anchor_penalty(b.hat))
  
  # Step 3. Create a problem to solve
  problem <- Problem(Minimize(objective))
  # Step 4. Solve it!
  result <- solve(problem)
  # Step 5. Extract solution and objective value
  b.AGLM <- result$getValue(b.hat)
  b.AGLM
  
  return(b.AGLM)
}

##########################################################################
# Playground

gamma <- 2 # Set gamma and xi
xi <- gamma-1

AGLM_normal(xi) # Run AGLM for normal relation

##########################################################################
# Using alabama, optim, optimize
xi=1
# Step 2. Define the objective to be optimized
loss <- function(b.hat){
  return(-sum(Y*(X%*%b.hat)-m*log(1+exp(X%*%b.hat))))
}

anchor_penalty <- function(b.hat){
  p.hat <- exp(X%*%b.hat)/(1+exp(X%*%b.hat)) # inverse of logit link
  
  special.case1 <- function(Y){
    ifelse(Y==0, 0, Y*log(Y/(m*p.hat)))
  }
  special.case2 <- function(Y){
    ifelse(Y==m, 0, (m-Y)*log((m-Y)/(m-m*p.hat)))
  }
  
  #r.D <- (Y/m-p.hat)/abs(Y/m-p.hat)*sqrt(2*(Y*log(Y/(m*p.hat))+(m-Y)*log((m-Y)/(m-m*p.hat)))) # deviance residuals
  r.D <- (Y/m-p.hat)/abs(Y/m-p.hat)*sqrt(2*(special.case1(Y)+special.case2(Y))) # deviance residuals
  return(t(r.D)%*%P.A%*%r.D)
}

objective <- function(b.hat){
  return(1/n*(loss(b.hat) + xi * anchor_penalty(b.hat)))
}

# For one dimensional unconstrained optimization
#ans1 <- optimize(f = objective, interval = c(-20,20))
#ans1

# For multidimensional unconstrained optimization
start.val <- c(1,1)
ans2 <- optim(par=start.val, fn=objective, hessian = TRUE)
ans2$par

# For constrained optimization
#ans3 <- auglag(par=c(1,1), fn=objective)
#ans3