# Anchor regression using GLM approach and CCXR, optim, alabama
# Binary logistic regression
#
# Author: Maic Rakitta
# Date: 07.10.20
##########################################################################
set.seed(1)

n <- 1000 # number of samples from unpertubed and pertubed distribution

# Initialize data
library(extraDistr) # for rademacher distribution

# Anchor coefficients
g1 <- rnorm(n=1)
g2 <- rnorm(n=1)
g3 <- -2

# initialize training data
A.train <- matrix(nrow = n, ncol = 1)
H.train <- matrix(nrow = n, ncol = 1)
X.train <- matrix(nrow = n, ncol = 10)
Y.train <- matrix(nrow = n, ncol = 1)

m <- 5 # number of trials for binary distribution

for (i in 1:n) {
  
  A.train[i] <- rsign(n=1)
  
  epsH.train <- rnorm(n=1, mean=0, sd=1)
  H.train[i] <- epsH.train
  
  epsX.train <- rnorm(n=10, mean=0, sd=1)
  X.train[i,] <- g1*A.train[i]+H.train[i]+epsX.train
  
  Y.train[i] <- rbinom(n=1, size=m, plogis(3*X.train[i,2]+H.train[i]))
}

# Objective data
X <- X.train[,2]
Y <- Y.train
H <- H.train
A <- A.train

# Orthogonal projection on column space of anchor A
P.A <- A%*%solve(t(A)%*%A)%*%t(A)

##########################################################################
# Anchor GLM for binary logistic regression

# library(CVXR)
# 
# AGLM_CVXR <- function(xi){
# 
#   # Step 1. Define the variable to be estimated
#   b.hat <- Variable(1)
#   # Step 2. Define the objective to be optimized
# 
#   #loss <- -1/n*sum(Y*X%*%b.hat-m*log(1+exp(X%*%b.hat)))
#   loss <- -sum(X[Y==1]%*%b.hat)+sum(m*logistic(X%*%b.hat))
# 
#   anchor_penalty <- function(b.hat){
#     p.hat <- exp(X%*%b.hat)/(1+exp(X%*%b.hat)) # inverse of logit link
#     r.D <- (Y/m-p.hat)/abs(Y/m-p.hat)*sqrt(2*(Y*log(Y/(m*p.hat))+(m-Y)*log((m-Y)/(m-m*p.hat)))) # deviance residuals
#     return(quad_form(r.D, P.A))
#   }
# 
#   #objective <- 1/n*loss
#   objective <- 1/n*(loss + xi * anchor_penalty(b.hat))
# 
#   # Step 3. Create a problem to solve
#   problem <- Problem(Minimize(objective))
#   # Step 4. Solve it!
#   result <- solve(problem)
#   # Step 5. Extract solution and objective value
#   b.AGLM <- result$getValue(b.hat)
#   b.AGLM
# 
#   return(b.AGLM)
# }

##########################################################################
# AGLM for Logistic regression data using optim as optimizer
AGLM <- function(xi){
  
  # Step 1. Define the objective loss
  loss <- function(b.hat){
    return(-sum(Y*(X*b.hat)-m*log(1+exp(X*b.hat))))
  }
  
  # Step 2. Define anchor penalty
  anchor_penalty <- function(b.hat){
    p.hat <- exp(X*b.hat)/(1+exp(X*b.hat)) # inverse of logit link
    
    special.case1 <- function(Y){
      ifelse(Y==0, 0, Y*log(Y/(m*p.hat)))
    }
    special.case2 <- function(Y){
      ifelse(Y==m, 0, (m-Y)*log((m-Y)/(m-m*p.hat)))
    }
    
    #r.D <- (Y/m-p.hat)/abs(Y/m-p.hat)*sqrt(2*(Y*log(Y/(m*p.hat))+(m-Y)*log((m-Y)/(m-m*p.hat)))) # deviance residuals
    r.D <- sign(Y/m-p.hat)*sqrt(2*(special.case1(Y)+special.case2(Y))) # deviance residuals
    return(t(r.D)%*%P.A%*%r.D)
  }
  
  # Step 3. Contruct objective by 1. and 2.
  objective <- function(b.hat){
    return(1/n*(loss(b.hat) + xi * anchor_penalty(b.hat)))
  }
  
  # Set start value for optimization
  start.val <- c(1)
  ans2 <- optimize(interval = c(-10,10), f= objective)
  return(ans2$minimum)
}

##########################################################################
# Playground
gamma <- 2 # Set gamma and xi
xi <- gamma-1
AGLM(xi)

##########################################################################
# Iterating over different hyper parameters for Rothenhaeusler e2 plot

g1.test <- 0 # 0, i.e. chooses MLE, 
g1.test <- 0.2 # close to xi=1 => gamma = 2
g1.test <- -0.2 # 0.1
g1.test <- -0.9 # 4.9

# Initialize test data
A.test <- matrix(nrow = n, ncol = 1)
H.test <- matrix(nrow = n, ncol = 1)
X.test <- matrix(nrow = n, ncol = 10)
Y.test <- matrix(nrow = n, ncol = 1)

for (i in 1:n) {
  A.test[i] <- rsign(n=1)
  epsH.test <- rnorm(n=1, mean=0, sd=1)
  H.test[i] <- epsH.test
  epsX.test <- rnorm(n=10, mean=0, sd=1)
  X.test[i,] <- g1.test*A.test[i]+H.test[i]+epsX.test
  Y.test[i] <- rbinom(n=1, size=m, plogis(3*X.test[i,2]+H.test[i]))
}

# Objective data
X.test <- X.test[,2]
Y.test <- Y.test
H.test <- H.test
A.test <- A.test

# Iterating over xi
xi.vec <- seq(-1,10,by=0.1)
b.AGLM.matrix <- matrix(nrow=length(xi.vec), ncol = 1)
deviance.vec <- numeric(length(xi.vec))

for (i in 1:length(xi.vec)) {
  xi <- xi.vec[i]
  b.AGLM.matrix[i] <- AGLM(xi)
  
  p.hat.test <- exp(X.test*b.AGLM.matrix[i])/(1+exp(X.test*b.AGLM.matrix[i])) # inverse of logit link
  
  special.case1 <- function(Y.test){
    ifelse(Y.test==0, 0, Y.test*log(Y.test/(m*p.hat.test)))
  }
  special.case2 <- function(Y.test){
    ifelse(Y.test==m, 0, (m-Y.test)*log((m-Y.test)/(m-m*p.hat.test)))
  }
  
  r.D.test <- sign(Y.test/m-p.hat.test)*sqrt(2*(special.case1(Y.test)+special.case2(Y.test))) # deviance residuals
  deviance.vec[i] <- 1/n*t(r.D.test)%*%r.D.test
}

# Plot like in ex2 of Rothenhaeusler
plot(xi.vec, deviance.vec, type = "l")
xi.vec[which.min(deviance.vec)]

plot(b.AGLM.matrix, deviance.vec)
