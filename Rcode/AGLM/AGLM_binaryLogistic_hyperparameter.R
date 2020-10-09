# Using test data
# Anchor regression using GLM approach and optim
# Binary logistic regression
#
# Author: Maic Rakitta
# Date: 09.10.20
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



# initialize training data
A.test <- matrix(nrow = n, ncol = 2)
H.test <- matrix(nrow = n, ncol = 1)
X.test <- matrix(nrow = n, ncol = 10)
Y.test <- matrix(nrow = n, ncol = 1)

for (i in 1:n) {
  
  A.test[i,] <- rsign(n=2)
  
  epsH.test <- rnorm(n=1, mean=0, sd=1)
  H.test[i] <- epsH.test
  
  epsX.test <- rnorm(n=10, mean=0, sd=1)
  X.test[i,] <- 4+4+H.test[i]+epsX.test
  
  Y.test[i] <- rbinom(n=1, size=m, plogis(3*X.test[i,2]+3*X.test[i,3]+H.test[i]-2*A.test[i,1]))
}

# Objective data
X.test <- X.test[,2:3]
Y.test <- Y.test
H.test <- H.test
A.test <- A.test

##########################################################################
# AGLM for Logistic regression data using optim as optimizer
AGLM <- function(xi){
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
  
  start.val <- c(1,1)
  ans2 <- optim(par=start.val, fn=objective, hessian = TRUE)
  return(ans2$par)
}

##########################################################################
# Playground
gamma <- 2 # Set gamma and xi
xi <- gamma-1
AGLM(xi)

# Iterating over different hyper parameters
xi.vec <- -1:10
b.AGLM.matrix <- matrix(nrow=length(xi.vec), ncol = 2)
for (i in 1:length(xi.vec)) {
  xi <- xi.vec[i]
  b.AGLM.matrix[i,] <- AGLM(xi)
}

plot(xi.vec, b.AGLM.matrix[,1], type = "l")
lines(xi.vec, b.AGLM.matrix[,2])
abline(h=3)
