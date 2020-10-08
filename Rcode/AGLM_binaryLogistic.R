# Anchor regression using GLM approach and CCXR package
#
# Author: Maic Rakitta
# Date: 07.10.20
##########################################################################
set.seed(1)

n <- 1000 # number of samples from unpertubed and pertubed distribution

# initialize data
library(extraDistr) # for rademacher distribution
A <- rsign(n) # anchor variable

epsH <- rnorm(n)
epsX <- rnorm(n)
H <- epsH # hidden confounder
X <- A+H+epsX # dependend variable

# Binomial Data
psi <- 0.5 # true b
#o <- plogis(1+0.8*X-0.39*H)
#Y <- rbinom(n, m, plogis(psi*X+log(o/(1-o)))) # if Y also depends on hidden H

m <- 5
Y <- rbinom(n, m, plogis(psi*X))

# Projection matrix onto column space of anchor A
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

