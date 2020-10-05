# Investigating anchor approach using likelihood
#
# Author: Maic Rakitta
# Date: 28.09.20
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
#m <- plogis(1+0.8*X-0.39*H)
#Y <- rbinom(n, 1, plogis(psi*X+log(m/(1-m)))) # if Y also depends on hidden H
Y <- rbinom(n, 1, plogis(psi*X))

# Projection matrix onto column space of anchor A
P.A <- A%*%solve(t(A)%*%A)%*%t(A)

# Inverse of logit
link.inv <- function(X,b){
  return(exp(b*X)/(1+exp(b*X)))
}

##########################################################################
# Fitting AGLM for binomial data with likelihood objective

A_GMM_obj <- function(b){
  
  error <- (Y*(1+exp(b*X))-exp(b*X))/sqrt(exp(b*X))
  
  sum_1 <- sum(Y*X-X*exp(b*X)/(exp(b*X)+1))
  sum_2 <- sum(A*error)
  
  return(sum_1+(gamma-1)*sum_2)
}

gamma <- 0
b.vec <- seq(-10,10,by=0.01)
objec.vec <- numeric(length(b.vec))
for (i in 1:length(b.vec)) {
  b <- b.vec[i]
  objec.vec[i] <- A_GMM_obj(b)
}
plot(b.vec, objec.vec)
index <- which(abs(objec.vec)==min(abs(objec.vec)))
b.vec[index]


A_GMM <- function(gamma){
  A_GMM_obj <- function(b){
    
    error <- (Y*(1+exp(b*X))-exp(b*X))/sqrt(exp(b*X))
    
    sum_1 <- sum(Y*X-X*exp(b*X)/(exp(b*X)+1))
    sum_2 <- sum(A*error)
    
    return(sum_1+(gamma-1)*sum_2)
  }
  
  b.vec <- seq(-10,10,by=0.01)
  objec.vec <- numeric(length(b.vec))
  for (i in 1:length(b.vec)) {
    b <- b.vec[i]
    objec.vec[i] <- A_GMM_obj(b)
  }
  plot(b.vec, objec.vec)
  index <- which(abs(objec.vec)==min(abs(objec.vec)))
  b.vec[index]
}

A_GMM(gamma)


objec.binom <- function(b){
  
  return(-sum(Y*log(link.inv(X,b))+(1-Y)*log(1-link.inv(X,b)))-
           (gamma-1)*sum(P.A%*%Y*log(link.inv(P.A%*%X,b))+(1-P.A%*%Y)*log(1-link.inv(P.A%*%X,b))))
}

# classic MLE
fit <- glm(Y~X-1, family="binomial")
coef(fit)

gamma <- 1
optim <- optimize(objec.binom, interval = c(-20,20))
optim$minimum

A_GMM(gamma)

# Instrumental variables with GLM
library(ivtools)
data <- data.frame(A, X, Y)
fitX.HA <- glm(formula=X~A-1, family="gaussian", data=data) #two-stage estimation
fitY.HX <- glm(formula=Y~X-1, family="binomial", data=data)
fitIV <- ivglm(estmethod="ts", fitX.LZ=fitX.HA, fitY.LX=fitY.HX, data=data,
               ctrl=F)
summary(fitIV)

gamma <- 1000000
optim <- optimize(objec.binom, interval = c(-20,20))
optim$minimum

A_GMM(gamma)

# General gamma
gamma <- 0.5
optim <- optimize(objec.binom, interval = c(-20,20))
optim$minimum

# Looking for true optimal gamma
gamma.vec <- seq(0,100,by=1)
b.vec <- numeric(length(gamma.vec))
for (g in 1:length(gamma.vec)) {
  gamma <- gamma.vec[g]
  b.vec[g] <- optimize(objec.binom, interval = c(-20,20))$minimum
}
plot(gamma.vec,b.vec, ylim=c(0.44,0.5))
abline(h=psi)
