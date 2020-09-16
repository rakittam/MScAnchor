# Re-implementation of example 2 in Anchor Regression - Rothenhaeuser 2020
# Author: Maic Rakitta
# Date: 15.09.20
##########################################################################
set.seed(1)

n <- 1000

library(extraDistr)

# initialize training data
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

##########################################################################
#Anchor regression

anchor.regression <- function(X, Y, A, gamma, n){
  
  P.A <- A%*%solve(t(A)%*%A)%*%t(A)
  W <- diag(n)-(1-sqrt(gamma))*P.A
  
  Y.tilde <- W%*%Y
  X.tilde <- W%*%X
  
  fit.AR <- fit.AR <- lm(Y.tilde~X.tilde-1)
}

fit.AR <- anchor.regression(X.train, Y.train, A, 2, n)
summary(fit.AR)


#PA
fit.PA <- anchor.regression(X.train, Y.train, A, 0, n)
b.PA <- coef(fit.PA)

#OLS
fit.OLS <- anchor.regression(X.train, Y.train, A, 1, n)
b.OLS <- coef(fit.OLS)

#IV
fit.2SLS.step <- lm(X.train~A-1)
X.train.hat <- fitted.values(fit.2SLS.step)
fit.IV <- lm(Y.train~X.train.hat-1)

#P.A <- A%*%solve(t(A)%*%A)%*%t(A) # manually
#X.train.tilde <- P.A%*%X.train
#b.IV <- solve(t(X.train.tilde)%*%X.train.tilde)%*%t(X.train.tilde)%*%Y.train
b.IV <- coef(fit.IV)

#training data fit
plot(X.train,Y.train)
abline(fit.OLS, col="2")
abline(fit.IV, col="3")
abline(fit.PA,col="4")
legend(0, -5, legend=c("OLS", "IV", "PA"),
       col=c(2, 3,4), lty=1, cex=0.8)

#test data fit
plot(X.test,Y.test)
abline(fit.OLS, col="2")
abline(fit.IV, col="3")
abline(fit.PA,col="4")
legend(2, -3, legend=c("OLS", "IV", "PA"),
       col=c(2, 3,4), lty=1, cex=0.8)

# Training MSE
MSE.PA <- mean((Y.test - t(X.test)*b.PA) ^ 2)
MSE.OLS <- mean((Y.test - t(X.test)*b.OLS) ^ 2)
MSE.IV <- mean((Y.test - t(X.test)*b.IV) ^ 2)

#Calculation ex2 figure 1
#gamma.vec <- seq(0,10,by=0.01)
gamma.vec <- seq(0,100,by=0.1)
b.vec <- numeric(length(gamma.vec))
MSE.vec <- numeric(length(gamma.vec))

for (i in 1:length(gamma.vec)){
  
  gamma <- gamma.vec[i]
  
  fit <- anchor.regression(X.train, Y.train, A, gamma, n)
  
  MSE.vec[i] <- mean((Y.test - t(X.test)*coef(fit)) ^ 2)
  b.vec[i] <- coef(fit)
}

# gamma limit for comparison to IV
gamma <- 10000
fit <- anchor.regression(X.train, Y.train, A, gamma, n)
MSE.limit <- mean((Y.test - t(X.test)*coef(fit)) ^ 2)
b.limit <- coef(fit)

b.vec.limit <- c(b.vec,b.limit)
MSE.vec.limit <- c(MSE.vec,MSE.limit)
gamma.vec.limit <- c(gamma.vec,gamma)

plot(b.vec.limit, MSE.vec.limit, type = "l", xlim = c(0.95,2.05))
points(b.OLS, MSE.OLS,col="2", pch=16)
points(b.IV, MSE.IV,col="3", pch=16)
points(b.PA, MSE.PA,col="4", pch=16)
legend(1, 6.5, legend=c("OLS", "IV", "PA"),
       col=c(2, 3,4), cex=0.8,pch=16)

plot(gamma.vec, MSE.vec, type = "l")
points(1, MSE.OLS,col="2", pch=16)
points(0, MSE.PA,col="4", pch=16)
abline(h=MSE.IV, col=3)
legend(60, 5.5, legend=c("OLS", "IV MSE" ,"PA"),
       col=c(2, 3, 4), cex=0.8,pch=16)

plot(gamma.vec.limit, MSE.vec.limit, type = "l")
points(1, MSE.OLS,col="2", pch=16)
points(0, MSE.PA,col="4", pch=16)
abline(h=MSE.IV)
abline(h=MSE.IV, col=3)
legend(2000, 5.5, legend=c("OLS", "IV MSE" ,"PA"),
       col=c(2, 3, 4), cex=0.8,pch=16)