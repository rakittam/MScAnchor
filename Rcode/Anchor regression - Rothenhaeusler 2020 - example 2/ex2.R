# Re-implementation of example 2 in Anchor Regression - Rothenhaeuser 2020
#
# Author: Maic Rakitta
# Date: 15.09.20
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

##########################################################################
# Anchor regression
anchor.regression <- function(X, Y, A, gamma, n){
  
  P.A <- A%*%solve(t(A)%*%A)%*%t(A)
  W <- diag(n)-(1-sqrt(gamma))*P.A
  
  Y.tilde <- W%*%Y
  X.tilde <- W%*%X
  
  fit <- lm(Y.tilde~X.tilde-1)
}

# PA = AR for gamma = 0
fit.PA <- anchor.regression(X.train, Y.train, A, 0, n)
b.PA <- coef(fit.PA)

# OLS = AR for gamma = 1
fit.OLS <- anchor.regression(X.train, Y.train, A, 1, n)
b.OLS <- coef(fit.OLS)

# IV - 2SLS
fit.2SLS.step <- lm(X.train~A-1)
X.train.hat <- fitted.values(fit.2SLS.step)
fit.IV <- lm(Y.train~X.train.hat-1)

#P.A <- A%*%solve(t(A)%*%A)%*%t(A) # manually
#X.train.tilde <- P.A%*%X.train
#b.IV <- solve(t(X.train.tilde)%*%X.train.tilde)%*%t(X.train.tilde)%*%Y.train
b.IV <- coef(fit.IV)

# training data plot
plot(X.train,Y.train)
abline(fit.OLS, col="2")
abline(fit.IV, col="3")
abline(fit.PA,col="4")
legend(0, -5, legend=c("OLS", "IV", "PA"),
       col=c(2, 3,4), lty=1, cex=0.8)

# test data plot
plot(X.test,Y.test)
abline(fit.OLS, col="2")
abline(fit.IV, col="3")
abline(fit.PA,col="4")
legend(2, -3, legend=c("OLS", "IV", "PA"),
       col=c(2, 3,4), lty=1, cex=0.8)

##########################################################################
# Calculation and plot of ex2 figure 1 in Rothenhaeusler 2020
gamma.vec <- seq(0,100,by=0.1)
b.vec <- numeric(length(gamma.vec))
MSE.vec <- numeric(length(gamma.vec))

for (i in 1:length(gamma.vec)){
  
  gamma <- gamma.vec[i]
  
  fit <- anchor.regression(X.train, Y.train, A, gamma, n)
  
  MSE.vec[i] <- mean((Y.test - t(X.test)*coef(fit)) ^ 2)
  b.vec[i] <- coef(fit)
}

# gamma big enough to compare with IV
gamma <- 10000
fit <- anchor.regression(X.train, Y.train, A, gamma, n)
MSE.limit <- mean((Y.test - t(X.test)*coef(fit)) ^ 2)
b.limit <- coef(fit)

b.vec.limit <- c(b.vec,b.limit)
MSE.vec.limit <- c(MSE.vec,MSE.limit)
gamma.vec.limit <- c(gamma.vec,gamma)

# test MSE for specific cases
MSE.PA <- mean((Y.test - t(X.test)*b.PA) ^ 2)
MSE.OLS <- mean((Y.test - t(X.test)*b.OLS) ^ 2)
MSE.IV <- mean((Y.test - t(X.test)*b.IV) ^ 2)

# plots
plot(gamma.vec.limit, MSE.vec.limit, type = "l")
points(1, MSE.OLS,col="2", pch=16)
points(0, MSE.PA,col="4", pch=16)
abline(h=MSE.IV)
abline(h=MSE.IV, col=3)
legend(2000, 5.5, legend=c("OLS", "IV MSE" ,"PA"),
       col=c(2, 3, 4), cex=0.8,pch=16)

plot(gamma.vec, MSE.vec, type = "l")
points(1, MSE.OLS,col="2", pch=16)
points(0, MSE.PA,col="4", pch=16)
abline(h=MSE.IV, col=3)
legend(60, 5.5, legend=c("OLS", "IV MSE" ,"PA"),
       col=c(2, 3, 4), cex=0.8,pch=16)

plot(b.vec.limit, MSE.vec.limit, type = "l", xlim = c(0.95,2.05))
points(b.OLS, MSE.OLS,col="2", pch=16)
points(b.IV, MSE.IV,col="3", pch=16)
points(b.PA, MSE.PA,col="4", pch=16)
legend(1, 6.5, legend=c("OLS", "IV", "PA"),
       col=c(2, 3,4), cex=0.8,pch=16)

# final plot: b and gamma vs test MSE
P.A <- A%*%solve(t(A)%*%A)%*%t(A) # first calulate coresponding gammas
X <- X.train
Y <- Y.train
b.gamma <- c(1,1.25,1.5,1.75,2)
b.gamma <- seq(1,2,by=1/9)
gamma.axis <- numeric(length(b.gamma))
for (i in 1:length(b.gamma)) {
  
  b <- b.gamma[i]
  
  gamma.axis[i] <- (-2*t(X)%*%X*b+t(X)%*%Y+t(Y)%*%X) * # solved by derivation of b optictive, set to 0 and solve for gamma
    solve(2*t(X)%*%P.A%*%X*b-t(X)%*%P.A%*%Y-t(Y)%*%P.A%*%X)+1
  
}

gamma.axis
c("inf",17.6987,7.894467, 4.542962, 2.851117, 1.8308, 1.148405, 0.659905, 0.2929409, 0.007174545)

library(ggplot2)
data <- data.frame(b=b.vec.limit, gamma = gamma.vec.limit, testMSE = MSE.vec.limit)
ggplot(data, aes(y=testMSE)) +
  
  geom_line( aes(x=b))+
  
  annotate("point", x = b.OLS, y = MSE.OLS, colour = 2)+
  annotate("text", x = b.OLS-0.02, y = MSE.OLS+0.1, label = "OLS")+
  annotate("point", x = b.IV, y = MSE.IV, colour = 3)+
  annotate("text", x = b.IV, y = MSE.IV+0.1, label = "IV")+
  annotate("point", x = b.PA, y = MSE.PA, colour = 4)+
  annotate("text", x = b.PA-0.04, y = MSE.PA, label = "PA")+
  
  scale_x_continuous(
    
    # Features of the first axis
    name = "b",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans= ~.,
                        name="gamma",
                        breaks=c(1.000000, 1.111111, 1.222222, 1.333333, 1.444444, 1.555556, 1.666667, 1.777778, 1.888889, 2.000000),
                        labels=c("inf",17.7,7.9, 4.5, 2.9, 1.8, 1.1, 0.7, 0.3, 0.0))
  )

# unobserved optimal gamma
b.opt <- b.vec.limit[which(MSE.vec.limit==min(MSE.vec.limit))]
gamma.opt <- (-2*t(X)%*%X*b.opt+t(X)%*%Y+t(Y)%*%X) * # solved by derivation of b optictive, set to 0 and solve for gamma
  solve(2*t(X)%*%P.A%*%X*b.opt-t(X)%*%P.A%*%Y-t(Y)%*%P.A%*%X)+1

gamma.opt
