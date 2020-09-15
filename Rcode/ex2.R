# Re-implementation of example 2 in Anchor Regression - Rothenhaeuser 2020
# Author: Maic Rakitta
# Date: 15.09.20

n <- 100

library(extraDistr)
A <- rsign(n)
epsH1 <- rnorm(n)
epsX1 <- rnorm(n)
epsY1 <- rnorm(n)

H1 <- epsH1
X1 <- A+H1+epsX1
Y1 <- X1+2*H1+epsY1



epsH2 <- rnorm(n)
epsX2 <- rnorm(n)
epsY2 <- rnorm(n)

H2 <- epsH2
X2 <- 1.8+H2+epsX2
Y2 <- X2+2*H2+epsY2

#OLS
fit.OLS <- lm(Y1~X1+H1-1)

#IV
fit.2SLS1 <- lm(X1~A-1)
X1.hat <- fitted.values(fit.2SLS1)
fit.IV <- lm(Y1~X1.hat-1)
summary(fit.OLS)
summary(fit.IV)

plot(X1,Y1)
abline(fit.OLS, col="2")
abline(fit.IV, col="3")

#PA
#r2 = residuals(lm(A ~ X1 + H))               
#r2




#Anchor regression
Xnew <- cbind(X1,H1)

gamma <- 20
P.A <- A%*%solve(t(A)%*%A)%*%t(A)
P.Asec <- diag(n)+(sqrt(gamma)-1)*P.A

Y.tilde <- P.Asec%*%Y1
X.tilde <- P.Asec%*%Xnew

fit.AR <- lm(Y.tilde~X.tilde-1)
summary(fit.AR)

#OLS
gamma <- 1
P.A <- A%*%solve(t(A)%*%A)%*%t(A)
P.Asec <- diag(n)+(sqrt(gamma)-1)*P.A

Y.tilde <- P.Asec%*%Y1
X.tilde <- P.Asec%*%Xnew

fit.OLS <- lm(Y.tilde~X.tilde-1)

#IV


#PA
gamma <- 0
P.A <- A%*%solve(t(A)%*%A)%*%t(A)
P.Asec <- diag(n)+(sqrt(gamma)-1)*P.A

Y.tilde <- P.Asec%*%Y1
X.tilde <- P.Asec%*%Xnew

fit.PA <- lm(Y.tilde~X.tilde)

#training data fit
plot(X1,Y1)
abline(fit.OLS, col="2")
abline(fit.IV, col="3")
abline(fit.PA,col="4")
abline(fit.AR,col="5")


#Calculation ex2 figure 1
gamma.vec <- seq(0,2,by=0.001)
Xnew2 <- data.frame(X1=X2,H1=H2)
residuals <- numeric(length(gamma.vec))
b.matrix <- matrix(nrow=length(gamma.vec),ncol=2)

for (i in 1:length(gamma.vec)){
  gamma <- gamma.vec[i]
  P.A <- A%*%solve(t(A)%*%A)%*%t(A)
  P.Asec <- diag(n)+(sqrt(gamma)-1)*P.A
  
  Y.tilde <- P.Asec%*%Y1
  X.tilde <- P.Asec%*%Xnew
  
  fit <- lm(Y.tilde~X.tilde-1)
  
  prediction <- predict(fit, Xnew2)
  residuals[i] <- sum((Y2-prediction)^2)
  b.matrix[i,] <- c(coef(fit))
}

plot(gamma.vec,residuals, type = "l")

prediction <- predict(fit.OLS, Xnew2)
residual <- sum((Y2-prediction)^2)
points(0,residual,col="2")

prediction <- predict(fit.IV, Xnew2)
residual <- sum((Y2-prediction)^2)
points(0,residual,col="3")



plot(b.matrix[,1],residuals, type = "l")

points(b.matrix[which(gamma.vec==1),1],residuals[which(gamma.vec==1)],col="2") #OLS
points(b.matrix[which(gamma.vec==0),1],residuals[which(gamma.vec==0)],col="3") #PA

prediction <- predict(fit.IV, Xnew2)
residual <- sum((Y2-prediction)^2)
points(coef(fit)[1],residual)

prediction <- predict(fit.OLS, Xnew2)
residual <- sum((Y2-prediction)^2)
points(coef(fit)[1],residual)

