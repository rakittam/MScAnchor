# Anchor regression function
#
# Input:  X Covectors
#         Y Responde
#         A Anchor variables
#         gamma hyperparameter
#
# Output: anchor regression fit
#
# Author: Maic Rakitta
# Date: 23.09.20
##########################################################################

anchor.regression <- function(X, Y, A, gamma){
  
  n <- length(Y)
  
  P.A <- A%*%solve(t(A)%*%A)%*%t(A)
  W <- diag(n)-(1-sqrt(gamma))*P.A
  
  Y.tilde <- W%*%Y
  X.tilde <- W%*%X
  
  fit <- lm(Y.tilde~X.tilde-1)
}
