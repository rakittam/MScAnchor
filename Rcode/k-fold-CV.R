

##### In construction #####


# template k-fold CV


##########################################################################
# Optimal gamma calculation with cross validation
data <- data.frame(Y=Y.train,X=X.train, A=A)

k <- 10
folds <- sample(1:k, nrow(data), replace=T)

CV.MSE <- numeric(k)
for (out in 1:k) {
  
  # Split the data into training and test set
  train <- data[folds!=out,]
  test <- data[folds==out,]
  
  # Build the model
  model <- anchor.regression(train$X, train$Y, train$A, gamma, nrow(train))
  
  # Make predictions and compute test MSE
  predictions <- test$X * model$coefficients
  CV.MSE[out] <- mean((test$Y - predictions) ^ 2)
  
}
mean(CV.MSE)