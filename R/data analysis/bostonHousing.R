# Using GLARE for Boston Housing Data
#
# Date: 04.01.21     Author: Maic Rakitta
###############################################################################

library(glare)
library(mlbench)

data("BostonHousing2")
summary(BostonHousing2)

# Create observed and unobserved environment
train_set <- BostonHousing2[BostonHousing2$town != "Cambridge", ]
test_set <- BostonHousing2[BostonHousing2$town == "Cambridge", ]

# Initialize and run GLARE
xi_values <- seq(0, 10, by = 0.1)
xi_len <- length(xi_values)

b <- matrix(nrow = xi_len, ncol = 16)
b_se <- matrix(nrow = xi_len, ncol = 16)
MSE_pert <- numeric(xi_len)

pb <- txtProgressBar(min = 0, max = xi_len, style = 3)
for (i in 1:xi_len) {
  
  xi <- xi_values[i]
  
  fit_temp <- glare(formula = medv ~ crim + zn + indus + chas + nox + rm + age +
                      dis + rad + tax + ptratio + b + lstat + lon + lat,
                    A_formula = ~ town, data = train_set, xi = xi,
                    family = gaussian, type = "pearson")
  
  b[i, ] <- as.numeric(coef(fit_temp))
  b_se[i, ] <- fit_temp$coef_se
  
  resid_pert <- test_set$medv - predict(fit_temp,
                                        type = "response",
                                        newdata = test_set)
  MSE_pert[i] <- mean((resid_pert)^2)
  
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
}
close(pb)


# Plots
gg_data <- data.frame(MSE_pert = MSE_pert, xi = xi_values)

ggplot(gg_data, aes(y = MSE_pert, x = xi)) +
  geom_line() +
  xlab(expression(xi))
