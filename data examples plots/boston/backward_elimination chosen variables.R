# GLARE for Boston Housing Data
#
# Date: 04.01.21     Author: Maic Rakitta
###############################################################################

library(glare)
library(mlbench)
library(ggplot2)
library(tidyr)

# Load data
data("BostonHousing2") # see help file why we use nr.2
bostonPrices <- BostonHousing2[, c("medv", "crim", "zn", "indus", "chas", "nox",
                                   "rm", "age", "dis", "rad", "tax", "ptratio",
                                   "b", "lstat", "town")]

# Explore data ----------------------------------------------------------------
plot(bostonPrices)
summary(bostonPrices)

library(Hmisc)
hist.data.frame(BostonHousing2[, c("medv", "crim", "zn", "indus", "nox",
                                   "rm", "age", "dis", "ptratio",
                                   "b", "lstat")])

hist(bostonPrices$medv) # right skewed and >0 suggests log-transform
hist(log(bostonPrices$medv))

hist(bostonPrices$crim)
hist(log(bostonPrices$crim))

hist(bostonPrices$zn)
hist(log(bostonPrices$zn))

hist(bostonPrices$nox)
hist(log(bostonPrices$nox))

hist(bostonPrices$rm)

hist(bostonPrices$age)
#hist(bostonPrices$age^2)

hist(bostonPrices$dis)
hist(log(bostonPrices$dis))

hist(bostonPrices$ptratio)
hist(bostonPrices$ptratio^2)

hist(bostonPrices$b)
hist(bostonPrices$b^2)

hist(bostonPrices$lstat)
hist(log(bostonPrices$lstat))

#suggested tranformations:
#rightskewed->log: medv, crim, dis, nox, zn        and lstat
#leftskwered->square: ptratio                      and b (even higher)
#non: indus, chas, rm, age (^2 results in U), rad, tax, town

# Correlation of continuous variables
cor(bostonPrices[, c("medv", "crim", "zn", "indus", "nox",
                     "rm", "age", "dis", "ptratio",
                     "b", "lstat")], bostonPrices$medv)
# -> lstat&ptratio high negative corr, rm high positive corr

bostonPrices$medv <- log(bostonPrices$medv)
colnames(bostonPrices)[1] <- "log_medv"

# Fit GLARE -------------------------------------------------------------------

# Towns with 15 and more observations
perturbation_towns <- c("Cambridge", "Boston Savin Hill", "Lynn",
                        "Boston Roxbury", "Newton", "Somerville",
                        "Boston South Boston", "Quincy", "Brookline",
                        "Boston East Boston", "Waltham", "Medford",
                        "Boston Dorchester", "Framingham")

fit_list <- list()

for (t in 1:length(perturbation_towns)) {
  
  print(paste(t, "of", length(perturbation_towns)))
  
  pert_town <- perturbation_towns[t]
  
  # Create observed and unobserved environment
  train_set <- bostonPrices[bostonPrices$town != pert_town, ]
  test_set <- bostonPrices[bostonPrices$town == pert_town, ]
  
  # Initialize and run GLARE
  xi_values <- c(0, 1, 5, 10, 50, 100, 10000)
  xi_len <- length(xi_values)
  
  b <- matrix(nrow = xi_len, ncol = 11)
  colnames(b) <- c("Intercept", "crim", "chas", "nox", "rm",
                   "dis", "rad", "tax", "ptratio", "b", "lstat")
  b_se <- matrix(nrow = xi_len, ncol = 11)
  colnames(b_se) <- colnames(b)
  
  predictions <- matrix(nrow = xi_len, ncol = nrow(test_set))
  resid_pert <- matrix(nrow = xi_len, ncol = nrow(test_set))
  RMSE <- numeric(xi_len)
  
  logLik_indiv <- matrix(nrow = xi_len, ncol = nrow(test_set))
  logLik <- numeric(xi_len)
  
  pb_xi <- txtProgressBar(min = 0, max = xi_len, style = 3)
  for (i in 1:xi_len) {
    
    xi <- xi_values[i]
    
    fit_temp <- glare(formula = log_medv ~ crim + chas + nox + rm +
                        dis + rad + tax + ptratio + b + lstat,
                      A_formula = ~ town, data = train_set, xi = xi,
                      family = gaussian, type = "pearson")
    
    # Coefficients
    b[i, ] <- as.numeric(coef(fit_temp))
    b_se[i, ] <- fit_temp$coef_se
    
    # Predictions, Residuals and MSE
    predictions[i, ] <- predict(fit_temp, type = "response", newdata = test_set)
    resid_pert[i, ] <- test_set$log_medv - predictions[i, ]
    RMSE[i] <- sqrt(mean((resid_pert[i, ])^2))
    
    # Log-likelihood
    logLik_indiv[i, ] <- as.numeric(logLik(fit_temp, newdata = test_set,
                                           indiv = TRUE))
    logLik[i] <- as.numeric(quantile(logLik_indiv[i, ], 0.9))
    
    Sys.sleep(0.1)
    setTxtProgressBar(pb_xi, i)
  }
  close(pb_xi)
  
  fit_list[[t]] <- list(b, b_se, predictions, RMSE, logLik_indiv, logLik)
}

path_name <- "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/data_examples/bostonHousing/"
dir.create(paste(path_name, Sys.Date(), sep = ""))
save(fit_list, file = paste(paste(path_name, Sys.Date(), sep = ""), "/fit_list_no_age_indus_zn.Rdata", sep =""))

# Plots -----------------------------------------------------------------------

# Prepare data
xi_values <- c(0, 1, 5, 10, 50, 100, 10000)
xi_len <- length(xi_values)
load("C:/Users/maicr/Desktop/Github/MScAnchor/data sets/data_examples/bostonHousing/2021-01-12/fit_list_no_age_indus_zn.Rdata")

gg_data <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(gg_data) <- c("test_town", "xi", "RMSE")

for (t in 1:length(perturbation_towns)) {
  
  test_town <- perturbation_towns[t]
  
  for (i in 1:xi_len) {
    
    xi <- xi_values[i]
    
    RMSE <- fit_list[[t]][[4]][i]
    
    data_temp <- data.frame(test_town, xi, RMSE)
    gg_data <- rbind(gg_data, data_temp)
  }
}

# Plot all
gg_data$xi <- factor(gg_data$xi, levels = c("0", "1", "5", "10", "50", "100", "10000"))

gg_data$test_town <- factor(gg_data$test_town)
names.arg <- levels(gg_data$test_town)
levels(gg_data$test_town) <- abbreviate(names.arg, minlength = 2, use.classes = TRUE,
                                        dot = FALSE, strict = FALSE,
                                        method = c("left.kept", "both.sides"), named = TRUE)

ggplot(gg_data, aes(y = RMSE, x = test_town, color = xi, group = xi)) +
  
  geom_line() +
  
  labs(color = expression(xi)) +
  ylab("RMSE") +
  xlab("Perturbation town")

# Plot 0, 5 , 50, 10000
gg_data2 <- gg_data[gg_data$xi %in% c(0,5,50,10000), ]

ggplot(gg_data2, aes(y = RMSE, x = test_town, color = xi, group = xi)) +
  
  geom_line() +
  
  labs(color = expression(xi)) +
  ylab("RMSE") +
  xlab("Perturbation town")

# citys with 50 as best:
# somerville, Boston East Boston, Medford, Boston Dorchester, Newton

# citys with xi_big as best:
# Boston South Boston, Brookline

# citys with 5 as best:
# Chambridge

# rest best with 0









# Create observed and unobserved environment
train_set <- bostonPrices[bostonPrices$town != "Quincy", ]
test_set <- bostonPrices[bostonPrices$town == "Quincy", ]

nrow(train_set)
nrow(test_set)

# Initialize and run GLARE
xi_values <- seq(0, 100, by = 1)

#xi_values <- c(0, 5, 10, 20, 50, 100, 200, 400, 700, 1000)

xi_len <- length(xi_values)

b <- matrix(nrow = xi_len, ncol = 11)
colnames(b) <- c("Intercept", "crim", "chas", "nox", "rm",
                 "dis", "rad", "tax", "ptratio", "b", "lstat")
b_se <- matrix(nrow = xi_len, ncol = 11)
colnames(b_se) <- colnames(b)

predictions <- matrix(nrow = xi_len, ncol = nrow(test_set))
resid_pert <- matrix(nrow = xi_len, ncol = nrow(test_set))
RMSE <- numeric(xi_len)

logLik_indiv <- matrix(nrow = xi_len, ncol = nrow(test_set))
logLik <- numeric(xi_len)

pb <- txtProgressBar(min = 0, max = xi_len, style = 3)
for (i in 1:xi_len) {
  
  xi <- xi_values[i]
  
  fit_temp <- glare(formula = log_medv ~ crim + chas + nox + rm +
                      dis + rad + tax + ptratio + b + lstat,
                    A_formula = ~ town, data = train_set, xi = xi,
                    family = gaussian, type = "pearson")
  # 
  # fit_temp <- glare(formula = medv ~ crim + zn + indus + chas + nox + rm + age +
  #                     dis + rad + tax + ptratio + b + lstat,
  #                   A_formula = ~ town, data = train_set, xi = xi,
  #                   family = gaussian, type = "pearson")
  
  # resid_pert <- test_set$medv - predict(fit_temp,
  #                                            type = "response",
  #                                            newdata = test_set)
  
  # Coefficients
  b[i, ] <- as.numeric(coef(fit_temp))
  b_se[i, ] <- fit_temp$coef_se
  
  # Predictions, Residuals and MSE
  predictions[i, ] <- predict(fit_temp, type = "response", newdata = test_set)
  resid_pert[i, ] <- test_set$log_medv - predictions[i, ]
  RMSE[i] <- sqrt(mean((resid_pert[i, ])^2))
  
  # Log-likelihood
  logLik_indiv[i, ] <- as.numeric(logLik(fit_temp, newdata = test_set,
                                         indiv = TRUE))
  logLik[i] <- as.numeric(quantile(logLik_indiv[i, ], 0.9))
  
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
}
close(pb)

list <- list(b, b_se, predictions, RMSE, logLik_indiv, logLik)
path_name <- "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/data_examples/bostonHousing/"
dir.create(paste(path_name, Sys.Date(), sep = ""))
save(list, file = paste(paste(path_name, Sys.Date(), sep = ""), "/list_Medford_all.Rdata", sep =""))

#load("C:/Users/maicr/Desktop/Github/MScAnchor/data sets/data_examples/bostonHousing/2021-01-11/list_Medford_all.Rdata")

b <- list[[1]]


# Fit GLARE for xi big
xi_big <- 10000
fit_big <- glare(formula = log_medv ~ crim + zn + indus + chas + nox + rm +
                   age + dis + rad + tax + ptratio + b + lstat,
                 A_formula = ~ town, data = train_set, xi = xi_big,
                 family = gaussian, type = "pearson")
b_big <- as.numeric(coef(fit_big))
b_se_big <- fit_big$coef_se
predictions_big <- predict(fit_big, type = "response", newdata = test_set)
resid_pert_big <- test_set$log_medv - predictions_big
RMSE_big <- sqrt(mean((resid_pert_big)^2))

b_big
RMSE_big



# Plots -----------------------------------------------------------------------

# Parameter plot
gg_data_b <- as.data.frame(cbind(xi_values, b))

library(dplyr)
df_b <- gg_data_b %>%
  select(xi_values, Intercept, crim, zn, indus, chas, nox, rm,
         age, dis, rad, tax, ptratio, b, lstat) %>%
  gather(key = "Variable", value = "value", -xi_values)
head(df_b)

ggplot(df_b, aes(x = xi_values, y = value)) + 
  geom_line(aes(color = Variable)) +
  ylab("Parameters") +
  xlab(expression(xi)) +
  geom_vline(xintercept = 0, linetype = "dashed")

# 0.9-quantile of log-likelihood plot
gg_data_log <- data.frame(logLik_pert = logLik, xi = xi_values)

ggplot(gg_data_log, aes(y = logLik_pert, x = xi)) +
  geom_line() +
  ylab("0.9-quantile of log-likelihood") +
  xlab(expression(xi)) +
  geom_vline(xintercept = 0, linetype = "dashed")

# MSE
gg_data <- data.frame(RMSE = RMSE, xi = xi_values)

ggplot(gg_data, aes(y = RMSE, x = xi)) +
  geom_line() +
  xlab(expression(xi))


gg_data2 <- gg_data[gg_data$xi <101,]

ggplot(gg_data2, aes(y = RMSE, x = xi)) +
  geom_line() +
  xlab(expression(xi))




# Single comparison -----------------------------------------------------------

# Create observed and unobserved environment
train_set <- bostonPrices[bostonPrices$town != "Medford", ]
test_set <- bostonPrices[bostonPrices$town == "Medford", ]

# Fit GLARE for xi fixed
xi <- 5
fit_glare <- glare(formula = log_medv ~ crim + zn + indus + chas + nox + rm +
                     age + dis + rad + tax + ptratio + b + lstat,
                   A_formula = ~ town, data = train_set, xi = xi,
                   family = gaussian, type = "pearson")
summary(fit_glare)


b_big <- as.numeric(coef(fit_big))
b_se_big <- fit_big$coef_se
predictions_big <- predict(fit_big, type = "response", newdata = test_set)
resid_pert_big <- test_set$log_medv - predictions_big
RMSE_big <- sqrt(mean((resid_pert_big)^2))

b_big
RMSE_big

# Fit GLM
fit_glm <- glm(formula = log_medv ~ crim + zn + indus + chas + nox + rm +
                 age + dis + rad + tax + ptratio + b + lstat,
               data = train_set, family = gaussian)
summary(fit_glm)

# no age
fit_glm <- glm(formula = log_medv ~ crim + zn + indus + chas + nox + rm +
                 dis + rad + tax + ptratio + b + lstat,
               data = train_set, family = gaussian)
summary(fit_glm)

# no indus
fit_glm <- glm(formula = log_medv ~ crim + zn + chas + nox + rm +
                 dis + rad + tax + ptratio + b + lstat,
               data = train_set, family = gaussian)
summary(fit_glm)

# no zn
fit_glm <- glm(formula = log_medv ~ crim + chas + nox + rm +
                 dis + rad + tax + ptratio + b + lstat,
               data = train_set, family = gaussian)
summary(fit_glm)
