# GLARE for Boston Housing Data
#
# Date: 08.02.21     Author: Maic Rakitta
###############################################################################

library(glare)
library(mlbench)
library(ggplot2)
library(tidyr)
library(ggpubr) # for theme and ggarrange

theme_set(theme_bw())

# Load data
data("BostonHousing2")
bostonPrices <- BostonHousing2[, c("cmedv", "crim", "zn", "indus", "nox",
                                   "rm", "age", "lstat", "town")]

# Explore data ----------------------------------------------------------------
plot(bostonPrices)
summary(bostonPrices)

table(bostonPrices$town)

library(Hmisc)
hist.data.frame(BostonHousing2[, c("cmedv", "crim", "zn", "indus", "nox",
                                   "rm", "age", "lstat")])

hist(bostonPrices$cmedv) # right skewed and >0 suggests log-transform
hist(log(bostonPrices$cmedv))

hist(bostonPrices$crim)
hist(log(bostonPrices$crim))

hist(bostonPrices$zn)
hist(log(bostonPrices$zn))

hist(bostonPrices$rm)

hist(bostonPrices$age)
#hist(bostonPrices$age^2)

hist(bostonPrices$lstat)
hist(log(bostonPrices$lstat))

#suggested tranformations:
#rightskewed->log: cmedv, crim, nox, zn        and lstat
#non: indus, rm, age (^2 results in U), town

# Correlation of continuous variables
cor(bostonPrices[, c("cmedv", "crim", "zn", "indus", "nox",
                     "rm", "age", "lstat")], bostonPrices$cmedv)
# -> lstat high negative corr, rm high positive corr

bostonPrices$cmedv <- log(bostonPrices$cmedv)
colnames(bostonPrices)[1] <- "log_cmedv"

# Fit GLARE -------------------------------------------------------------------

# Initialize
xi_values <- c(0, 1, 5, 10, 50, 100, 10000)
xi_len <- length(xi_values)

bostonPrices$town <- factor(bostonPrices$town)
names.arg <- levels(bostonPrices$town)
levels(bostonPrices$town) <- abbreviate(names.arg, minlength = 2, use.classes = TRUE,
                                        dot = FALSE, strict = FALSE,
                                        method = c("left.kept", "both.sides"), named = TRUE)

perturbation_towns <- names(table(bostonPrices$town))

fit_list <- list()

for (t in 1:length(perturbation_towns)) {
  
  print(paste(t, "of", length(perturbation_towns)))
  
  pert_town <- perturbation_towns[t]
  
  # Create observed and unobserved environment
  train_set <- bostonPrices[bostonPrices$town != pert_town, ]
  test_set <- bostonPrices[bostonPrices$town == pert_town, ]

  # Run GLAREs
  b <- matrix(nrow = xi_len, ncol = 8)
  colnames(b) <- c("Intercept", "crim", "zn", "indus",
                   "nox", "rm", "age","lstat")
  b_se <- matrix(nrow = xi_len, ncol = 8)
  colnames(b_se) <- colnames(b)
  
  predictions <- matrix(nrow = xi_len, ncol = nrow(test_set))
  resid_pert <- matrix(nrow = xi_len, ncol = nrow(test_set))
  RMSE <- numeric(xi_len)
  
  pb_xi <- txtProgressBar(min = 0, max = xi_len, style = 3)
  for (i in 1:xi_len) {
    
    xi <- xi_values[i]
    
    fit_temp <- glare(formula = log_cmedv ~ crim + zn + indus + nox + rm +
                        age + lstat,
                      A_formula = ~ town, data = train_set, xi = xi,
                      family = gaussian, type = "pearson")
    
    # Coefficients
    b[i, ] <- as.numeric(coef(fit_temp))
    b_se[i, ] <- fit_temp$coef_se
    
    # Predictions, Residuals and MSE
    predictions[i, ] <- predict(fit_temp, type = "response", newdata = test_set)
    resid_pert[i, ] <- test_set$log_cmedv - predictions[i, ]
    RMSE[i] <- sqrt(mean((resid_pert[i, ])^2))
    
    Sys.sleep(0.1)
    setTxtProgressBar(pb_xi, i)
  }
  close(pb_xi)
  
  fit_list[[t]] <- list(b, b_se, predictions, RMSE)
}

# path_name <- "./data sets/data_examples/bostonHousing/"
# dir.create(paste(path_name, Sys.Date(), sep = ""))
# save(fit_list, file = paste(paste(path_name, Sys.Date(), sep = ""), "/fit_list.Rdata", sep =""))

# Plots -----------------------------------------------------------------------

# Prepare data
# xi_values <- c(0, 1, 5, 10, 50, 100, 10000)
# xi_len <- length(xi_values)
# load("./data sets/data_examples/bostonHousing/2021-02-08/fit_list.Rdata")

gg_data <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(gg_data) <- c("test_town", "xi", "RMSE", "b")

b_matrix <- matrix(0, nrow = 0, ncol = 10)
colnames(b_matrix) <- c("xi", "Town", "Intercept", "crim", "zn", "indus",
                 "nox", "rm", "age","lstat")
b_matrix <- data.frame(b_matrix)
b_matrix <- transform(b_matrix, xi = as.character(xi))
b_matrix <- transform(b_matrix, Town = as.character(Town))

for (t in 1:length(perturbation_towns)) {
  
  test_town <- perturbation_towns[t]
  
  for (i in 1:xi_len) {
    
    xi <- xi_values[i]
    
    RMSE <- fit_list[[t]][[4]][i]
    b <- fit_list[[t]][[1]][i, ]
  
    data_temp <- data.frame(test_town, xi, RMSE, b)
    gg_data <- rbind(gg_data, data_temp)
    
    b_frame <- data.frame(t(b))
    b_frame <- cbind(xi, test_town, b_frame)
    b_matrix <- rbind(b_matrix, b_frame)
  }
}

# Plot all
gg_data$xi <- factor(gg_data$xi, levels = c(0, 1, 5, 10, 50, 100, 10000))

ggplot(gg_data, aes(y = RMSE, x = test_town, color = xi, group = xi)) +
  
  geom_line() +
  
  labs(color = expression(xi)) +
  ylab("RMSE") +
  xlab("Perturbation town") +
  scale_color_viridis_d(option = "C", end = 0.95)

# Plot 0, 5 , 50, 10000
gg_data2 <- gg_data[gg_data$xi %in% c(0,5,50,10000), ]

gg <- ggplot(gg_data2, aes(y = RMSE, x = test_town, color = xi, group = xi)) +
  
  geom_line() +
  
  labs(color = expression(xi)) +
  ylab("RMSE") +
  xlab("Perturbation town") +
  scale_color_viridis_d(option = "C", end = 0.95)

ggsave("all.pdf", gg, width = 30, height = 9)

# Only certain towns
gg_data_towns <- gg_data2[gg_data2$test_town %in% c("Ar", "As", "Bd",
                                                    "Bl", "Bv", "BA",
                                                    "BBB", "BBH", "BC",
                                                    "BstnDr", "BstnDw", "BEB",
                                                    "BFH", "BHP", "BM",
                                                    "BNE"), ]

ggplot(gg_data_towns, aes(y = RMSE, x = test_town, color = xi, group = xi)) +
  
  geom_line() +
  
  labs(color = expression(xi)) +
  ylab("RMSE") +
  xlab("Perturbation town") +
  scale_color_viridis_d(option = "C", end = 0.95)

# Upper towns
upper_RMSE <- sort(gg_data$RMSE, decreasing = TRUE)[1:500]
upper_index <- which(gg_data$RMSE %in% upper_RMSE)
upper_table <- table(gg_data$test_town[upper_index])
upper_table_t <- upper_table[upper_table != 0]
upper_towns <- names(upper_table_t)

gg_data_upper <- gg_data2[gg_data2$test_town %in% upper_towns, ]

gg_upper <- ggplot(gg_data_upper, aes(y = RMSE, x = test_town, color = xi, group = xi)) +
  
  geom_line() +
  
  labs(color = expression(xi)) +
  ylab("RMSE") +
  xlab("Perturbation town") +
  scale_color_viridis_d(option = "C", end = 0.95)

ggsave("upper.pdf", gg_upper, width = 9, height = 9)

# Towns with more then 10 obs
table_10 <- table(bostonPrices$town)[table(bostonPrices$town) >= 10]
towns_10 <- names(table_10)

gg_data_10 <- gg_data2[gg_data2$test_town %in% towns_10, ]

gg_10 <- ggplot(gg_data_10, aes(y = RMSE, x = test_town, color = xi, group = xi)) +
  
  geom_line() +
  
  labs(color = expression(xi)) +
  ylab("RMSE") +
  xlab("Perturbation town") +
  scale_color_viridis_d(option = "C", end = 0.95)

gg_10

ggsave("10obs.pdf", gg_10, width = 9, height = 9)

# -----------------------------------------------------------------------------
# Parameter plot
gg_b_data <- b_matrix
gg_b_data$xi <- factor(gg_b_data$xi, levels = c(0, 1, 5, 10, 50, 100, 10000))
gg_b_data$test_town <- factor(gg_b_data$test_town)
names.arg <- levels(gg_b_data$test_town)
levels(gg_b_data$test_town) <- abbreviate(names.arg, minlength = 2, use.classes = TRUE,
                                        dot = FALSE, strict = FALSE,
                                        method = c("left.kept", "both.sides"), named = TRUE)


gg_b_data2 <- gg_b_data[gg_b_data$xi %in% c(0,5,50,10000), ]

ggplot(gg_b_data2, aes(y = indus, x = test_town, color = xi, group = xi)) +
  
  geom_line() +
  
  labs(color = expression(xi)) +
  xlab("Perturbation town") +
  scale_color_viridis_d(option = "C", end = 0.95)





