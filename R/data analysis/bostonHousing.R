# GLARE for Boston Housing Data
#
# Date: 12.02.21     Author: Maic Rakitta
###############################################################################

library(glare)
library(mlbench)
library(ggplot2)
library(tidyr)
library(ggpubr) # for theme and ggarrange
library(dplyr) # for labels in boxplot
library(tibble) # for labels in boxplot

theme_set(theme_bw())

# Load data
data("BostonHousing2")

# standartize covariates and normalize response
std_covariates <- scale(BostonHousing2[, c("crim", "zn", "indus", "nox",
                                           "rm", "age", "lstat")])
nor_response <- (BostonHousing2[, "cmedv"] - min(BostonHousing2[, "cmedv"])) /
  (max(BostonHousing2[, "cmedv"]) - min(BostonHousing2[, "cmedv"]))

df_temp <- data.frame(cmedv = nor_response, town = BostonHousing2[, "town"])

bostonPrices <- data.frame(std_covariates[, c("crim", "zn", "indus", "nox",
                                              "rm", "age", "lstat")])
bostonPrices <- cbind(df_temp, bostonPrices)

# Create shortcuts for towns
names.arg <- levels(bostonPrices$town)
town_lookup_table <- abbreviate(names.arg, minlength = 2, use.classes = TRUE,
                                dot = FALSE, strict = FALSE,
                                method = c("left.kept", "both.sides"), named = TRUE)
levels(bostonPrices$town) <- town_lookup_table

# Explore data ----------------------------------------------------------------
plot(bostonPrices)
summary(bostonPrices)

table(bostonPrices$town)

library(Hmisc)
hist.data.frame(BostonHousing2[, c("cmedv", "crim", "zn", "indus", "nox",
                                   "rm", "age", "lstat")])

hist(bostonPrices$cmedv) # right skewed and >0 suggests log-transform
hist(log(bostonPrices$cmedv))

# Correlation of continuous variables
cor(bostonPrices[, c("cmedv", "crim", "zn", "indus", "nox",
                     "rm", "age", "lstat")], bostonPrices$cmedv)
# -> lstat high negative corr, rm high positive corr

# Fit GLARE -------------------------------------------------------------------

# Initialize
perturbation_towns <- names(table(bostonPrices$town))
pert_town_len <- length(perturbation_towns)

xi_values <- c(0, 1, 5, 10)
# xi_values <- c(0, 1, 5, 10, 50, 100, 10000)
xi_len <- length(xi_values)

fit_list <- list()
convergence_message <- matrix(nrow = pert_town_len,
                              ncol = xi_len)

for (t in 1:pert_town_len) {
  
  print(paste(t, "of", pert_town_len))
  
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
  
  logLik <- numeric(xi_len)
  
  pb_xi <- txtProgressBar(min = 0, max = xi_len, style = 3)
  for (i in 1:xi_len) {
    
    xi <- xi_values[i]
    
    fit_temp <- glare(formula = cmedv ~ crim + zn + indus + nox + rm +
                        age + lstat,
                      A_formula = ~ town, data = train_set, xi = xi,
                      family = gaussian(link = "log"), type = "deviance")
    
    convergence_message[t, i] <- fit_temp$optim$convergence
    
    # Coefficients
    b[i, ] <- as.numeric(coef(fit_temp))
    b_se[i, ] <- fit_temp$coef_se
    
    # Predictions, Residuals and MSE
    predictions[i, ] <- predict(fit_temp, type = "response", newdata = test_set)
    resid_pert[i, ] <- test_set$cmedv - predictions[i, ]
    RMSE[i] <- sqrt(mean((resid_pert[i, ])^2))
    
    # Likelihood
    logLik_indiv <- logLik(fit_temp, newdata = test_set, indiv = TRUE)
    logLik[i] <- mean(logLik_indiv)
    
    Sys.sleep(0.1)
    setTxtProgressBar(pb_xi, i)
  }
  close(pb_xi)
  
  fit_list[[t]] <- list(b, b_se, predictions, RMSE, logLik = logLik)
}

# path_name <- "./data sets/data_examples/bostonHousing/"
# dir.create(paste(path_name, Sys.Date(), sep = ""))
# save(fit_list, file = paste(paste(path_name, Sys.Date(), sep = ""), "/fit_list.Rdata", sep =""))

# Plots -----------------------------------------------------------------------

# Prepare data
load("./data sets/data_examples/bostonHousing/2021-02-11/fit_list.Rdata")

gg_data <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(gg_data) <- c("test_town", "xi", "logLik", "RMSE")

b_matrix <- matrix(0, nrow = 0, ncol = 10)
colnames(b_matrix) <- c("xi", "Town", "Intercept", "crim", "zn", "indus",
                        "nox", "rm", "age","lstat")
b_matrix <- data.frame(b_matrix)
b_matrix <- transform(b_matrix, xi = as.character(xi))
b_matrix <- transform(b_matrix, Town = as.character(Town))

for (t in 1:pert_town_len) {
  
  test_town <- perturbation_towns[t]
  
  for (i in 1:xi_len) {
    
    xi <- xi_values[i]
    
    logLik <- fit_list[[t]][["logLik"]][i]
    RMSE <- fit_list[[t]][[4]][i]
    b <- fit_list[[t]][[1]][i, ]
    
    data_temp <- data.frame(test_town, xi, logLik, RMSE)
    gg_data <- rbind(gg_data, data_temp)
    
    b_frame <- data.frame(t(b))
    b_frame <- cbind(xi, test_town, b_frame)
    b_matrix <- rbind(b_matrix, b_frame)
  }
}

# Boxplot without labels
gg_data$xi <- factor(gg_data$xi, levels = xi_values)

gg_boxplot <- ggplot(gg_data, aes(y = logLik, x = xi, color = xi, group = xi)) +
  geom_boxplot()
gg_boxplot

ggsave(filename = "boxplot.pdf", plot = gg_boxplot, height = 4, width = 6)

# Boxplot with labels
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

dat <- gg_data %>% mutate(outlier = test_town) %>% group_by(xi) %>% mutate(is_outlier=ifelse(is_outlier(logLik), test_town, as.numeric(NA)))
dat$outlier[which(is.na(dat$is_outlier))] <- as.numeric(NA)

gg_boxplot_label <- ggplot(dat, aes(y = logLik, x = xi, color = xi, group = xi)) +
  geom_boxplot() +
  geom_text(aes(label = outlier), na.rm = TRUE, nudge_x = 0.3, nudge_y = 0)
gg_boxplot_label

ggsave(filename = "boxplot_labels.pdf", plot = gg_boxplot_label, height = 4, width = 6)

# Quantile Plots
alpha_values <- seq(0.8, 1, by = 0.01)

logLik_quant <- data.frame(matrix(nrow = 0, ncol = 3)) 
                           
for (i in 1:xi_len) {
  xi <- xi_values[i]
  
  logLik_xi <- gg_data$logLik[gg_data$xi == xi]
  logLik_quant_temp <- as.numeric(quantile(logLik_xi, probs = alpha_values))
  logLik_quant <- rbind(logLik_quant, cbind(xi, alpha_values, logLik_quant_temp))
}
colnames(logLik_quant) <- c("xi", "alpha", "logLik")
logLik_quant$xi <- factor(logLik_quant$xi) 

gg_quantiles <- ggplot(logLik_quant, aes(y = logLik, x = xi, color = alpha, group = alpha)) +
  geom_line() +
  labs(color = expression(alpha)) +
  ylab("logLik") +
  xlab(expression(xi)) 
  #scale_color_viridis_d(option = "C", end = 0.95)
gg_quantiles 

ggsave(filename = "quantiles.pdf", plot = gg_quantiles, height = 4, width = 6)












# old -------------------------------------------------------------------------
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

# old -----------------------------------------------------------------------------
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





