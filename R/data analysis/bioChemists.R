# GLARE for bioChemists data
#
# Date: 04.01.21     Author: Maic Rakitta
###############################################################################

library(glare)
library(pscl)
library(ggplot2)
library(tidyverse) # for data manipulation 

data("bioChemists")
summary(bioChemists)
plot(bioChemists)

hist(bioChemists$art)
hist(bioChemists$kid5)
hist(bioChemists$phd)
hist(bioChemists$ment)

table(bioChemists$fem)
table(bioChemists$mar)

# Data preparations -----------------------------------------------------------

# Create observed and unobserved environment
# train_set <- bioChemists[bioChemists$phd >= 2, ]
# test_set <- bioChemists[bioChemists$phd < 2, ]

# Create observed and unobserved environment
train_set <- bioChemists[bioChemists$phd <= 4, ]
test_set <- bioChemists[bioChemists$phd > 4, ]

# Create observed and unobserved environment
# index <- which(bioChemists$phd < 4 & bioChemists$phd > 1)
# train_set <- bioChemists[index, ]
# test_set <- bioChemists[-index, ]

# Initialize and run GLARE
xi_values <- seq(0, 20, by = 0.1)
xi_len <- length(xi_values)

# Fit GLARE -------------------------------------------------------------------
b <- matrix(nrow = xi_len, ncol = 5)
colnames(b) <- c("Intercept", "fem", "mar", "kid5", "ment")
b_se <- matrix(nrow = xi_len, ncol = 5)
colnames(b_se) <- colnames(b)

logLik_pert_indiv_matrix <- matrix(nrow = xi_len, ncol = nrow(test_set))
logLik_pert <- numeric(xi_len)

pb <- txtProgressBar(min = 0, max = xi_len, style = 3)
for (i in 1:xi_len) {
  
  xi <- xi_values[i]
  
  fit_temp <- glare(formula = art ~ fem + mar + kid5 + ment,
                    A_formula = ~ phd, data = train_set, xi = xi,
                    family = poisson, type = "pearson")
  
  b[i, ] <- as.numeric(coef(fit_temp))
  b_se[i, ] <- fit_temp$coef_se
  
  logLik_pert_indiv <- logLik(fit_temp, newdata = test_set, indiv = TRUE)
  logLik_pert[i] <- as.numeric(quantile(logLik_pert_indiv, 0.9))
  
  logLik_pert_indiv_matrix[i, ] <- as.numeric(logLik_pert_indiv)
  
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
}
close(pb)

# xi big
xi_big <- 10000

fit_big <- glare(formula = art ~ fem + mar + kid5 + ment,
                 A_formula = ~ phd, data = train_set, xi = xi_big,
                 family = poisson, type = "pearson")

b_big <- as.numeric(coef(fit_big))
b_se_big <- fit_big$coef_se

logLik_pert_indiv_big <- logLik(fit_big, newdata = test_set, indiv = TRUE)
logLik_pert_big <- as.numeric(quantile(logLik_pert_indiv_big, 0.9))

logLik_pert_big
b_big

# Plots -----------------------------------------------------------------------

# Parameter plot
gg_data_b <- as.data.frame(cbind(xi_values, b))

df_b <- gg_data_b %>%
  select(xi_values, Intercept, fem, mar, kid5, ment) %>%
  gather(key = "Variable", value = "value", -xi_values)
head(df_b)

ggplot(df_b, aes(x = xi_values, y = value)) + 
  geom_line(aes(color = Variable)) +
  ylab("Parameters") +
  xlab(expression(xi)) +
  geom_vline(xintercept = 0, linetype = "dashed")

# Single quantile 0.5 - 0.9
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 5
cols = gg_color_hue(n)
alpha_values <- seq(0, 1, by = 0.1)
quantile_matrix <- apply(logLik_pert_indiv_matrix, 1, quantile, probs = alpha_values)

alp <- 0.5
gg_data_log2 <- data.frame(logLik_pert = quantile_matrix[alp * 10 + 1, ], xi = xi_values)
ggplot(gg_data_log2, aes(y = logLik_pert, x = xi)) +
  geom_line(color = cols[alp * 10 - 4]) +
  ylab(paste(alp, "-quantile of log-likelihood")) +
  xlab(expression(xi)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(color = NULL)

# All quantiles plot
df_temp <- as.data.frame(quantile_matrix)
colnames(df_temp) <- xi_values
df_temp <- df_temp[6:10, ]
df_temp <- cbind(rownames(df_temp), df_temp)
rownames(df_temp) <- NULL
colnames(df_temp)[1] <- 'alpha'

library(reshape2)
df_melted = melt(df_temp, id.vars = 'alpha')
df_melted$variable <- as.numeric(as.character(df_melted$variable))

ggplot(df_melted, aes(y = value, x = variable, color = alpha, group = alpha)) +
  geom_line() +
  ylab("Quantiles of log-likelihood") +
  xlab(expression(xi)) +
  labs(color = expression(alpha))





aa <- alpha_values

alpha_values <- seq(0,1,by=0.01)


# Individual comparison -------------------------------------------------------

# Fit GLM
fit_glm_orig <- glm(formula = art ~ fem + mar + kid5 + ment,
               data = train_set, family = poisson)
summary(fit_glm_orig)

xi <- 0
fit_glm <- glare(formula = art ~ fem + mar + kid5 + ment,
                 A_formula = ~ phd, data = train_set, xi = xi,
                 family = poisson, type = "pearson")
summary(fit_glm)

logLik_pert_indiv_glm <- as.numeric(logLik(fit_glm, newdata = test_set,
                                             indiv = TRUE))
quantiles_glm <-apply(as.matrix(logLik_pert_indiv_glm), 2, quantile,
                        probs = alpha_values)
quantiles_glm


# Fit GLARE
xi <- 100
fit_glare <- glare(formula = art ~ fem + mar + kid5 + ment,
             A_formula = ~ phd, data = train_set, xi = xi,
             family = poisson, type = "pearson")
summary(fit_glare)

logLik_pert_indiv_glare <- as.numeric(logLik(fit_glare, newdata = test_set,
                                             indiv = TRUE))
quantiles_glare <-apply(as.matrix(logLik_pert_indiv_glare), 2, quantile,
                        probs = alpha_values)
quantiles_glare

# Plot
gg_data_quantiles <- data.frame(log_quantiles = c(as.numeric(quantiles_glare),
                                                  as.numeric(quantiles_glm)),
                                alpha = c(alpha_values, alpha_values),
                                group = c(rep("GLARE", length(alpha_values)),
                                          rep("GLM", length(alpha_values))))
ggplot(gg_data_quantiles, aes(y = log_quantiles, x = alpha,
                              color = group, group = group)) +
  geom_line() +
  ylab("Quantiles of log-likelihood") +
  xlab(expression(alpha)) +
  labs(color = "Fit")

gg_data_quantiles2 <- gg_data_quantiles[gg_data_quantiles$alpha >= 0.5, ]
ggplot(gg_data_quantiles2, aes(y = log_quantiles, x = alpha,
                              color = group, group = group)) +
  geom_line() +
  ylab("Quantiles of log-likelihood") +
  xlab(expression(alpha)) +
  labs(color = "Fit")
