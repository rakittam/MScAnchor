# Using GLARE for bioChemists data
#
# Date: 04.01.21     Author: Maic Rakitta
###############################################################################

library(glare)
library(pscl)

data("bioChemists")
summary(bioChemists)

# Create observed and unobserved environment
train_set <- bioChemists[bioChemists$phd >= 2, ]
test_set <- bioChemists[bioChemists$phd < 2, ]

# Create observed and unobserved environment
train_set <- bioChemists[bioChemists$phd <= 4, ]
test_set <- bioChemists[bioChemists$phd > 4, ]

# Initialize and run GLARE
xi_values <- seq(0, 100, by = 1)
xi_len <- length(xi_values)

b <- matrix(nrow = xi_len, ncol = 5)
b_se <- matrix(nrow = xi_len, ncol = 5)
logLik_pert <- numeric(xi_len)

pb <- txtProgressBar(min = 0, max = xi_len, style = 3)
for (i in 1:xi_len) {
  
  xi <- xi_values[i]
  
  fit_temp <- glare(formula = art ~ fem + mar + kid5 + ment,
                    A_formula = ~ phd, data = train_set, xi = xi,
                    family = poisson, type = "pearson")
  
  b[i, ] <- as.numeric(coef(fit_temp))
  b_se[i, ] <- fit_temp$coef_se
  
  logLik_pert_indiv <- logLik(fit_temp, newdata = train_set, indiv = TRUE)
  logLik_pert[i] <- as.numeric(quantile(logLik_pert_indiv, 0.9))
  
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
}
close(pb)


# Plots
gg_data <- data.frame(logLik_pert = logLik_pert, xi = xi_values)

ggplot(gg_data, aes(y = logLik_pert, x = xi)) +
  geom_line() +
  xlab(expression(xi))
