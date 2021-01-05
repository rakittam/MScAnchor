# Define, run and store simulation data sets
#
# Date: 14.12.20     Author: Maic Rakitta
###############################################################################

library(glare)
library(simstudy)

# -----------------------------------------------------------------------------
### Function Definitions
# Function definition for anchor regression comparison ------------------------

onerep_rot <- function(rep, nobs = 300, data_table, data_pert_table, 
                       formula, A_formula,
                       xi_values = seq(-1, 2, by = 0.1), xi_big = 10000,
                       family) {
  
  # Generate data set with nobs observations
  dd <- genData(nobs, data_table)
  dd_pert <- genData(nobs, data_pert_table)
  
  # Hyperparameter values
  xi_values <- c(xi_values, xi_big)
  gamma_values <- 2 * xi_values + 1
  
  # Fit
  b_ar <- matrix(nrow = length(xi_values), ncol = 1)
  b_ar_se <- matrix(nrow = length(xi_values), ncol = 1)
  
  b_dev <- matrix(nrow = length(xi_values), ncol = 1)
  b_dev_se <- matrix(nrow = length(xi_values), ncol = 1)
  
  b_pea <- matrix(nrow = length(xi_values), ncol = 1)
  b_pea_se <- matrix(nrow = length(xi_values), ncol = 1)
  
  MSE_pert_ar <- numeric(length(xi_values))
  MSE_pert_dev <- numeric(length(xi_values))
  MSE_pert_pea <- numeric(length(xi_values))
  
  for (i in 1:length(xi_values)) {
    
    xi <- xi_values[i]
    gamma <- gamma_values[i]
    
    # Fit anchor regression
    fit_ar <- anchor_regression(formula = formula,
                                A_formula = A_formula,
                                data = dd,
                                gamma = gamma)
    
    b_ar[i] <- as.numeric(coef(fit_ar))
    b_ar_se[i] <- coef(summary(fit_ar))[, 2]
    resid_pert_ar <- dd_pert$Y - dd_pert$X * b_ar[i]
    MSE_pert_ar[i] <- mean((resid_pert_ar)^2)
    
    if (xi < 0) {
      b_dev[i] <- NA
      b_dev_se[i] <- NA
      b_dev[i] <- NA
      b_dev_se[i] <- NA
      MSE_pert_dev[i] <- NA
      MSE_pert_pea[i] <- NA
    } else {
      # Fit glare with deviance
      fit_dev <- glare(formula = formula,
                       A_formula = A_formula,
                       data = dd,
                       xi = xi,
                       family = family,
                       type = "deviance")
      
      b_dev[i] <- as.numeric(coef(fit_dev))
      b_dev_se[i] <- fit_dev$coef_se
      resid_pert_dev <- dd_pert$Y - predict(fit_dev,
                                            type = "response",
                                            newdata = dd_pert)
      MSE_pert_dev[i] <- mean((resid_pert_dev)^2)
      
      # Fit glare with pearson
      fit_pea <- glare(formula = formula,
                       A_formula = A_formula,
                       data = dd,
                       xi = xi,
                       family = family,
                       type = "pearson")
      
      b_pea[i] <- as.numeric(coef(fit_pea))
      b_pea_se[i] <- fit_pea$coef_se
      resid_pert_pea <- dd_pert$Y - predict(fit_pea,
                                            type = "response",
                                            newdata = dd_pert)
      MSE_pert_pea[i] <- mean((resid_pert_pea)^2)
    }
  }
  
  # Return
  data.frame(rep = rep,
             xi_values = xi_values,
             gamma_values = gamma_values,
             b_ar = b_ar,
             b_ar_se = b_ar_se,
             b_dev = b_dev,
             b_dev_se = b_dev_se,
             b_pea = b_pea,
             b_pea_se = b_pea_se,
             MSE_pert_ar = MSE_pert_ar,
             MSE_pert_dev = MSE_pert_dev,
             MSE_pert_pea = MSE_pert_pea)
}

# Define simulation function
simulate_rot <- function(nsim, nobs = 300, xi_values, xi_big = 10000,
                         data_table, data_pert_table, 
                         formula, A_formula, family) {

  xi_len <- length(xi_values) + 1
  
  states <- array(dim = c(nsim, 626))
  sim_data <- data.frame(matrix(nrow = nsim * xi_len, ncol = 12))
  colnames(sim_data) <- c("rep", "xi", "gamma",
                          "b_ar", "se_ar",
                          "b_dev", "se_dev",
                          "b_pea", "se_pea",
                          "MSE_pert_ar", "MSE_pert_dev", "MSE_pert_pea")
  
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  for (r in 1:nsim) {
    
    states[r, ] <- .Random.seed
    sim_data[(xi_len * (r - 1) + 1):(xi_len * r), ] <- 
      onerep_rot(rep = r, nobs = 300,
                 data_table = data_table, data_pert_table = data_pert_table, 
                 formula = formula, A_formula = A_formula,
                 xi_values = xi_values, xi_big = xi_big,
                 family = family)
    
    Sys.sleep(0.1)
    setTxtProgressBar(pb, r)
  }
  close(pb)
  
  list(states = states, sim_data = sim_data)
}

# Function definition for fixed v simulation for Rothenhaeusler ---------------

onerep_rot_fivi <- function(rep, nobs = 300, data_table, data_pert_table, 
                            formula, A_formula, xi_values, xi_big = 10000,
                            family, quant_value = 0.9) {
  
  # Generate data set with nobs observations
  dd <- genData(nobs, data_table)
  dd_pert <- genData(nobs, data_pert_table)
  
  xi_len <- length(xi_values)
  
  b_dev <- matrix(nrow = xi_len, ncol = 1)
  b_dev_se <- matrix(nrow = xi_len, ncol = 1)
  
  logLik_dev_pert <- numeric(xi_len)
  
  b_pea <- matrix(nrow = xi_len, ncol = 1)
  b_pea_se <- matrix(nrow = xi_len, ncol = 1)
  
  logLik_pea_pert <- numeric(xi_len)
  
  for (i in 1:xi_len) {
    
    xi <- xi_values[i]
    
    # Fit glare with deviance
    fit_temp <- glare(formula = formula,
                      A_formula = A_formula,
                      data = dd,
                      xi = xi,
                      family = family,
                      type = "deviance")
    
    b_dev[i] <- as.numeric(coef(fit_temp))
    b_dev_se[i] <- fit_temp$coef_se
    
    logLik_dev_pert_indiv <- logLik(fit_temp, newdata = dd_pert, indiv = TRUE)
    logLik_dev_pert[i] <-
      as.numeric(quantile(logLik_dev_pert_indiv, quant_value))
    
    # Fit glare with pearson
    fit_temp <- glare(formula = formula,
                      A_formula = A_formula,
                      data = dd,
                      xi = xi,
                      family = family,
                      type = "pearson")
    
    b_pea[i] <- as.numeric(coef(fit_temp))
    b_pea_se[i] <- fit_temp$coef_se
    
    logLik_pea_pert_indiv <- logLik(fit_temp, newdata = dd_pert, indiv = TRUE)
    logLik_pea_pert[i] <-
      as.numeric(quantile(logLik_pea_pert_indiv, quant_value))
  }
  
  # Fit GLM, i.e. glare with deviance
  fit_glm <- glare(formula = formula,
                   A_formula = A_formula,
                   data = dd,
                   xi = 0,
                   family = family,
                   type = "deviance")
  
  b_glm <- as.numeric(coef(fit_glm))
  b_glm_se <- fit_glm$coef_se
  
  logLik_glm_pert_indiv <- logLik(fit_glm, newdata = dd_pert, indiv = TRUE)
  logLik_glm_pert <- as.numeric(quantile(logLik_glm_pert_indiv, quant_value))
  
  # Fit glare with xi big for deviance
  fit_dev_big <- glare(formula = formula,
                       A_formula = A_formula,
                       data = dd,
                       xi = xi_big,
                       family = family,
                       type = "deviance")
  
  b_dev_big <- as.numeric(coef(fit_dev_big))
  b_dev_big_se <- fit_dev_big$coef_se
  
  logLik_dev_big_pert_indiv <- logLik(fit_dev_big, newdata = dd_pert, indiv = TRUE)
  logLik_dev_big_pert <- as.numeric(quantile(logLik_dev_big_pert_indiv, quant_value))
  
  # Fit glare for xi = 0 for pearson
  fit_pea_0 <- glare(formula = formula,
                     A_formula = A_formula,
                     data = dd,
                     xi = 0,
                     family = family,
                     type = "pearson")
  
  b_pea_0 <- as.numeric(coef(fit_pea_0))
  b_pea_0_se <- fit_pea_0$coef_se
  
  logLik_pea_0_pert_indiv <- logLik(fit_pea_0, newdata = dd_pert, indiv = TRUE)
  logLik_pea_0_pert <- as.numeric(quantile(logLik_pea_0_pert_indiv, quant_value))
  
  # Fit glare with xi big for pearson
  fit_pea_big <- glare(formula = formula,
                       A_formula = A_formula,
                       data = dd,
                       xi = xi_big,
                       family = family,
                       type = "pearson")
  
  b_pea_big <- as.numeric(coef(fit_pea_big))
  b_pea_big_se <- fit_pea_big$coef_se
  
  logLik_pea_big_pert_indiv <- logLik(fit_pea_big, newdata = dd_pert, indiv = TRUE)
  logLik_pea_big_pert <- as.numeric(quantile(logLik_pea_big_pert_indiv, quant_value))
  
  # Return
  data.frame(rep = rep,
             xi = c(0, xi_big, xi_values),
             b_dev = c(b_glm, b_dev_big, b_dev),
             b_dev_se = c(b_glm_se, b_dev_big_se, b_dev_se),
             logLik_dev_pert = c(logLik_glm_pert,
                                 logLik_dev_big_pert,
                                 logLik_dev_pert),
             b_pea = c(b_pea_0, b_pea_big, b_pea),
             b_pea_se = c(b_pea_0_se, b_pea_big_se, b_pea_se),
             logLik_pea_pert = c(logLik_pea_0_pert,
                                 logLik_pea_big_pert,
                                 logLik_pea_pert))
}

# Define simulation function
simulate_rot_fivi <- function(nsim, nobs = 300, xi_values, xi_big = 10000,
                              data_table, data_pert_table, 
                              formula, A_formula, family,
                              quant_value = 0.9) {
  
  xi_len <- length(xi_values)
  
  states <- array(dim = c(nsim, 626))
  
  sim_data <- data.frame(matrix(nrow = nsim * (xi_len + 2), ncol = 8))
  colnames(sim_data) <- c("rep", "xi",
                          "b_dev", "se_dev", "logLik_pert_dev",
                          "b_pea", "se_pea", "logLik_pert_pea")
  
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  for (r in 1:nsim) {
    states[r, ] <- .Random.seed
    sim_data[((xi_len + 2) * (r - 1) + + 1):((xi_len + 2) * r), ] <- 
      onerep_rot_fivi(rep = r,
                      data_table = data_table, data_pert_table = data_pert_table, 
                      formula = formula, A_formula = A_formula,
                      xi_values = xi_values, xi_big = xi_big,
                      family = family,
                      quant_value = quant_value)
    
    Sys.sleep(0.1)
    setTxtProgressBar(pb, r)
  }
  close(pb)
  
  list(states = states, sim_data = sim_data)
}

# Function definition for fixed xi --------------------------------------------

onerep_fixi <- function(nobs = 300, data_table, data_pert_table, 
                        formula, A_formula, xi_values, xi_big = 10000,
                        family, type, quant_value = 0.9, states) {
  
  
  # Generate data with nobs observations
  set.seed(6423) # use the same seeds for each v to use same errors
  dd <- genData(nobs, data_table)
  set.seed(6423)
  dd_pert <- genData(nobs, data_pert_table)
  
  # Fit GLM
  fit_glm <- glare(formula = formula,
                   A_formula = A_formula,
                   data = dd,
                   xi = 0,
                   family = family,
                   type = type)
  
  b_glm <- as.numeric(coef(fit_glm))
  b_glm_se <- fit_glm$coef_se
  
  logLik_glm_pert_indiv <- logLik(fit_glm, newdata = dd_pert, indiv = TRUE)
  logLik_glm_pert <- as.numeric(quantile(logLik_glm_pert_indiv, quant_value))
  

  
  # Fit glares
  xi_len <- length(xi_values)
  b_glare <- matrix(nrow = xi_len, ncol = 1)
  b_glare_se <- matrix(nrow = xi_len, ncol = 1)
  logLik_glare_pert <- matrix(nrow = xi_len, ncol = 1)
    
  for (j in 1:xi_len) {
    xi <- xi_values[j]
    
    fit_glare <- glare(formula = formula,
                       A_formula = A_formula,
                       data = dd,
                       xi = xi,
                       family = family,
                       type = type)
    
    b_glare[j] <- as.numeric(coef(fit_glare))
    b_glare_se[j] <- fit_glare$coef_se
    
    logLik_glare_pert_indiv <- logLik(fit_glare, newdata = dd_pert, indiv = TRUE)
    logLik_glare_pert[j] <- as.numeric(quantile(logLik_glare_pert_indiv, quant_value))
  }
  
  # Fit glare with xi big
  fit_big <- glare(formula = formula,
                   A_formula = A_formula,
                   data = dd,
                   xi = xi_big,
                   family = family,
                   type = type)
  
  b_big <- as.numeric(coef(fit_big))
  b_big_se <- fit_big$coef_se
  
  logLik_big_pert_indiv <- logLik(fit_big, newdata = dd_pert, indiv = TRUE)
  logLik_big_pert <- as.numeric(quantile(logLik_big_pert_indiv, quant_value))
  
  # Return
  data.frame(xi = c(0, xi_values, xi_big),
             b = c(b_glm, b_glare, b_big),
             se = c(b_glm_se, b_glare_se, b_big_se),
             logLik_pert = c(logLik_glm_pert, logLik_glare_pert, logLik_big_pert))
}


# Define simulation function
simulate_fixi <- function(nobs = 300, xi_values, xi_big = 10000, v_values, 
                          data_table, data_pert_table, 
                          formula, A_formula, family, type,
                          quant_value = 0.9) {
  
  v_len <- length(v_values)
  xi_len <- length(xi_values)

  sim_data <- data.frame(matrix(nrow = v_len * (xi_len + 2), ncol = 5))
  colnames(sim_data) <- c("v", "xi",
                          "b", "se", "logLik_pert")
  

  pb <- txtProgressBar(min = 0, max = v_len, style = 3)
  for (s in 1:v_len) {
    data_pert_table_temp <- updateDef(data_pert_table,
                                      changevar = "v",
                                      newformula = v_values[s])
    sim_data[((xi_len + 2) * (s - 1) + 1):(s * (xi_len + 2)), 1] <- v_values[s]
    
    sim_data[((xi_len + 2) * (s - 1) + 1):(s * (xi_len + 2)), 2:5] <-
      onerep_fixi(nobs = nobs,
                  data_table = data_table, data_pert_table = data_pert_table_temp, 
                  formula = formula, A_formula = A_formula,
                  xi_values = xi_values, xi_big = xi_big, family = family, type = type,
                  quant_value = quant_value, states = states[r, ])
    Sys.sleep(0.1)
    setTxtProgressBar(pb, s)
  }
  close(pb)
  
  
  
  list(sim_data = sim_data)
}

# Function definition for fixed v simulation ----------------------------------

onerep_fivi <- function(rep, nobs = 300, data_table, data_pert_table, 
                        formula, A_formula, xi_values, xi_big = 10000,
                        family, type, quant_value = 0.9) {
  
  # Generate data set with nobs observations
  dd <- genData(nobs, data_table)
  dd_pert <- genData(nobs, data_pert_table)
  
  xi_len <- length(xi_values)
  
  b_glare <- matrix(nrow = xi_len, ncol = 1)
  b_glare_se <- matrix(nrow = xi_len, ncol = 1)
  
  logLik_glare_pert <- numeric(xi_len)

  for (i in 1:xi_len) {
    
    xi <- xi_values[i]
    
    fit_temp <- glare(formula = formula,
                      A_formula = A_formula,
                      data = dd,
                      xi = xi,
                      family = family,
                      type = type)
    
    b_glare[i] <- as.numeric(coef(fit_temp))
    b_glare_se[i] <- fit_temp$coef_se
    
    logLik_glare_pert_indiv <- logLik(fit_temp, newdata = dd_pert, indiv = TRUE)
    logLik_glare_pert[i] <-
      as.numeric(quantile(logLik_glare_pert_indiv, quant_value))
  }
  
  # Fit GLM
  fit_glm <- glare(formula = formula,
                   A_formula = A_formula,
                   data = dd,
                   xi = 0,
                   family = family,
                   type = type)
  
  b_glm <- as.numeric(coef(fit_glm))
  b_glm_se <- fit_glm$coef_se
  
  logLik_glm_pert_indiv <- logLik(fit_glm, newdata = dd_pert, indiv = TRUE)
  logLik_glm_pert <- as.numeric(quantile(logLik_glm_pert_indiv, quant_value))
  
  
  # Fit glare with xi big
  fit_big <- glare(formula = formula,
                   A_formula = A_formula,
                   data = dd,
                   xi = xi_big,
                   family = family,
                   type = type)
  
  b_big <- as.numeric(coef(fit_big))
  b_big_se <- fit_big$coef_se
  
  logLik_big_pert_indiv <- logLik(fit_big, newdata = dd_pert, indiv = TRUE)
  logLik_big_pert <- as.numeric(quantile(logLik_big_pert_indiv, quant_value))
  
  # Return
  data.frame(rep = rep,
             xi = c(0, xi_big, xi_values),
             b = c(b_glm, b_big, b_glare),
             b_se = c(b_glm_se, b_big_se, b_glare_se),
             logLik_pert = c(logLik_glm_pert,
                             logLik_big_pert,
                             logLik_glare_pert))
}

# Define simulation function
simulate_fivi <- function(nsim, nobs = 300, xi_values, xi_big = 10000,
                          data_table, data_pert_table, 
                          formula, A_formula, family, type,
                          quant_value = 0.9) {
  
  xi_len <- length(xi_values)
  
  states <- array(dim = c(nsim, 626))
  
  sim_data <- data.frame(matrix(nrow = nsim * (xi_len + 2), ncol = 5))
  colnames(sim_data) <- c("rep", "xi",
                          "b", "se", "logLik_pert")
  
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  for (r in 1:nsim) {
    states[r, ] <- .Random.seed
    sim_data[((xi_len + 2) * (r - 1) + + 1):((xi_len + 2) * r), ] <- 
      onerep_fivi(rep = r,
                  data_table = data_table, data_pert_table = data_pert_table, 
                  formula = formula, A_formula = A_formula,
                  xi_values = xi_values, xi_big = xi_big,
                  family = family, type = type,
                  quant_value = quant_value)
    
    Sys.sleep(0.1)
    setTxtProgressBar(pb, r)
  }
  close(pb)
  
  list(states = states, sim_data = sim_data)
}

# Function definition for fixed v simulation with quantiles -------------------

onerep_fivi_list <- function(rep, nobs = 300, data_table, data_pert_table, 
                             formula, A_formula, xi_values, xi_big = 10000,
                             family, type) {
  
  # Generate data set with nobs observations
  dd <- genData(nobs, data_table)
  dd_pert <- genData(nobs, data_pert_table)
  
  xi_len <- length(xi_values)
  
  b_glare <- matrix(nrow = xi_len, ncol = 1)
  b_glare_se <- matrix(nrow = xi_len, ncol = 1)
  
  logLik_glare_pert <- matrix(nrow = xi_len, ncol = nobs)

  for (i in 1:xi_len) {  
    
    xi <- xi_values[i]
    
    fit_temp <- glare(formula = formula,
                      A_formula = A_formula,
                      data = dd,
                      xi = xi,
                      family = family,
                      type = type)
    
    b_glare[i] <- as.numeric(coef(fit_temp))
    b_glare_se[i] <- fit_temp$coef_se
    
    logLik_glare_pert[i, ] <- logLik(fit_temp, newdata = dd_pert, indiv = TRUE)

  }
  
  # Fit GLM
  fit_glm <- glare(formula = formula,
                   A_formula = A_formula,
                   data = dd,
                   xi = 0,
                   family = family,
                   type = type)
  
  b_glm <- as.numeric(coef(fit_glm))
  b_glm_se <- fit_glm$coef_se

  logLik_glm_pert <- logLik(fit_glm, newdata = dd_pert, indiv = TRUE)
  
  # Fit glare with xi big
  fit_big <- glare(formula = formula,
                   A_formula = A_formula,
                   data = dd,
                   xi = xi_big,
                   family = family,
                   type = type)
  
  b_big <- as.numeric(coef(fit_big))
  b_big_se <- fit_big$coef_se
  
  logLik_big_pert <- logLik(fit_big, newdata = dd_pert, indiv = TRUE)
  
  # Return
  coeffs <- data.frame(rep = rep,
                             xi = c(0, xi_big, xi_values),
                             b = c(b_glm, b_big, b_glare),
                             b_se = c(b_glm_se, b_big_se, b_glare_se))
  
  logLik_pert <- rbind(logLik_glm_pert, logLik_glare_pert, logLik_big_pert)
  logLiks <- cbind(c(0, xi_values, xi_big), logLik_pert)
  
  list(coeffs, logLiks)
}

# Define simulation function
simulate_fivi_list <- function(nsim, nobs = 300, xi_values, xi_big = 10000,
                               data_table, data_pert_table, 
                               formula, A_formula, family, type) {
  
  xi_len <- length(xi_values)
  
  states <- array(dim = c(nsim, 626))
  
  sim_data <- list()
  
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  for (r in 1:nsim) {
    states[r, ] <- .Random.seed
    sim_data[[r]] <- 
      onerep_fivi_list(rep = r,
                       data_table = data_table, data_pert_table = data_pert_table, 
                       formula = formula, A_formula = A_formula,
                       xi_values = xi_values, xi_big = xi_big,
                       family = family, type = type)
    
    Sys.sleep(0.1)
    setTxtProgressBar(pb, r)
  }
  close(pb)
  
  list(states = states, sim_data = sim_data)
}

# -----------------------------------------------------------------------------
### Generate data and store
# sim1: Rothenhaeusler Comparison ---------------------------------------------

# Define variables for unperturbed and perturbed data set
def_rot <- defData(varname = "randA", dist = "binary",
                   formula = 0.5)
def_rot <- defData(def_rot, varname = "H", dist = "normal", variance = 1,
                   formula = 0)
def_rot <- defData(def_rot, varname = "A", 
                   formula = "2 * randA - 1")
def_rot <- defData(def_rot, varname = "X", dist = "normal", variance = 1,
                   formula = "A + H")
def_rot <- defData(def_rot, varname = "Y", dist = "normal", variance = 1,
                   formula = "X + 2 * H")

def_rot_pert <- defData(varname = "randA", dist = "binary",
                        formula = 0.5)
def_rot_pert <- defData(def_rot_pert, varname = "H", dist = "normal",
                        variance = 1,
                        formula = 0)
def_rot_pert <- defData(def_rot_pert, varname = "A", 
                        formula = "2 * randA - 1")
def_rot_pert <- defData(def_rot_pert, varname = "v", 
                        formula = 1.8)
def_rot_pert <- defData(def_rot_pert, varname = "X", dist = "normal", 
                        variance = 1,
                        formula = "v + H")
def_rot_pert <- defData(def_rot_pert, varname = "Y", dist = "normal", 
                        variance = 1,
                        formula = "X + 2 * H")

dd <- genData(300, def_rot)
dd_pert <- genData(300, def_rot_pert)
hist(dd$Y)

# Initialize for Rothenhaeusler comparison
set.seed(1813)

# Simulate
sim_data_rot <- simulate_rot(nsim = 100, nobs = 1000,
                             xi_values = seq(-0.5, 50, by = 0.1), xi_big = 10000,
                             data_table = def_rot, data_pert_table = def_rot_pert, 
                             formula = Y ~ X - 1, A_formula = ~ A - 1,
                             family = gaussian)

data_rot_states <- sim_data_rot$states
data_rot <- sim_data_rot$sim_data

# Store states and data
# path_name <- "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/sim1/sim1_comp_states"
# write.table(data_rot_states, file = paste(paste(path_name, Sys.Date(), sep = "_"), ".csv", sep =""))
# path_name <- "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/sim1/sim1_comp_data"
# write.table(data_rot, file = paste(paste(path_name, Sys.Date(), sep = "_"), ".csv", sep =""))

# sim2: Anchor on X normal Rothenhaeusler (IV setting) ------------------------

def_rot <- defData(varname = "randA", dist = "binary",
                   formula = 0.5)
def_rot <- defData(def_rot, varname = "H", dist = "normal", variance = 1,
                   formula = 0)
def_rot <- defData(def_rot, varname = "A", 
                   formula = "2 * randA - 1")
def_rot <- defData(def_rot, varname = "X", dist = "normal", variance = 1,
                   formula = "A + H")
def_rot <- defData(def_rot, varname = "Y", dist = "normal", variance = 1,
                   formula = "X + 2 * H")

def_rot_pert <- defData(varname = "randA", dist = "binary",
                        formula = 0.5)
def_rot_pert <- defData(def_rot_pert, varname = "H", dist = "normal",
                        variance = 1,
                        formula = 0)
def_rot_pert <- defData(def_rot_pert, varname = "A", 
                        formula = "2 * randA - 1")
def_rot_pert <- defData(def_rot_pert, varname = "v", 
                        formula = 1.8)
def_rot_pert <- defData(def_rot_pert, varname = "X", dist = "normal", 
                        variance = 1,
                        formula = "v + H")
def_rot_pert <- defData(def_rot_pert, varname = "Y", dist = "normal", 
                        variance = 1,
                        formula = "X + 2 * H")

dd <- genData(300, def_rot)
dd_pert <- genData(300, def_rot_pert)
hist(dd$Y)

# Initialize for fixed v
set.seed(8653)
nsim <- 100

xi_values <- seq(0, 10, by = 0.1)
xi_values <- xi_values[-1]
xi_big <- 10000

data_table <- def_rot
data_pert_table <- def_rot_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- gaussian

# Simulate
sim_data_rot_X_fivi <- simulate_rot_fivi(nsim = nsim, nobs = 300,
                                         xi_values = xi_values, xi_big = xi_big,
                                         data_table = data_table,
                                         data_pert_table = data_pert_table, 
                                         formula = formula, A_formula = A_formula,
                                         family = family,
                                         quant_value = 0.9) 

data_rot_X_states_fivi <- sim_data_rot_X_fivi$states
data_rot_X_fivi <- sim_data_rot_X_fivi$sim_data

# Store states and data
path_name <- "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/sim2/sim2_fivi_states"
write.table(data_rot_X_states_fivi, file = paste(paste(path_name, Sys.Date(), sep = "_"), ".csv", sep =""))
path_name <- "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/sim2/sim2_fivi_data"
write.table(data_rot_X_fivi, file = paste(paste(path_name, Sys.Date(), sep = "_"), ".csv", sep =""))

# Initialize for fixed xi
set.seed(7863)

xi_values <- c(0.2, 0.5, 0.8)
xi_big <- 10000
v_values <- seq(0, 10, by = 0.1)

data_table <- def_rot
data_pert_table <- def_rot_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- gaussian
type <- "pearson"

# Simulate
sim_data_rot_X_fixi <- simulate_fixi(nobs = 300,
                                     xi_values = xi_values, xi_big = xi_big,
                                     v_values = v_values, 
                                     data_table = data_table,
                                     data_pert_table = data_pert_table, 
                                     formula = formula, A_formula = A_formula,
                                     family = family, type = type,
                                     quant_value = 0.9)

data_rot_X_fixi <- sim_data_rot_X_fixi$sim_data

# Store states and data
path_name <- "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/sim2/sim2_fixi_data"
write.table(data_rot_X_fixi, file = paste(paste(path_name, Sys.Date(), sep = "_"), ".csv", sep =""))

# sim3: Anchor on X poisson (IV setting) --------------------------------------

# Define variables for unperturbed and perturbed data set
def_poi_X <- defData(varname = "A", dist = "normal", 
                     formula = 0, variance = 0.25)
def_poi_X <- defData(def_poi_X, varname = "H", dist = "normal",
                     formula = 0, variance = 0.25)
def_poi_X <- defData(def_poi_X, varname = "X", dist = "normal", 
                     formula = "1.2 * H + A", variance = 0.25)
def_poi_X <- defData(def_poi_X, varname = "Y", dist = "poisson", link = "log", 
                     formula = "0.4 * X + 0.5 * H", variance = 1)

def_poi_X_pert <- defData(varname = "A", dist = "normal", 
                          formula = 0, variance = 0.25)
def_poi_X_pert <- defData(def_poi_X_pert, varname = "v", 
                          formula = 4) # set perturbation strength
def_poi_X_pert <- defData(def_poi_X_pert, varname = "H", dist = "normal",
                          formula = 0, variance = 0.25)
def_poi_X_pert <- defData(def_poi_X_pert, varname = "X", dist = "normal", 
                          formula = "1.3 * H + v", variance = 0.25)
def_poi_X_pert <- defData(def_poi_X_pert, varname = "Y", dist = "poisson", link = "log",
                          formula = "0.4 * X + 0.5 * H", variance = 1)

dd <- genData(300, def_poi_X)
dd_pert <- genData(300, def_poi_X_pert)
hist(dd$Y)
hist(dd_pert$Y)

# Initialize for fixed v
set.seed(14981)
nsim <- 100

xi_values <- seq(0, 10, by = 0.1)
xi_values <- xi_values[-1]
xi_big <- 10000

data_table <- def_poi_X
data_pert_table <- def_poi_X_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- poisson
type <- "pearson"

# nobs <- 300
# rep <- 1
# data_table <- def_poi_X
# data_pert_table <- def_poi_X_pert

# Simulate
sim_data_poi_X_fivi <- simulate_fivi(nsim = nsim, nobs = 300,
                                     xi_values = xi_values, xi_big = xi_big,
                                     data_table = data_table,
                                     data_pert_table = data_pert_table, 
                                     formula = formula, A_formula = A_formula,
                                     family = family, type = type,
                                     quant_value = 0.9) 

data_poi_X_states_fivi <- sim_data_poi_X_fivi$states
data_poi_X_fivi <- sim_data_poi_X_fivi$sim_data

# Store states and data
path_name <- "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/sim3/sim3_fivi_states"
write.table(data_poi_X_states_fivi, file = paste(paste(path_name, Sys.Date(), sep = "_"), ".csv", sep =""))
path_name <- "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/sim3/sim3_fivi_data"
write.table(data_poi_X_fivi, file = paste(paste(path_name, Sys.Date(), sep = "_"), ".csv", sep =""))

# Initialize for fixed xi
set.seed(33524)

xi <- 1
xi_big <- 10000
v_values <- seq(-5, 5, by = 1)

data_table <- def_poi_X
data_pert_table <- def_poi_X_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- poisson
type <- "pearson"

# Simulate
sim_data_poi_X_fixi <- simulate_fixi(nobs = 300,
                                     xi_values = xi, xi_big = xi_big,
                                     v_values = v_values, 
                                     data_table = data_table,
                                     data_pert_table = data_pert_table, 
                                     formula = formula, A_formula = A_formula,
                                     family = family, type = type,
                                     quant_value = 0.9)

data_poi_X_fixi <- sim_data_poi_X_fixi$sim_data

# Store states and data
path_name <- "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/sim3/sim3_fixi_data"
write.table(data_poi_X_fixi, file = paste(paste(path_name, Sys.Date(), sep = "_"), ".csv", sep =""))


# sim3.5: Anchor on X poisson (IV setting) plotting quantiles -----------------

# Define variables for unperturbed and perturbed data set
def_poi_X <- defData(varname = "A", dist = "normal", 
                     formula = 0, variance = 0.25)
def_poi_X <- defData(def_poi_X, varname = "H", dist = "normal",
                     formula = 0, variance = 0.25)
def_poi_X <- defData(def_poi_X, varname = "X", dist = "normal", 
                     formula = "1.2 * H + A", variance = 0.25)
def_poi_X <- defData(def_poi_X, varname = "Y", dist = "poisson", link = "log", 
                     formula = "0.4 * X + 0.5 * H", variance = 1)

def_poi_X_pert <- defData(varname = "A", dist = "normal", 
                          formula = 0, variance = 0.25)
def_poi_X_pert <- defData(def_poi_X_pert, varname = "v", 
                          formula = 4) # set perturbation strength
def_poi_X_pert <- defData(def_poi_X_pert, varname = "H", dist = "normal",
                          formula = 0, variance = 0.25)
def_poi_X_pert <- defData(def_poi_X_pert, varname = "X", dist = "normal", 
                          formula = "1.3 * H + v", variance = 0.25)
def_poi_X_pert <- defData(def_poi_X_pert, varname = "Y", dist = "poisson", link = "log",
                          formula = "0.4 * X + 0.5 * H", variance = 1)

dd <- genData(300, def_poi_X)
dd_pert <- genData(300, def_poi_X_pert)
hist(dd$Y)
hist(dd_pert$Y)

# Initialize for fixed v
set.seed(14981)
nsim <- 10

xi_values <- seq(0, 10, by = 2)
xi_values <- xi_values[-1]
xi_big <- 10000

data_table <- def_poi_X
data_pert_table <- def_poi_X_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- poisson
type <- "pearson"

# nobs <- 300
# rep <- 1
# data_table <- def_poi_X
# data_pert_table <- def_poi_X_pert

# Simulate
sim_data_poi_X_fivi_list <- simulate_fivi_list(nsim = nsim, nobs = 300,
                                               xi_values = xi_values, xi_big = xi_big,
                                               data_table = data_table,
                                               data_pert_table = data_pert_table, 
                                               formula = formula, A_formula = A_formula,
                                               family = family, type = type) 

# data_poi_X_states_fivi_list <- sim_data_poi_X_fivi_list$states
# data_poi_X_fivi_list <- sim_data_poi_X_fivi_list$sim_data

path_name <- "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/sim3/"
dir.create(paste(path_name, Sys.Date(), sep = ""))
save(sim_data_poi_X_fivi_list, file = paste(paste(path_name, Sys.Date(), sep = ""), "/sim_data_poi_X_fivi_list.Rdata", sep =""))

# Initialize for fixed xi
set.seed(33524)

xi <- 1
xi_big <- 10000
v_values <- seq(-5, 5, by = 1)

data_table <- def_poi_X
data_pert_table <- def_poi_X_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- poisson
type <- "pearson"

# Simulate
sim_data_poi_X_fixi <- simulate_fixi(nobs = 300,
                                     xi_values = xi, xi_big = xi_big,
                                     v_values = v_values, 
                                     data_table = data_table,
                                     data_pert_table = data_pert_table, 
                                     formula = formula, A_formula = A_formula,
                                     family = family, type = type,
                                     quant_value = 0.9)

data_poi_X_fixi <- sim_data_poi_X_fixi$sim_data

# Store states and data
path_name <- "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/sim3/sim3_fixi_data"
write.table(data_poi_X_fixi, file = paste(paste(path_name, Sys.Date(), sep = "_"), ".csv", sep =""))


# sim4: Anchor on X, H and Y poisson ------------------------------------------

# Define variables for unperturbed and perturbed data set
def_poi_XHY <- defData(varname = "A", dist = "normal", variance = 0.25,
                       formula = 0)
def_poi_XHY <- defData(def_poi_XHY, varname = "H", dist = "normal",
                       variance = 0.25,
                       formula = "0.6 + A",)
def_poi_XHY <- defData(def_poi_XHY, varname = "X", dist = "normal",
                       variance = 0.25,
                       formula = "H + A")
def_poi_XHY <- defData(def_poi_XHY, varname = "Y",
                       dist = "poisson", link = "log", 
                       formula = "0.8 * X - H - A")
# HOW TO INTERVERNE?
def_poi_XHY_pert <- defData(varname = "A", dist = "normal", 
                            variance = 0.25,
                            formula = 0)
def_poi_XHY_pert <- defData(def_poi_XHY_pert, varname = "H", dist = "normal",
                            variance = 0.25,
                            formula = "0.6 + 0.2 * A") # set perturbation
def_poi_XHY_pert <- defData(def_poi_XHY_pert, varname = "X", dist = "normal", 
                            variance = 0.25,
                            formula = "H - 1 * A") # set perturbation 
def_poi_XHY_pert <- defData(def_poi_XHY_pert, varname = "Y",
                            dist = "poisson", link = "log", 
                            formula = "0.8 * X - H - 2 * A") # set perturbation 

dd <- genData(300, def_poi_XHY)
dd_pert <- genData(300, def_poi_XHY_pert)
hist(dd$Y)

# Initialize
set.seed(32416)
nsim <- 100

xi_values <- seq(0, 10, by = 2)
xi_values <- xi_values[-1]
xi_big <- 10000

data_table <- def_poi_XHY
data_pert_table <- def_poi_XHY_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- poisson
type <- "pearson"

# Simulate
sim_data_poi_XHY_fivi <- simulate_fivi(nsim = nsim, nobs = 300,
                                       xi_values = xi_values, xi_big = xi_big,
                                       data_table = data_table,
                                       data_pert_table = data_pert_table, 
                                       formula = formula, A_formula = A_formula,
                                       family = family, type = type,
                                       quant_value = 0.9)

data_poi_XHY_states_fivi <- sim_data_poi_XHY_fivi$states
data_poi_XHY_fivi <- sim_data_poi_XHY_fivi$sim_data

# Store states and data
path_name <- "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/sim4/sim4_fivi_states"
write.table(data_poi_XHY_states_fivi, file = paste(paste(path_name, Sys.Date(), sep = "_"), ".csv", sep =""))
path_name <- "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/sim4/sim4_fivi_data"
write.table(data_poi_XHY_fivi, file = paste(paste(path_name, Sys.Date(), sep = "_"), ".csv", sep =""))

# varying quantiles
sim_data_poi_XHY_fivi_list <- simulate_fivi_list(nsim = nsim, nobs = 300,
                                                 xi_values = xi_values, xi_big = xi_big,
                                                 data_table = data_table,
                                                 data_pert_table = data_pert_table, 
                                                 formula = formula, A_formula = A_formula,
                                                 family = family, type = type)

path_name <- "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/sim4/"
dir.create(paste(path_name, Sys.Date(), sep = ""))
save(sim_data_poi_XHY_fivi_list, file = paste(paste(path_name, Sys.Date(), sep = ""), "/sim_data_poi_XHY_fivi_list.Rdata", sep =""))
