# Define, run and store simulation data sets
#
# Date: 26.02.21     Author: Maic Rakitta
###############################################################################
# Description -----------------------------------------------------------------
#
# The file simulation_functions.R provides the simulations for each simulation
#  example in the thesis. First load the function definitions and packages and 
#  then head to the corresponding example of interest.
#
# The file plot_functions.R stores the plot functions used in this script. For
#  visualizing results only, simply execute the plot_functions.R file and the
#  call the corresponding plot function. If interests in specific graphics with
#  descriptions, load the data manually in the plot function and use the commen-
#  ted code.
#
# The file run_analysis.R is the main file that calls stored or loaded data and
#  executes plot functions or makes analysis.
#
# -----------------------------------------------------------------------------
### Libraries
library(glare)
library(simstudy)

# -----------------------------------------------------------------------------
### Function Definitions

# Function definition for Anchor Regression Comparison ------------------------

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
  
  max_y_hat_diff <- matrix(nrow = length(xi_values), ncol = 1)
  
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
      
      
      # estimate comparance:
      y_hat_ar <- b_ar[i] * dd_pert$X
      y_hat_dev <- b_dev[i] * dd_pert$X
      
      max_y_hat_diff[i] <- max(abs(y_hat_ar - y_hat_dev))
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
             MSE_pert_pea = MSE_pert_pea,
             max_y_hat_diff = max_y_hat_diff)
}

# Define simulation function
simulate_rot <- function(nsim, nobs = 300, xi_values, xi_big = 10000,
                         data_table, data_pert_table, 
                         formula, A_formula, family) {
  
  xi_len <- length(xi_values) + 1
  
  states <- array(dim = c(nsim, 626))
  sim_data <- data.frame(matrix(nrow = nsim * xi_len, ncol = 13))
  colnames(sim_data) <- c("rep", "xi", "gamma",
                          "b_ar", "se_ar",
                          "b_dev", "se_dev",
                          "b_pea", "se_pea",
                          "MSE_pert_ar", "MSE_pert_dev", "MSE_pert_pea",
                          "max_y_hat_diff")
  
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

# Function definition for fixed v Simulation for Rothenhaeusler ---------------

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
  
  logLik_dev_big_pert_indiv <- logLik(fit_dev_big, newdata = dd_pert,
                                      indiv = TRUE)
  logLik_dev_big_pert <- as.numeric(quantile(logLik_dev_big_pert_indiv,
                                             quant_value))
  
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
  
  logLik_pea_big_pert_indiv <- logLik(fit_pea_big, newdata = dd_pert,
                                      indiv = TRUE)
  logLik_pea_big_pert <- as.numeric(quantile(logLik_pea_big_pert_indiv,
                                             quant_value))
  
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



# Function definition for fixed v ---------------------------------------------

onerep_fivi <- function(rep, nobs = 300, data_table, data_pert_table, 
                        formula, A_formula, xi_values,
                        family, type, quant_value = 0.9) {
  
  # Generate data set with nobs observations
  dd <- genData(nobs, data_table)
  dd_pert <- genData(nobs, data_pert_table)
  
  xi_len <- length(xi_values)
  
  b <- matrix(nrow = xi_len, ncol = 1)
  b_se <- matrix(nrow = xi_len, ncol = 1)
  
  logLik_pert <- numeric(xi_len)
  
  for (i in 1:xi_len) {
    
    xi <- xi_values[i]
    
    fit_temp <- glare(formula = formula,
                      A_formula = A_formula,
                      data = dd,
                      xi = xi,
                      family = family,
                      type = type)
    
    b[i] <- as.numeric(coef(fit_temp))
    b_se[i] <- fit_temp$coef_se
    
    logLik_pert_indiv <- -logLik(fit_temp, newdata = dd_pert, indiv = TRUE)
    logLik_pert[i] <-
      as.numeric(quantile(logLik_pert_indiv, quant_value))
  }
  
  # Return
  data.frame(rep = rep,
             xi = xi_values,
             b = b,
             b_se = b_se,
             logLik_pert = logLik_pert)
}

# Define simulation function
simulate_fivi <- function(nsim, nobs = 300, xi_values,
                          data_table, data_pert_table, 
                          formula, A_formula, family, type,
                          quant_value = 0.9) {
  
  xi_len <- length(xi_values)
  
  states <- array(dim = c(nsim, 626))
  
  sim_data <- data.frame(matrix(nrow = 0, ncol = 5))
  colnames(sim_data) <- c("rep", "xi",
                          "b", "se", "logLik_pert")
  
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  for (r in 1:nsim) {
    states[r, ] <- .Random.seed
    sim_data_temp <- 
      onerep_fivi(rep = r,
                  data_table = data_table, data_pert_table = data_pert_table, 
                  formula = formula, A_formula = A_formula,
                  xi_values = xi_values,
                  family = family, type = type,
                  quant_value = quant_value)
    
    sim_data <- rbind(sim_data, sim_data_temp)
    
    Sys.sleep(0.1)
    setTxtProgressBar(pb, r)
  }
  close(pb)
  
  list(states = states, sim_data = sim_data)
}

# Function definition for fixed xi --------------------------------------------

onerep_fixi <- function(rep, nobs = 300, data_table, data_pert_table, 
                        formula, A_formula, xi_values, v_values,
                        family, type, quant_value = 0.9, rep_seed) {
  
  v_len <- length(v_values)
  
  data_frame <- data.frame()
  indv_list <- list()
  
  for (s in 1:v_len) {
    data_pert_table_temp <- updateDef(data_pert_table,
                                      changevar = "v",
                                      newformula = v_values[s])
    
    # Generate data with nobs observations
    set.seed(rep_seed) # use the same seeds for each v to use same errors
    dd <- genData(nobs, data_table)
    set.seed(rep_seed + 200)
    dd_pert <- genData(nobs, data_pert_table_temp)
    
    # Fit glares
    xi_len <- length(xi_values)
    b <- matrix(nrow = xi_len, ncol = 1)
    b_se <- matrix(nrow = xi_len, ncol = 1)
    logLik_pert_indiv <- matrix(nrow = xi_len, ncol = nobs)
    logLik_pert <- matrix(nrow = xi_len, ncol = 1)
    
    devres <- matrix(nrow = xi_len, ncol = nobs)
    deviance <- matrix(nrow = xi_len, ncol = 1)
    
    peares <- matrix(nrow = xi_len, ncol = nobs)
    pearson <- matrix(nrow = xi_len, ncol = 1)
    
    for (j in 1:xi_len) {
      xi <- xi_values[j]
      
      fit <- glare(formula = formula,
                   A_formula = A_formula,
                   data = dd,
                   xi = xi,
                   family = family,
                   type = type)
      
      b[j] <- as.numeric(coef(fit))
      b_se[j] <- fit$coef_se
      
      logLik_pert_indiv[j, ] <- -logLik(fit, newdata = dd_pert, indiv = TRUE)
      logLik_pert[j] <- as.numeric(quantile(logLik_pert_indiv[j, ], quant_value))
      
      if(family == "poisson") {
        
        link_inv <- poisson()$linkinv
        vari <- poisson()$variance
        mu <- link_inv(dd_pert$X * b[j])
        
        r <- mu
        p <- which(dd_pert$Y > 0)
        r[p] <- (dd_pert$Y * log(dd_pert$Y / mu) - (dd_pert$Y - mu))[p]
        
        devres[j, ] <- sign(dd_pert$Y - mu) * sqrt(2 * r)
        
      } else if(family == "binomial") {
        
        link_inv <- binomial()$linkinv
        vari <- binomial()$variance
        mu <- link_inv(dd_pert$X * b[j])
        
        r <- dd_pert$Y * log(dd_pert$Y / (mu)) + (1 - dd_pert$Y) *
          log((1 - dd_pert$Y) / (1 - mu))
        p <- which(dd_pert$Y == 0)
        r[p] <- ((1 - dd_pert$Y) * log((1 - dd_pert$Y) / (1 - mu)))[p]
        q <- which(dd_pert$Y == 1)
        r[q] <- (dd_pert$Y * log(dd_pert$Y / (mu)))[q]
        
        devres[j, ] <- sign(dd_pert$Y / 1 - mu) * sqrt(2 * r)
        
      } else if(family == "gaussian"){
        
        link_inv <- gaussian()$linkinv
        vari <- gaussian()$variance
        mu <- link_inv(dd_pert$X * b[j])
        
        n <- length(dd_pert$Y)
        s_hat <- sqrt(1/n * sum((dd_pert$Y-mu)^2))
        
        devres[j, ] <- 1 / s_hat * (dd_pert$Y - mu)
        
      } else {
        stop("No valid family! Use either 'poisson', 'binomial' or 'gaussian'.")
      }
      
      # Deviance
      deviance[j] <- sqrt(mean(devres[j, ]^2))
      
      # Pearson
      V <- vari(mu)
      peares[j, ] <- (dd_pert$Y - mu)/sqrt(V)
      pearson[j] <-  sqrt(mean(peares[j, ]^2))
    }
    
    # Return
    data_frame_temp <- data.frame(rep = rep, v = v_values[s], xi = xi_values,
                                  b = b,
                                  se = b_se,
                                  logLik_pert = logLik_pert,
                                  deviance = deviance,
                                  pearson = pearson)
    
    data_frame <- rbind(data_frame, data_frame_temp)
    
    indv_list[[s]] <- list(logLik_pert_indiv, devres, peares,
                           xi_values, v = v_values[s])
  }
  
  list(data = data_frame,
       indv = indv_list)
}


# Define simulation function
simulate_fixi <- function(nsim = 100, nobs = 300, xi_values, v_values, 
                          data_table, data_pert_table, 
                          formula, A_formula, family, type,
                          quant_value = 0.9, start_seed) {
  
  data_list <- list()
  
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  for (r in 1:nsim) {
    
    data_temp <-
      onerep_fixi(rep = r, nobs = nobs,
                  data_table = data_table,
                  data_pert_table = data_pert_table, 
                  formula = formula, A_formula = A_formula,
                  xi_values = xi_values, v_values = v_values,
                  family = family, type = type,
                  quant_value = quant_value, rep_seed = r + start_seed)
    
    data_list[[r]] <- list(rep = r,
                           data = data_temp)
    
    Sys.sleep(0.1)
    setTxtProgressBar(pb, r)
  }
  close(pb)
  
  data_list
}

# Function definition for IV causal parameter estimation ----------------------

onerep_IV <- function(rep, nobs_values = 300, data_table, 
                      formula, A_formula, xi_values,
                      family, type) {
  
  data <- as.data.frame(matrix(ncol = 5, nrow = 0))
  colnames(data) <- c("rep", "nobs", "xi", "b", "b_se")
  
  for (k in 1:length(nobs_values)) {
    
    nobs <- nobs_values[k]
    
    # Generate data set with nobs observations
    dd <- genData(nobs, data_table)
    
    xi_len <- length(xi_values)
    
    b <- matrix(nrow = xi_len, ncol = 1)
    b_se <- matrix(nrow = xi_len, ncol = 1)
    
    for (i in 1:xi_len) {
      
      xi <- xi_values[i]
      
      fit_temp <- glare(formula = formula,
                        A_formula = A_formula,
                        data = dd,
                        xi = xi,
                        family = family,
                        type = type)
      
      b[i] <- as.numeric(coef(fit_temp))
      b_se[i] <- fit_temp$coef_se
      
    }
    
    # Return
    data_temp <- data.frame(rep = rep, nobs = nobs,
                            xi = xi_values,
                            b = b,
                            b_se = b_se)
    
    data <- rbind(data, data_temp)
  }
  
  data
}

# Define simulation function
simulate_IV <- function(nsim, nobs_values = 300, xi_values,
                        data_table, formula, A_formula, family, type) {
  
  xi_len <- length(xi_values)
  
  states <- array(dim = c(nsim, 626))
  
  b_data <- as.data.frame(matrix(ncol = 5, nrow = 0))
  colnames(b_data) <- c("rep", "nobs", "xi", "b", "b_se")
  
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  for (r in 1:nsim) {
    states[r, ] <- .Random.seed
    
    data_temp <- onerep_IV(rep = r, nobs_values = nobs_values,
                           data_table = data_table, 
                           formula = formula, A_formula = A_formula,
                           xi_values = xi_values,
                           family = family, type = type)
    
    b_data <- rbind(b_data, data_temp)
    
    Sys.sleep(0.1)
    setTxtProgressBar(pb, r)
  }
  close(pb)
  
  list(states = states, sim_data = b_data)
}


# Function definition for fixed v simulation with quantiles -------------------

onerep_fivi_quant <- function(rep, nobs = 300, data_table, data_pert_table, 
                              formula, A_formula, xi_values, xi_big = 10000,
                              family, type) {
  
  # Generate data set with nobs observations
  dd <- genData(nobs, data_table)
  dd_pert <- genData(nobs, data_pert_table)
  
  # number of parameters
  n_b <- ncol(model.matrix(formula, data = dd))
  
  xi_len <- length(xi_values)
  
  b_glare <- matrix(nrow = xi_len, ncol = n_b)
  b_glare_se <- matrix(nrow = xi_len, ncol = n_b)
  
  logLik_glare_pert <- matrix(nrow = xi_len, ncol = nobs)
  
  for (i in 1:xi_len) {  
    
    xi <- xi_values[i]
    
    fit_temp <- glare(formula = formula,
                      A_formula = A_formula,
                      data = dd,
                      xi = xi,
                      family = family,
                      type = type)
    
    b_glare[i, ] <- as.numeric(coef(fit_temp))
    b_glare_se[i, ] <- fit_temp$coef_se
    
    logLik_glare_pert[i, ] <- -logLik(fit_temp, newdata = dd_pert, indiv = TRUE)
    
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
  
  logLik_glm_pert <- -logLik(fit_glm, newdata = dd_pert, indiv = TRUE)
  
  # Fit glare with xi big
  fit_big <- glare(formula = formula,
                   A_formula = A_formula,
                   data = dd,
                   xi = xi_big,
                   family = family,
                   type = type)
  
  b_big <- as.numeric(coef(fit_big))
  b_big_se <- fit_big$coef_se
  
  logLik_big_pert <- -logLik(fit_big, newdata = dd_pert, indiv = TRUE)
  
  # Return
  coeffs <- list(rep = rep,
                 xi = c(0, xi_big, xi_values),
                 b = rbind(b_glm, b_big, b_glare),
                 b_se = rbind(b_glm_se, b_big_se, b_glare_se))
  
  logLik_pert <- rbind(logLik_glm_pert, logLik_glare_pert, logLik_big_pert)
  logLiks <- cbind(c(0, xi_values, xi_big), logLik_pert)
  
  list(coeffs, logLiks)
}

# Define simulation function
simulate_fivi_quant <- function(nsim, nobs = 300, xi_values, xi_big = 10000,
                                data_table, data_pert_table, 
                                formula, A_formula, family, type) {
  
  xi_len <- length(xi_values)
  
  states <- array(dim = c(nsim, 626))
  
  sim_data <- list()
  
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  for (r in 1:nsim) {
    states[r, ] <- .Random.seed
    sim_data[[r]] <- 
      onerep_fivi_quant(rep = r,
                        data_table = data_table,
                        data_pert_table = data_pert_table, 
                        formula = formula, A_formula = A_formula,
                        xi_values = xi_values, xi_big = xi_big,
                        family = family, type = type)
    
    Sys.sleep(0.1)
    setTxtProgressBar(pb, r)
  }
  close(pb)
  
  list(states = states, sim_data = sim_data)
}

# Function definition for fixed xi for glm with A -----------------------------

onerep_fixi_glm <- function(rep, nobs = 300, data_table, data_pert_table, 
                            formula, A_formula, xi_values, v_values,
                            family, type, quant_value = 0.9, rep_seed) {
  
  v_len <- length(v_values)
  
  data_frame <- data.frame()
  indv_list <- list()
  
  for (s in 1:v_len) {
    data_pert_table_temp <- updateDef(data_pert_table,
                                      changevar = "v",
                                      newformula = v_values[s])
    
    # Generate data with nobs observations
    set.seed(rep_seed) # use the same seeds for each v to use same errors
    dd <- genData(nobs, data_table)
    set.seed(rep_seed + 200)
    dd_pert <- genData(nobs, data_pert_table_temp)
    
    # Fit glares
    xi_len <- length(xi_values)
    b <- matrix(nrow = xi_len, ncol = 2)
    b_se <- matrix(nrow = xi_len, ncol = 2)
    logLik_pert_indiv <- matrix(nrow = xi_len, ncol = nobs)
    logLik_pert <- matrix(nrow = xi_len, ncol = 1)
    
    devres <- matrix(nrow = xi_len, ncol = nobs)
    deviance <- matrix(nrow = xi_len, ncol = 1)
    
    peares <- matrix(nrow = xi_len, ncol = nobs)
    pearson <- matrix(nrow = xi_len, ncol = 1)
    
    for (j in 1:xi_len) {
      xi <- xi_values[j]
      
      fit <- glare(formula = formula,
                   A_formula = A_formula,
                   data = dd,
                   xi = xi,
                   family = family,
                   type = type)
      
      b[j, ] <- as.numeric(coef(fit))
      b_se[j, ] <- fit$coef_se
      
      logLik_pert_indiv[j, ] <- logLik(fit, newdata = dd_pert, indiv = TRUE)
      logLik_pert[j] <- as.numeric(quantile(logLik_pert_indiv[j, ], quant_value))
      
      if(family == "poisson") {
        
        link_inv <- poisson()$linkinv
        vari <- poisson()$variance
        mu <- link_inv(dd_pert$X * b[j, ])
        
        r <- mu
        p <- which(dd_pert$Y > 0)
        r[p] <- (dd_pert$Y * log(dd_pert$Y / mu) - (dd_pert$Y - mu))[p]
        
        devres[j, ] <- sign(dd_pert$Y - mu) * sqrt(2 * r)
        
      } else if(family == "binomial") {
        
        link_inv <- binomial()$linkinv
        vari <- binomial()$variance
        mu <- link_inv(dd_pert$X * b[j, ])
        
        r <- dd_pert$Y * log(dd_pert$Y / (mu)) + (1 - dd_pert$Y) *
          log((1 - dd_pert$Y) / (1 - mu))
        p <- which(dd_pert$Y == 0)
        r[p] <- ((1 - dd_pert$Y) * log((1 - dd_pert$Y) / (1 - mu)))[p]
        q <- which(dd_pert$Y == 1)
        r[q] <- (dd_pert$Y * log(dd_pert$Y / (mu)))[q]
        
        devres[j, ] <- sign(dd_pert$Y / 1 - mu) * sqrt(2 * r)
        
      } else if(family == "gaussian"){
        
        link_inv <- gaussian()$linkinv
        vari <- gaussian()$variance
        mu <- link_inv(dd_pert$X * b[j, ])
        
        n <- length(dd_pert$Y)
        s_hat <- sqrt(1/n * sum((dd_pert$Y-mu)^2))
        
        devres[j, ] <- 1 / s_hat * (dd_pert$Y - mu)
        
      } else {
        stop("No valid family! Use either 'poisson', 'binomial' or 'gaussian'.")
      }
      
      # Deviance
      deviance[j] <- sqrt(mean(devres[j, ]^2))
      
      # Pearson
      V <- vari(mu)
      peares[j, ] <- (dd_pert$Y - mu)/sqrt(V)
      pearson[j] <-  sqrt(mean(peares[j, ]^2))
    }
    
    # Return
    data_frame_temp <- data.frame(rep = rep, v = v_values[s], xi = xi_values,
                                  logLik_pert = logLik_pert,
                                  deviance = deviance,
                                  pearson = pearson)
    
    data_frame <- rbind(data_frame, data_frame_temp)
    
    indv_list[[s]] <- list(logLik_pert_indiv, devres, peares,
                           xi_values, v = v_values[s], b = b, b_se =b_se)
  }
  
  list(data = data_frame,
       indv = indv_list)
}


# Define simulation function
simulate_fixi_glm <- function(nsim = 100, nobs = 300, xi_values, v_values, 
                              data_table, data_pert_table, 
                              formula, A_formula, family, type,
                              quant_value = 0.9, start_seed) {
  
  data_list <- list()
  
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  for (r in 1:nsim) {
    
    data_temp <-
      onerep_fixi_glm(rep = r, nobs = nobs,
                      data_table = data_table,
                      data_pert_table = data_pert_table, 
                      formula = formula, A_formula = A_formula,
                      xi_values = xi_values, v_values = v_values,
                      family = family, type = type,
                      quant_value = quant_value, rep_seed = r + start_seed)
    
    data_list[[r]] <- list(rep = r,
                           data = data_temp)
    
    Sys.sleep(0.1)
    setTxtProgressBar(pb, r)
  }
  close(pb)
  
  data_list
}

# Function definition for fixed v simulation with quantiles for glm with A ----

onerep_fivi_quant_glm <- function(rep, nobs = 300, data_table, data_pert_table, 
                                  formula, A_formula, xi_values, xi_big = 10000,
                                  family, type) {
  
  # Generate data set with nobs observations
  dd <- genData(nobs, data_table)
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
  
  logLik_glm_pert <- -logLik(fit_glm, newdata = dd_pert, indiv = TRUE)
  
  # Return
  coeffs <- data.frame(b = b_glm,
                       b_se = b_glm_se)
  
  logLik_pert <- logLik_glm_pert
  logLiks <- cbind(-1, logLik_pert)
  
  list(coeffs, logLiks)
}

# Define simulation function
simulate_fivi_quant_glm <- function(nsim, nobs = 300, xi_values, xi_big = 10000,
                                    data_table, data_pert_table, 
                                    formula, A_formula, family, type) {
  
  states <- array(dim = c(nsim, 626))
  
  sim_data <- list()
  
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  for (r in 1:nsim) {
    states[r, ] <- .Random.seed
    sim_data[[r]] <- 
      onerep_fivi_quant_glm(rep = r,
                            data_table = data_table,
                            data_pert_table = data_pert_table, 
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

# ex1: Rothenhaeusler Example -------------------------------------------------

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

# Parameters of linear predictor
b <- 1
h <- 2

# Distribution histogram
dd <- genData(300, def_rot)
dd_pert <- genData(300, def_rot_pert)
par(mfrow = c(2,2))
hist(b*dd$X + h * dd$H, breaks = 100)
hist(dd$Y, breaks = 100)
hist(b*dd_pert$X + h * dd_pert$H, breaks = 100)
hist(dd_pert$Y, breaks = 100)

# Initialize
data_table <- def_rot
data_pert_table <- def_rot_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- "gaussian"
type <- "deviance"

# sim1: Rothenhaeusler Comparison ---

set.seed(1813)

nsim <- 100

xi_values <- seq(-0.5, 50, by = 0.1)
xi_big <-  10000

# Simulate
sim_data_rot <- simulate_rot(nsim = nsim, nobs = 1000,
                             xi_values = xi_values, xi_big = xi_big,
                             data_table = data_table,
                             data_pert_table = data_pert_table, 
                             formula = formula, A_formula = A_formula,
                             family = family)

# yhat diff 
temp_data <- sim_data_rot$sim_data
# # Save simulated data list
# path_name <- "./data sets/simulation_study/ex1/"
# dir.create(paste(path_name, Sys.Date(), sep = ""))
# save(temp_data, file = paste(paste(path_name, Sys.Date(), sep = ""),
#                              "/temp_data_yhat_diff.Rdata", sep =""))

xi_vec <- c(xi_values, xi_big)
max_y_diff <- matrix(nrow = length(xi_vec))
b_vec <- matrix(nrow = length(xi_vec))
for (i in 1:length(xi_vec)) {
  xi <- xi_vec[i]
  
  b_vec[i] <- mean(temp_data$b_ar[temp_data$xi==xi])
  max_y_diff[i] <- mean(temp_data$max_y_hat_diff[temp_data$xi==xi])
  
}

plot(b_vec, max_y_diff)


# # Save simulated data list
# path_name <- "./data sets/simulation_study/ex1/"
# dir.create(paste(path_name, Sys.Date(), sep = ""))
# save(sim_data_rot, file = paste(paste(path_name, Sys.Date(), sep = ""),
#                                 "/sim1.Rdata", sep =""))

# sim2: Rothenhaeusler with fixed v ---

set.seed(8653)
nsim <- 100

xi_values <- seq(0, 10, by = 0.1)
xi_values <- xi_values[-1]
xi_big <- 10000

# Simulate
sim_data_rot_X_fivi <- simulate_rot_fivi(nsim = nsim, nobs = 1000,
                                         xi_values = xi_values, xi_big = xi_big,
                                         data_table = data_table,
                                         data_pert_table = data_pert_table, 
                                         formula = formula,
                                         A_formula = A_formula,
                                         family = family,
                                         quant_value = 0.9) 

# # Save simulated data list
# path_name <- "./data sets/simulation_study/ex1/"
# dir.create(paste(path_name, Sys.Date(), sep = ""))
# save(sim_data_rot_X_fivi, file = paste(paste(path_name, Sys.Date(), sep = ""),
                                       # "/sim2.Rdata", sep =""))

# sim3: Rothenhaeusler with fixed xi ---

start_seed <- 7863

nsim <- 100

xi_values <- c(0, 1, 2, 10000)
v_values <- seq(0, 10, by = 0.1)

# Simulate
sim_data_rot_X_fixi <- simulate_fixi(nsim = nsim, nobs = 1000,
                                     xi_values = xi_values,
                                     v_values = v_values, 
                                     data_table = data_table,
                                     data_pert_table = data_pert_table, 
                                     formula = formula, A_formula = A_formula,
                                     family = family, type = type,
                                     quant_value = 0.9, start_seed = start_seed)

# # Save simulated data list
# path_name <- "./data sets/simulation_study/ex1/"
# dir.create(paste(path_name, Sys.Date(), sep = ""))
# save(sim_data_rot_X_fixi, file = paste(paste(path_name, Sys.Date(), sep = ""),
#                                        "/sim3.Rdata", sep =""))

# ex2: Poisson IV Example -----------------------------------------------------

# Define variables for unperturbed and perturbed data set
def_poi_X <- defData(varname = "A", dist = "normal", 
                     formula = 0, variance = 0.25)
def_poi_X <- defData(def_poi_X, varname = "H", dist = "normal",
                     formula = 0, variance = 0.25)
def_poi_X <- defData(def_poi_X, varname = "X", dist = "normal", 
                     formula = "2 * H + A", variance = 0.25)
def_poi_X <- defData(def_poi_X, varname = "Y", dist = "poisson", link = "log", 
                     formula = "0.4 * X + 1 * H", variance = 1)

def_poi_X_pert <- defData(varname = "A", dist = "normal", 
                          formula = 0, variance = 0.25)
def_poi_X_pert <- defData(def_poi_X_pert, varname = "v", 
                          formula = 3) # set perturbation strength
def_poi_X_pert <- defData(def_poi_X_pert, varname = "H", dist = "normal",
                          formula = 0, variance = 0.25)
def_poi_X_pert <- defData(def_poi_X_pert, varname = "X", dist = "normal", 
                          formula = "2 * H + v", variance = 0.25)
def_poi_X_pert <- defData(def_poi_X_pert, varname = "Y", dist = "poisson",
                          link = "log",
                          formula = "0.4 * X + 1 * H", variance = 1)

# Parameters of linear predictor
b <- 0.4
h <- 1

# Distribution histogram
dd <- genData(300, def_poi_X)
dd_pert <- genData(300, def_poi_X_pert)
par(mfrow = c(2,2))
hist(b*dd$X + h * dd$H, breaks = 100)
hist(dd$Y, breaks = 100)
hist(b*dd_pert$X + h * dd_pert$H, breaks = 100)
hist(dd_pert$Y, breaks = 100)

# Initialize
data_table <- def_poi_X
data_pert_table <- def_poi_X_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- "poisson"
type <- "deviance"

# sim1: Poisson IV for fixed v ---

set.seed(35732)

nsim <- 100

xi_values <- seq(0, 10, by = 0.1)
#xi_values <- c(xi_values, 10000)
xi_values2 <- seq(20, 250, by = 10)
xi_values <- c(xi_values, xi_values2, 10000)

# Simulate
sim_data_poi_X_fivi <- simulate_fivi(nsim = nsim, nobs = 1000,
                                     xi_values = xi_values,
                                     data_table = data_table,
                                     data_pert_table = data_pert_table, 
                                     formula = formula, A_formula = A_formula,
                                     family = family, type = type,
                                     quant_value = 0.9) 

# # Save data
# path_name <- "./data sets/simulation_study/ex2/"
# dir.create(paste(path_name, Sys.Date(), sep = ""))
# save(sim_data_poi_X_fivi, file = paste(paste(path_name, Sys.Date(), sep = ""), 
#                                        "/sim1.Rdata", sep =""))

# sim2: Poisson IV for fixed xi ---

start_seed <- 96574

nsim <- 100

xi_values <- c(0, 1, 3, 5, 8, 10, 20, 50, 10000)
v_values <- seq(-10, 10, by = 0.1)

# Simulate
sim_data_poi_X_fixi <- simulate_fixi(nsim = nsim, nobs = 1000,
                                     xi_values = xi_values,
                                     v_values = v_values, 
                                     data_table = data_table,
                                     data_pert_table = data_pert_table, 
                                     formula = formula, A_formula = A_formula,
                                     family = family, type = type,
                                     quant_value = 0.9, start_seed = start_seed)

# # Save data
# path_name <- "./data sets/simulation_study/ex2/"
# dir.create(paste(path_name, Sys.Date(), sep = ""))
# save(sim_data_poi_X_fixi, file = paste(paste(path_name, Sys.Date(), sep = ""), 
#                                        "/sim2.Rdata", sep =""))

# sim3: IV investigation for poisson ---

set.seed(88435)

nsim <- 100

nobs_values <- c(seq(10, 100, by = 10),
                 seq(100, 5000, by = 100)[-1])

# nobs_values <- c(300, 5000)

# xi_values <- seq(0, 100, by = 10)
xi_values <- c(0, 1, 3, 10, 10000)

# Simulate
IV_b_poi <- simulate_IV(nsim = nsim, nobs_values = nobs_values,
                        xi_values = xi_values,
                        data_table = data_table,
                        formula = formula, A_formula = A_formula,
                        family = family, type = type)

# # Save data
# path_name <- "./data sets/simulation_study/ex2/"
# dir.create(paste(path_name, Sys.Date(), sep = ""))
# save(IV_b_poi, file = paste(paste(path_name, Sys.Date(), sep = ""), 
#                             "/sim3.Rdata", sep =""))

# sim4: Poisson IV several quantiles ---

set.seed(77781)

nsim <- 100

xi_values <- seq(0, 10, by = 0.1)
xi_values <- xi_values[-1]
xi_big <- 10000

# Simulate
sim_data_poi_X_quant <- simulate_fivi_quant(nsim = nsim, nobs = 1000,
                                            xi_values = xi_values,
                                            xi_big = xi_big,
                                            data_table = data_table,
                                            data_pert_table = data_pert_table, 
                                            formula = formula,
                                            A_formula = A_formula,
                                            family = family, type = type) 

# # Save data
# path_name <- "./data sets/simulation_study/ex2/"
# dir.create(paste(path_name, Sys.Date(), sep = ""))
# save(sim_data_poi_X_quant, file = paste(paste(path_name, Sys.Date(), sep = ""),
#                                         "/sim4.Rdata", sep =""))




# ex3: Binomial IV Example ----------------------------------------------------

# Define variables for unperturbed and perturbed data set
def_bin_X <- defData(varname = "A", dist = "normal", 
                     formula = 0, variance = 0.25)
def_bin_X <- defData(def_bin_X, varname = "H", dist = "normal",
                     formula = 0, variance = 0.25)
def_bin_X <- defData(def_bin_X, varname = "X", dist = "normal", 
                     formula = "2 * H + A", variance = 0.25)
def_bin_X <- defData(def_bin_X, varname = "Y", dist = "binomial", link = "logit", 
                     formula = "0.4 * X + 1 * H", variance = 1)

def_bin_X_pert <- defData(varname = "A", dist = "normal", 
                          formula = 0, variance = 0.25)
def_bin_X_pert <- defData(def_bin_X_pert, varname = "v", 
                          formula = 8) # set perturbation strength
def_bin_X_pert <- defData(def_bin_X_pert, varname = "H", dist = "normal",
                          formula = 0, variance = 0.25)
def_bin_X_pert <- defData(def_bin_X_pert, varname = "X", dist = "normal", 
                          formula = "2 * H + v", variance = 0.25)
def_bin_X_pert <- defData(def_bin_X_pert, varname = "Y", dist = "binomial",
                          link = "logit",
                          formula = "0.4 * X + 1 * H", variance = 1)

# Parameters of linear predictor
b <- 0.4
h <- 1

v <- 8



fit_temp <- glare(formula = formula,
                  A_formula = A_formula,
                  data = dd,
                  xi = 0,
                  family = family,
                  type = type)
fit_temp2 <- glare(formula = formula,
                   A_formula = A_formula,
                   data = dd,
                   xi = 5,
                   family = family,
                   type = type)

predict(fit_temp2, type = "response")

# Distribution histogram
dd <- genData(300, def_bin_X)
dd_pert <- genData(300, def_bin_X_pert)

eta <- b*dd$X + h * dd$H
eta_pert <- b*dd_pert$X + h * dd_pert$H
par(mfrow = c(2,3))

hist(eta, breaks = 100)
hist(exp(eta)/(1+exp(eta)), breaks = 100)
hist(dd$Y, breaks = 100)

hist(eta_pert, breaks = 100)
hist(exp(eta_pert)/(1+exp(eta_pert)), breaks = 100)
hist(dd_pert$Y, breaks = 100)

# Initialize
data_table <- def_bin_X
data_pert_table <- def_bin_X_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- "binomial"
type <- "deviance"

# sim1: Binomial IV for fixed v ---

set.seed(35732)

nsim <- 100

xi_values <- seq(0, 50, by = 1)
xi_values <- c(xi_values, 10000)

# Simulate
sim_data_bin_X_fivi <- simulate_fivi(nsim = nsim, nobs = 1000,
                                     xi_values = xi_values,
                                     data_table = data_table,
                                     data_pert_table = data_pert_table, 
                                     formula = formula, A_formula = A_formula,
                                     family = family, type = type,
                                     quant_value = 0.9) 

# # Save simulated data list
# path_name <- "./data sets/simulation_study/ex3/"
# dir.create(paste(path_name, Sys.Date(), sep = ""))
# save(sim_data_bin_X_fivi, file = paste(paste(path_name, Sys.Date(), sep = ""), 
#                                        "/sim1.Rdata", sep =""))

# sim2: Binomial IV for fixed xi ---

start_seed <- 19336

nsim <- 100

xi_values <- c(0, 1, 2, 5, 10, 20, 50, 100, 10000)
v_values <- seq(-12, 12, by = 0.1)

# Simulate
sim_data_bin_X_fixi <- simulate_fixi(nsim = nsim, nobs = 1000,
                                     xi_values = xi_values,
                                     v_values = v_values, 
                                     data_table = data_table,
                                     data_pert_table = data_pert_table, 
                                     formula = formula, A_formula = A_formula,
                                     family = family, type = type,
                                     quant_value = 0.9, start_seed)
 
# # Save simulated data list
# path_name <- "./data sets/simulation_study/ex3/"
# dir.create(paste(path_name, Sys.Date(), sep = ""))
# save(sim_data_bin_X_fixi, file = paste(paste(path_name, Sys.Date(), sep = ""), 
#                                        "/sim2.Rdata", sep =""))

# sim3: IV investigation for binomial ---

set.seed(58435)

nsim <- 100


nobs_values <- seq(100, 5000, by = 100)

#nobs_values <- c(100, 1000, 5000, 10000)

# xi_values <- seq(0, 100, by = 10)
xi_values <- c(0, 1, 3, 10, 10000)

# Simulate
IV_b_bin <- simulate_IV(nsim = nsim, nobs_values = nobs_values,
                        xi_values = xi_values,
                        data_table = data_table,
                        formula = formula, A_formula = A_formula,
                        family = family, type = type)

IV_b_data <- IV_b_bin$sim_data
plot_IV(IV_b_data, causal_parameter = b)

# # Save simulated data list
# path_name <- "./data sets/simulation_study/ex3/"
# dir.create(paste(path_name, Sys.Date(), sep = ""))
# save(IV_b_bin, file = paste(paste(path_name, Sys.Date(), sep = ""), 
#                             "/sim3.Rdata", sep =""))
# 




IV_b_bin <- simulate_IV_mult(nsim = nsim, nobs_values = nobs_values,
                             xi_values = xi_values,
                             data_table = data_table,
                             formula = formula, A_formula = A_formula,
                             family = family, type = type)
IV_b_data <- IV_b_bin$sim_data
plot_IV_mult(IV_b_data,
             causal_intercept = c,
             causal_parameter = b)

# sim4: Binomial IV several quantiles ---

set.seed(333381)

nsim <- 100

xi_values <- seq(0, 10, by = 0.1)
xi_values <- xi_values[-1]
xi_big <- 10000

# Simulate
sim_data_bin_X_quant <- simulate_fivi_quant(nsim = nsim, nobs = 1000,
                                            xi_values = xi_values,
                                            xi_big = xi_big,
                                            data_table = data_table,
                                            data_pert_table = data_pert_table, 
                                            formula = formula,
                                            A_formula = A_formula,
                                            family = family, type = type) 

# # Save data
# path_name <- "./data sets/simulation_study/ex3/"
# dir.create(paste(path_name, Sys.Date(), sep = ""))
# save(sim_data_bin_X_quant, file = paste(paste(path_name, Sys.Date(), sep = ""),
#                                         "/sim4.Rdata", sep =""))
# 
# 


# ex4: Poisson with Anchor on X, H and Y --------------------------------------

# Define variables for unperturbed and perturbed data set
def_poi_XHY <- defData(varname = "A", dist = "normal", 
                     formula = 0, variance = 0.25)
def_poi_XHY <- defData(def_poi_XHY, varname = "H", dist = "normal",
                     formula = "0.5 * A", variance = 0.25)
def_poi_XHY <- defData(def_poi_XHY, varname = "X", dist = "normal", 
                     formula = "2 * H + A", variance = 0.25)
def_poi_XHY <- defData(def_poi_XHY, varname = "Y", dist = "poisson",
                       link = "log", 
                     formula = "0.4 * X + 1 * H - A", variance = 1)

def_poi_XHY_pert <- defData(varname = "A", dist = "normal", 
                          formula = 0, variance = 0.25)
def_poi_XHY_pert <- defData(def_poi_XHY_pert, varname = "v", 
                          formula = 2) # set perturbation strength
def_poi_XHY_pert <- defData(def_poi_XHY_pert, varname = "H", dist = "normal",
                          formula = "0.5", variance = 0.25)
def_poi_XHY_pert <- defData(def_poi_XHY_pert, varname = "X", dist = "normal", 
                          formula = "2 * H + v", variance = 0.25)
def_poi_XHY_pert <- defData(def_poi_XHY_pert, varname = "Y", dist = "poisson",
                            link = "log",
                          formula = "0.4 * X + 1 * H - 1.5", variance = 1)

# Parameters of linear predictor
b <- 0.4
h <- 1
a <- -1

a_pert <- -1.5

# Distribution histogram
dd <- genData(300, def_poi_XHY)
dd_pert <- genData(300, def_poi_XHY_pert)
par(mfrow = c(2,3))

hist(b*dd$X + h * dd$H + a * dd$A, breaks = 100)
hist(exp(b*dd$X + h * dd$H + a * dd$A), breaks = 100)
hist(dd$Y, breaks = 100)

hist(b*dd_pert$X + h * dd_pert$H + a_pert, breaks = 100)
hist(exp(b*dd_pert$X + h * dd_pert$H + a_pert), breaks = 100)
hist(dd_pert$Y, breaks = 100)

# Initialize
data_table <- def_poi_XHY
data_pert_table <- def_poi_XHY_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- "poisson"
type <- "deviance"

# sim1: Poisson XHY for fixed v ---

set.seed(22416)

nsim <- 100

xi_values <- seq(0, 50, by = 0.1)
xi_values <- c(xi_values, 10000)


# Simulate
sim_data_poi_XHY_fivi <- simulate_fivi(nsim = nsim, nobs = 1000,
                                       xi_values = xi_values,
                                       data_table = data_table,
                                       data_pert_table = data_pert_table, 
                                       formula = formula, A_formula = A_formula,
                                       family = family, type = type,
                                       quant_value = 0.9)

# # Save simulated data list
# path_name <- "./data sets/simulation_study/ex4/"
# dir.create(paste(path_name, Sys.Date(), sep = ""))
# save(sim_data_poi_XHY_fivi, file = paste(paste(path_name, Sys.Date(), sep = ""),
#                                          "/sim1.Rdata", sep =""))

# sim4: Poisson IV several quantiles ---

set.seed(93671)

nsim <- 100

xi_values <- seq(0, 10, by = 0.1)
xi_values <- xi_values[-1]
xi_big <- 10000

# Simulate
sim_data_poi_XHY_quant <- simulate_fivi_quant(nsim = nsim, nobs = 1000,
                                              xi_values = xi_values,
                                              xi_big = xi_big,
                                              data_table = data_table,
                                              data_pert_table = data_pert_table, 
                                              formula = formula,
                                              A_formula = A_formula,
                                              family = family, type = type) 

# # Save data
# path_name <- "./data sets/simulation_study/ex4/"
# dir.create(paste(path_name, Sys.Date(), sep = ""))
# save(sim_data_poi_XHY_quant, file = paste(paste(path_name, Sys.Date(), sep = ""),
#                                           "/sim4.Rdata", sep =""))


# ex5: Label Noise Example ----------------------------------------------------

# Define variables for unperturbed and perturbed data set
def_LN <- defData(varname = "A", dist = "normal", 
                  formula = 0, variance = 0.01)
def_LN <- defData(def_LN, varname = "B", dist = "binomial", link = "logit",
                  formula = "-2 + A", variance = 1)
def_LN <- defData(def_LN, varname = "H", dist = "normal",
                  formula = 0, variance = 0.01)
def_LN <- defData(def_LN, varname = "X", dist = "normal", 
                  formula = "2 * H", variance = 0.01)
def_LN <- defData(def_LN, varname = "Yorg", dist = "binomial", link = "logit", 
                  formula = "-1 + 0.4 * X + 1 * H", variance = 1)
def_LN <- defData(def_LN, varname = "Y", dist = "binary", 
                  formula = "abs(Yorg - B)")

def_LN_pert <- defData(varname = "A", dist = "normal", 
                       formula = 0, variance = 0.01)
def_LN_pert <- defData(def_LN_pert, varname = "v", 
                       formula = 0.8) # set perturbation strength
def_LN_pert <- defData(def_LN_pert, varname = "B", dist = "binomial",
                       link = "logit",
                       formula = "-2 + v", variance = 1)
def_LN_pert <- defData(def_LN_pert, varname = "H", dist = "normal",
                       formula = 0, variance = 0.01)
def_LN_pert <- defData(def_LN_pert, varname = "X", dist = "normal", 
                       formula = "2 * H", variance = 0.01)
def_LN_pert <- defData(def_LN_pert, varname = "Yorg", dist = "binomial", 
                       link = "logit", 
                       formula = "-1 + 0.4 * X + 1 * H", variance = 1)
def_LN_pert <- defData(def_LN_pert, varname = "Y", dist = "binary", 
                       formula = "abs(Yorg - B)")

# Parameters of linear predictor
b <- 0.4
h <- 1
c <- -1

t <- -2
v_pert <- 0.8

# Distribution histogram
dd <- genData(1000, def_LN)
dd_pert <- genData(1000, def_LN_pert)

eta <- c + b * dd$X + h * dd$H
eta_pert <- c + b * dd_pert$X + h * dd_pert$H

tau <- t + dd$A
tau_pert <- t + v_pert

par(mfrow = c(3,3))
hist(tau, breaks = 100)
hist(exp(tau)/(1+exp(tau)), breaks = 100)
hist(dd$B, breaks = 100)
hist(eta, breaks = 100)
hist(exp(eta)/(1+exp(eta)), breaks = 100)
hist(dd$Yorg, breaks = 100)
hist(dd$Y, breaks = 100)

par(mfrow = c(3,3))
hist(tau_pert, breaks = 100)
hist(exp(tau_pert)/(1+exp(tau_pert)), breaks = 100)
hist(dd_pert$B, breaks = 100)
hist(eta_pert, breaks = 100)
hist(exp(eta_pert)/(1+exp(eta_pert)), breaks = 100)
hist(dd_pert$Yorg, breaks = 100)
hist(dd_pert$Y, breaks = 100)

exp(tau_pert)/(1+exp(tau_pert))

# Initialize
data_table <- def_LN
data_pert_table <- def_LN_pert
formula <- Y ~ X
A_formula <- ~ A
family <- "binomial"
type <- "deviance"

# sim1: Label Noise for fixed v ---

set.seed(665401)

nsim <- 10

xi_values <- seq(0, 100, by = 1)
xi_values <- c(xi_values, 10000)

# Simulate
sim_data_LN_fivi <- simulate_fivi(nsim = nsim, nobs = 100,
                                  xi_values = xi_values,
                                  data_table = data_table,
                                  data_pert_table = data_pert_table, 
                                  formula = formula, A_formula = A_formula,
                                  family = family, type = type,
                                  quant_value = 0.9) 

# # Save simulated data list
# path_name <- "./data sets/simulation_study/ex5/"
# dir.create(paste(path_name, Sys.Date(), sep = ""))
# save(sim_data_LN_fivi, file = paste(paste(path_name, Sys.Date(), sep = ""), 
#                                     "/sim1.Rdata", sep =""))

# sim2: Label Noise for fixed xi ---

start_seed <- 78957

nsim <- 10

xi_values <- c(0, 1, 3, 5, 8, 10, 20, 50, 10000)
v_values <- seq(-10, 10, by = 1)

# Simulate
sim_data_LN_fixi <- simulate_fixi(nsim = nsim, nobs = 100,
                                  xi_values = xi_values,
                                  v_values = v_values, 
                                  data_table = data_table,
                                  data_pert_table = data_pert_table, 
                                  formula = formula, A_formula = A_formula,
                                  family = family, type = type,
                                  quant_value = 0.9, start_seed)

# # Save simulated data list
# path_name <- "./data sets/simulation_study/ex5/"
# dir.create(paste(path_name, Sys.Date(), sep = ""))
# save(sim_data_LN_fixi, file = paste(paste(path_name, Sys.Date(), sep = ""), 
#                                     "/sim2.Rdata", sep =""))

sim_data_LN_fixi_glm <- simulate_fixi(nsim = nsim, nobs = 1000,
                                      xi_values = 0,
                                      v_values = v_values, 
                                      data_table = data_table,
                                      data_pert_table = data_pert_table, 
                                      formula = Y ~ X + A - 1,
                                      A_formula = A_formula,
                                      family = family, type = type,
                                      quant_value = 0.9, start_seed)

# # Save simulated data list
# path_name <- "./data sets/simulation_study/ex5/"
# dir.create(paste(path_name, Sys.Date(), sep = ""))
# save(sim_data_LN_fixi_glm, file = paste(paste(path_name, Sys.Date(), sep = ""), 
#                                         "/sim2_glm.Rdata", sep =""))



# The following is code for tryout settings for LN with anchor only on tau
# # sim4: Label Noise several quantiles ---
# 
# set.seed(16541)
# 
# nsim <- 100
# 
# # xi_values <- c(0, 0.3, 0.5, 0.8, seq(1, 100, by = 5))
# 
# xi_values <- seq(0, 10, by = 0.1)
# 
# xi_big <- 10000
# 
# # Simulate for evaluation Y with noise
# sim_data_LN_quant <- simulate_fivi_quant(nsim = nsim, nobs = 1000,
#                                          xi_values = xi_values, xi_big = xi_big,
#                                          data_table = data_table,
#                                          data_pert_table = data_pert_table, 
#                                          formula = formula,
#                                          A_formula = A_formula,
#                                          family = family, type = type) 
# 
# # # Save data
# # path_name <- "./data sets/simulation_study/ex5/"
# # dir.create(paste(path_name, Sys.Date(), sep = ""))
# # save(sim_data_LN_quant, file = paste(paste(path_name, Sys.Date(), sep = ""),
# #                                      "/sim4.Rdata", sep =""))
# 
# data_LN_quant <- sim_data_LN_quant$sim_data
# q_values <- seq(0, 1, by = 0.01)
# plot_quant(data_LN_quant, q_values, xi_big = 10000)
# 
# # Prediction on Y_org
# data_pert_table_raw <- updateDef(data_pert_table,
#                                  changevar = "Y",
#                                  newformula = "Yorg")
# set.seed(44652)
# # Simulate
# sim_data_LN_quant_raw <- simulate_fivi_quant(nsim = nsim, nobs = 1000,
#                                              xi_values = xi_values,
#                                              xi_big = xi_big,
#                                              data_table = data_table,
#                                              data_pert_table = data_pert_table_raw, 
#                                              formula = formula,
#                                              A_formula = A_formula,
#                                              family = family, type = type) 
# # # Save data
# # path_name <- "./data sets/simulation_study/ex5/"
# # dir.create(paste(path_name, Sys.Date(), sep = ""))
# # save(sim_data_LN_quant_raw, file = paste(paste(path_name, Sys.Date(), sep = ""), 
# #                                          "/sim4_raw.Rdata", sep =""))
# 
# data_LN_quant_raw <- sim_data_LN_quant_raw$sim_data
# 
# # head(data_LN_quant_raw)
# # str(data_LN_quant_raw)
# 
# q_values <- seq(0, 1, by = 0.01)
# plot_quant(data_LN_quant_raw, q_values, xi_big = 10000)
# 
# 
# # # Save data
# # path_name <- "./data sets/simulation_study/ex5/"
# # dir.create(paste(path_name, Sys.Date(), sep = ""))
# # save(sim_data_LN_quant, file = paste(paste(path_name, Sys.Date(), sep = ""), "/sim4.Rdata", sep =""))
# # 
# set.seed(55612)
# 
# # Simulation for glm(Y ~ X + A - 1)
# sim_data_LN_quant_glm <- simulate_fivi_quant_glm(nsim = nsim, nobs = 1000,
#                                                  xi_values = 0, xi_big = xi_big,
#                                                  data_table = data_table,
#                                                  data_pert_table = data_pert_table, 
#                                                  formula =  Y ~ X + A - 1,
#                                                  A_formula = A_formula,
#                                                  family = family, type = type) 
# 
# # # Save data
# # path_name <- "./data sets/simulation_study/ex5/"
# # dir.create(paste(path_name, Sys.Date(), sep = ""))
# # save(sim_data_LN_quant_glm, file = paste(paste(path_name, Sys.Date(), sep = ""), "/sim4_glm.Rdata", sep =""))
# 
# data_LN_quant <- sim_data_LN_quant$sim_data
# 
# head(data_LN_quant)
# str(data_LN_quant)
# 
# q_values <- seq(0, 1, by = 0.01)
# #plot_quant(data_LN_quant, q_values, xi_big = 10000)
# 
# # Plot with glm(Y~X+A-1)
# data_LN_quant_glm <- sim_data_LN_quant_glm$sim_data
# 
# q_len <- length(q_values)
# 
# quantile_array <- matrix(nrow = nsim, ncol = q_len)
# 
# for (i in 1:nsim) {
#   for (q in 1:q_len) {
#     quant <- q_values[q]
#     
#     quantile_array[i, q] <- 
#       as.numeric(quantile(data_LN_quant_glm[[i]][[2]][, 2], quant))
#     
#   }
# }
# glm_quantile_vector <- apply(quantile_array, 2, mean)
# 
# plot_quant(data_LN_quant, q_values, xi_big = 10000, glm_quantile_vector)
# 
# 
# # Prediction on Y without noise
# data_pert_table_raw <- updateDef(data_pert_table,
#                                  changevar = "Y",
#                                  newformula = "Yorg")
# 
# set.seed(44652)
# 
# nsim <- 100
# 
# xi_values <- seq(0, 10, by = 0.1)
# xi_values <- xi_values[-1]
# xi_big <- 10000
# 
# 
# # Simulate
# sim_data_LN_quant_raw <- simulate_fivi_quant(nsim = nsim, nobs = 1000,
#                                              xi_values = xi_values, xi_big = xi_big,
#                                              data_table = data_table,
#                                              data_pert_table = data_pert_table_raw, 
#                                              formula = formula, A_formula = A_formula,
#                                              family = family, type = type) 
# 
# # # Save data
# # path_name <- "./data sets/simulation_study/ex5/"
# # dir.create(paste(path_name, Sys.Date(), sep = ""))
# # save(sim_data_LN_quant_raw, file = paste(paste(path_name, Sys.Date(), sep = ""), "/sim4_raw.Rdata", sep =""))
# 
# set.seed(55612)
# 
# # Simulate for glm predict on Y without noise
# sim_data_LN_quant_glm_raw <- simulate_fivi_quant_glm(nsim = nsim, nobs = 1000,
#                                                      xi_values = 0, xi_big = xi_big,
#                                                      data_table = data_table,
#                                                      data_pert_table = data_pert_table_raw, 
#                                                      formula =  Y ~ X + A - 1,
#                                                      A_formula = A_formula,
#                                                      family = family, type = type) 
# 
# # # Save data
# # path_name <- "./data sets/simulation_study/ex5/"
# # dir.create(paste(path_name, Sys.Date(), sep = ""))
# # save(sim_data_LN_quant_glm_raw, file = paste(paste(path_name, Sys.Date(), sep = ""), "/sim4_glm_raw.Rdata", sep =""))
# 
# data_LN_quant_raw <- sim_data_LN_quant_raw$sim_data
# 
# head(data_LN_quant_raw)
# str(data_LN_quant_raw)
# 
# q_values <- seq(0, 1, by = 0.01)
# #plot_quant(data_LN_quant_raw, q_values, xi_big = 10000)
# 
# # Plot with glm(Y~X+A-1)
# data_LN_quant_glm_raw <- sim_data_LN_quant_glm_raw$sim_data
# 
# q_len <- length(q_values)
# 
# quantile_array <- matrix(nrow = nsim, ncol = q_len)
# 
# for (i in 1:nsim) {
#   for (q in 1:q_len) {
#     quant <- q_values[q]
#     
#     quantile_array[i, q] <- 
#       as.numeric(quantile(data_LN_quant_glm_raw[[i]][[2]][, 2], quant))
#     
#   }
# }
# glm_quantile_vector_raw <- apply(quantile_array, 2, mean)
# 
# plot_quant(data_LN_quant_raw, q_values, xi_big = 10000, glm_quantile_vector_raw)

# ex6: Label Noise Example Involved -------------------------------------------

# Define variables for unperturbed and perturbed data set
def_LN <- defData(varname = "A", dist = "normal", 
                  formula = 0, variance = 0.01)
def_LN <- defData(def_LN, varname = "B", dist = "binomial", link = "logit",
                  formula = "-2 + A", variance = 1)
def_LN <- defData(def_LN, varname = "H", dist = "normal",
                  formula = "A", variance = 0.01)
def_LN <- defData(def_LN, varname = "X", dist = "normal", 
                  formula = "2 * H + A", variance = 0.01)
def_LN <- defData(def_LN, varname = "Yorg", dist = "binomial", link = "logit", 
                  formula = "-1 + 0.4 * X + 1 * H", variance = 1)
def_LN <- defData(def_LN, varname = "Y", dist = "binary", 
                  formula = "abs(Yorg - B)")

def_LN_pert <- defData(varname = "A", dist = "normal", 
                       formula = 0, variance = 0.01)
def_LN_pert <- defData(def_LN_pert, varname = "v", 
                       formula = 0.6) # set perturbation strength
def_LN_pert <- defData(def_LN_pert, varname = "B", dist = "binomial",
                       link = "logit",
                       formula = "-2 + v", variance = 1)
def_LN_pert <- defData(def_LN_pert, varname = "vH", 
                       formula = 0.4) # set perturbation strength
def_LN_pert <- defData(def_LN_pert, varname = "H", dist = "normal",
                       formula = "vH", variance = 0.01)
def_LN_pert <- defData(def_LN_pert, varname = "vX", 
                       formula = -0.5) # set perturbation strength
def_LN_pert <- defData(def_LN_pert, varname = "X", dist = "normal", 
                       formula = "2 * H + vX", variance = 0.01)
def_LN_pert <- defData(def_LN_pert, varname = "Yorg", dist = "binomial",
                       link = "logit", 
                       formula = "-1 + 0.4 * X + 1 * H", variance = 1)
def_LN_pert <- defData(def_LN_pert, varname = "Y", dist = "binary", 
                       formula = "abs(Yorg - B)")

# Parameters of linear predictor
b <- 0.4
h <- 1
c <- -1

t <- -2
v_pert <- 0.6

# Distribution histogram
dd <- genData(1000, def_LN)
dd_pert <- genData(1000, def_LN_pert)

eta <- c + b * dd$X + h * dd$H
eta_pert <- c + b * dd_pert$X + h * dd_pert$H

tau <- t + dd$A
tau_pert <- t + v_pert

par(mfrow = c(3,3))
hist(tau, breaks = 100)
hist(exp(tau)/(1+exp(tau)), breaks = 100)
hist(dd$B, breaks = 100)
hist(eta, breaks = 100)
hist(exp(eta)/(1+exp(eta)), breaks = 100)
hist(dd$Yorg, breaks = 100)
hist(dd$Y, breaks = 100)

par(mfrow = c(3,3))
hist(tau_pert, breaks = 100)
hist(exp(tau_pert)/(1+exp(tau_pert)), breaks = 100)
hist(dd_pert$B, breaks = 100)
hist(eta_pert, breaks = 100)
hist(exp(eta_pert)/(1+exp(eta_pert)), breaks = 100)
hist(dd_pert$Yorg, breaks = 100)
hist(dd_pert$Y, breaks = 100)

exp(tau_pert)/(1+exp(tau_pert))

# Initialize
data_table <- def_LN
data_pert_table <- def_LN_pert
formula <- Y ~ X
A_formula <- ~ A
family <- "binomial"
type <- "deviance"

# sim1: Label Noise for fixed v ---

set.seed(93441)

nsim <- 100

xi_values <- seq(0, 10, by = 0.1)
xi_values <- c(xi_values, 10000)

# Simulate
sim_data_LN_fivi <- simulate_fivi(nsim = nsim, nobs = 1000,
                                  xi_values = xi_values,
                                  data_table = data_table,
                                  data_pert_table = data_pert_table, 
                                  formula = formula, A_formula = A_formula,
                                  family = family, type = type,
                                  quant_value = 0.9) 

# # Save simulated data list
# path_name <- "./data sets/simulation_study/ex6/"
# dir.create(paste(path_name, Sys.Date(), sep = ""))
# save(sim_data_LN_fivi, file = paste(paste(path_name, Sys.Date(), sep = ""),
#                                     "/sim1.Rdata", sep =""))

# sim4: Label Noise several quantiles ---

set.seed(321459)

nsim <- 100

# xi_values <- c(0, 0.3, 0.5, 0.8, seq(1, 100, by = 5))

xi_values <- seq(0, 10, by = 0.1)

xi_big <- 10000

# Simulate
sim_data_LN_quant <- simulate_fivi_quant(nsim = nsim, nobs = 1000,
                                         xi_values = xi_values,
                                         xi_big = xi_big,
                                         data_table = data_table,
                                         data_pert_table = data_pert_table, 
                                         formula = formula,
                                         A_formula = A_formula,
                                         family = family, type = type) 

# # Save data
# path_name <- "./data sets/simulation_study/ex6/"
# dir.create(paste(path_name, Sys.Date(), sep = ""))
# save(sim_data_LN_quant, file = paste(paste(path_name, Sys.Date(), sep = ""),
#                                      "/sim4.Rdata", sep =""))

# Prediction on Y_org
data_pert_table_raw <- updateDef(data_pert_table,
                                 changevar = "Y",
                                 newformula = "Yorg")
set.seed(120938)

# Simulate
sim_data_LN_quant_raw <- simulate_fivi_quant(nsim = nsim, nobs = 1000,
                                             xi_values = xi_values,
                                             xi_big = xi_big,
                                             data_table = data_table,
                                             data_pert_table = data_pert_table_raw, 
                                             formula = formula,
                                             A_formula = A_formula,
                                             family = family, type = type) 
# # Save data
# path_name <- "./data sets/simulation_study/ex6/"
# dir.create(paste(path_name, Sys.Date(), sep = ""))
# save(sim_data_LN_quant_raw, file = paste(paste(path_name, Sys.Date(), sep = ""), 
#                                          "/sim4_raw.Rdata", sep =""))

data_LN_quant_raw <- sim_data_LN_quant_raw$sim_data
