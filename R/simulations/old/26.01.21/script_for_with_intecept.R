# Function definition for IV parameter identifiability ------------------------

onerep_IV_mult <- function(rep, nobs_values = 300, data_table, 
                           formula, A_formula, xi_values,
                           family, type, dim = 2) {
  
  data <- as.data.frame(matrix(ncol = 3 + (2 * dim), nrow = 0))
  colnames(data) <- c("rep", "nobs", "xi", "intercept", "int_se", "b", "b_se")
  
  for (k in 1:length(nobs_values)) {
    
    nobs <- nobs_values[k]
    
    # Generate data set with nobs observations
    dd <- genData(nobs, data_table)
    # dd_pert <- genData(nobs, data_pert_table)
    
    xi_len <- length(xi_values)
    
    b <- matrix(nrow = xi_len, ncol = dim)
    b_se <- matrix(nrow = xi_len, ncol = dim)
    
    # logLik_pert <- numeric(xi_len)
    
    for (i in 1:xi_len) {
      
      xi <- xi_values[i]
      
      fit_temp <- glare(formula = formula,
                        A_formula = A_formula,
                        data = dd,
                        xi = xi,
                        family = family,
                        type = type)
      
      b[i, ] <- as.numeric(coef(fit_temp))
      b_se[i, ] <- fit_temp$coef_se
      
      # logLik_pert_indiv <- logLik(fit_temp, newdata = dd_pert, indiv = TRUE)
      # logLik_pert[i] <-
      #   as.numeric(quantile(logLik_pert_indiv, quant_value))
    }
    
    # Return
    data_temp <- data.frame(rep = rep, nobs = nobs,
                            xi = xi_values,
                            intercept = b[, 1],
                            int_se = b_se[, 1],
                            b = b[, 2],
                            b_se = b_se[, 2])
    
    data <- rbind(data, data_temp)
  }
  
  data
}

# Define simulation function
simulate_IV_mult <- function(nsim, nobs_values = 300, xi_values,
                             data_table, formula, A_formula, family, type,
                             dim =2) {
  
  xi_len <- length(xi_values)
  
  states <- array(dim = c(nsim, 626))
  
  b_data <- as.data.frame(matrix(ncol = 3 + (2 * dim), nrow = 0))
  colnames(b_data) <- c("rep", "nobs", "xi", "intercept", "int_se", "b", "b_se")
  
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  for (r in 1:nsim) {
    states[r, ] <- .Random.seed
    
    data_temp <- onerep_IV_mult(rep = r, nobs_values = nobs_values,
                                data_table = data_table, 
                                formula = formula, A_formula = A_formula,
                                xi_values = xi_values,
                                family = family, type = type, dim = dim)
    
    b_data <- rbind(b_data, data_temp)
    
    Sys.sleep(0.1)
    setTxtProgressBar(pb, r)
  }
  close(pb)
  
  list(states = states, sim_data = b_data)
}


# Plot IV function ------------------------------------------------------------

plot_IV_mult <- function(data,
                         causal_intercept = 1,
                         causal_parameter = 0.4,
                         xi_big = 10000) {
  
  # Average over simulations
  mean_data <- data %>%
    group_by(xi, nobs) %>%
    summarise(mean_intercept = mean(intercept), se_int = sd(intercept),
              mean_b = mean(b), se_b = sd(b))
  
  # Calculate CI
  upper_int <- mean_data$mean_intercept +
    qt(0.975, 100 - 1) * mean_data$se_int / sqrt(100)
  lower_int <- mean_data$mean_intercept -
    qt(0.975, 100 - 1) * mean_data$se_int / sqrt(100)
  mean_data$upper_int <- upper_int
  mean_data$lower_int <- lower_int
  upper_b <- mean_data$mean_b +
    qt(0.975, 100 - 1) * mean_data$se_b / sqrt(100)
  lower_b <- mean_data$mean_b -
    qt(0.975, 100 - 1) * mean_data$se_b / sqrt(100)
  mean_data$upper_b <- upper_b
  mean_data$lower_b <- lower_b
  
  # Intercept
  gg_data1 <- mean_data
  gg_data1$xi <- as.factor(gg_data1$xi)
  gg1 <- ggplot(data = gg_data1, aes(y = mean_intercept, x = nobs, group = xi, color = xi)) +
    geom_point() +
    geom_line() + 
    geom_errorbar(aes(x = nobs, ymin = lower_int, ymax = upper_int)) +
    
    geom_hline(yintercept = causal_intercept, linetype = "dashed") +
    
    xlab("Number of Observation") +
    ylab("Intercept Estimate")
  
  # Causal Parameter
  gg_data2 <- gg_data1
  gg2 <- ggplot(data = gg_data2, aes(y = mean_b, x = nobs, group = xi, color = xi)) +
    geom_point() +
    geom_line() + 
    geom_errorbar(aes(x = nobs, ymin = lower_b, ymax = upper_b)) +
    
    geom_hline(yintercept = causal_parameter, linetype = "dashed") +
    
    xlab("Number of Observation") +
    ylab("Causal Parameter Estimate")
  
  
  # nobs 100-5000  
  # gg_data_big2 <- gg_data_big[gg_data_big$nobs >= 100 & gg_data_big$nobs <= 5000, ]
  
  
  print(gg1)
  print(gg2)
  
  invisible(mean_data)
}
