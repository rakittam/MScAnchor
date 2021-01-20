# Define plot functions
#
# Date: 14.12.20     Author: Maic Rakitta
###############################################################################

#library(extrafont)
library(tidyverse)
library(rsimsum)
library(ggplot2)
# library(gridExtra)
library(reshape2)
library(ggpubr) # for theme and ggarrange
theme_set(theme_bw())

# -----------------------------------------------------------------------------
### Plot Function Definitions
# Plot function for anchor regression comparison ------------------------------
plot_rot <- function(sim_data, xi_big = 10000) {
  
  mean_data <- sim_data %>%
    group_by(xi, gamma) %>%
    summarise(mean_b_ar = mean(b_ar), mean_se_ar = mean(se_ar),
              mean_MSE_pert_ar = mean(MSE_pert_ar),
              mean_b_dev = mean(b_dev), mean_se_dev = mean(se_dev),
              mean_MSE_pert_dev = mean(MSE_pert_dev),
              mean_b_pea = mean(b_pea), mean_se_pea = mean(se_pea),
              mean_MSE_pert_pea = mean(MSE_pert_pea))
  
  # Plot parameter comparison
  gg_parameter_data <- data.frame(x = c(mean_data$mean_b_ar[!is.na(mean_data$mean_b_dev)],
                                        mean_data$mean_b_ar[!is.na(mean_data$mean_b_pea)]),
                                  y = c(mean_data$mean_b_dev[!is.na(mean_data$mean_b_dev)],
                                        mean_data$mean_b_pea[!is.na(mean_data$mean_b_pea)]),
                                  Method = c(rep("DEV", length(mean_data$mean_b_dev[!is.na(mean_data$mean_b_dev)])),
                                             rep("PEA", length(mean_data$mean_b_pea[!is.na(mean_data$mean_b_pea)]))))
  
  gg_parameter <- ggplot(gg_parameter_data, aes(y = y, x = x, color = Method, group = Method)) +
    
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype ="dashed") +
    geom_vline(xintercept = mean_data$mean_b_ar[mean_data$xi == 0], linetype = "dashed", col = "gray") +
    geom_vline(xintercept = mean_data$mean_b_ar[which.max(mean_data$xi)], linetype = "dashed", col = "gray") +
    
    xlim(mean_data$mean_b_ar[which.max(mean_data$xi)],
         mean_data$mean_b_ar[mean_data$xi == 0]) +
    
    ylab("GLARE parameters") +
    xlab("AR parameter") 
  
  # Plot ex 2 from Rothenhaeusler
  b_OLS_ar <- mean_data$mean_b_ar[mean_data$gamma == 1]
  b_IV_ar <- mean_data$mean_b_ar[mean_data$xi == xi_big]
  b_PA_ar <- mean_data$mean_b_ar[mean_data$gamma == 0]
  
  MSE_OLS_ar <- mean_data$mean_MSE_pert_ar[mean_data$gamma == 1]
  MSE_IV_ar <- mean_data$mean_MSE_pert_ar[mean_data$xi == xi_big]
  MSE_PA_ar <- mean_data$mean_MSE_pert_ar[mean_data$gamma == 0]
  
  gg_data <- data.frame(x = c(mean_data$mean_b_ar,
                              mean_data$mean_b_dev[!is.na(mean_data$mean_b_dev)],
                              mean_data$mean_b_pea[!is.na(mean_data$mean_b_pea)]),
                        y = c(mean_data$mean_MSE_pert_ar,
                              mean_data$mean_MSE_pert_dev[!is.na(mean_data$mean_b_dev)],
                              mean_data$mean_MSE_pert_pea[!is.na(mean_data$mean_b_pea)]),
                        Method = c(rep("AR", length(mean_data$mean_b_ar)),
                                   rep("DEV", length(mean_data$mean_b_dev[!is.na(mean_data$mean_b_dev)])),
                                   rep("PEA", length(mean_data$mean_b_pea[!is.na(mean_data$mean_b_pea)]))))
  
  # Use this data if only AR vs DEV wanted
  gg_data2 <- data.frame(x = c(mean_data$mean_b_ar,
                               mean_data$mean_b_dev[!is.na(mean_data$mean_b_dev)]),
                         y = c(mean_data$mean_MSE_pert_ar,
                               mean_data$mean_MSE_pert_dev[!is.na(mean_data$mean_b_dev)]),
                         Method = c(rep("AR", length(mean_data$mean_b_ar)),
                                    rep("DEV", length(mean_data$mean_b_dev[!is.na(mean_data$mean_b_dev)]))))
  
  gg <- ggplot(gg_data2, aes(y = y, x = x, color = Method, group = Method)) +
    
    geom_line(aes(linetype = Method)) +
    
    ylab("MSE of perturbed data") +
    
    annotate("point", colour = 1, x = b_OLS_ar, y = MSE_OLS_ar) +
    annotate("text", label = "OLS", x = b_OLS_ar - 0.02, y = MSE_OLS_ar + 0.1) +
    annotate("point", colour = 1, x = b_IV_ar, y = MSE_IV_ar) +
    annotate("text", label = "IV", x = b_IV_ar + 0.02, y = MSE_IV_ar + 0.1) +
    annotate("point", colour = 1, x = b_PA_ar, y = MSE_PA_ar) +
    annotate("text", label = "PA", x = b_PA_ar - 0.02, y = MSE_PA_ar + 0.1) +
    
    scale_x_continuous(
      
      # Features of the first axis
      name = "b",
      
      # Add a second axis and specify its features
      sec.axis = sec_axis(trans= ~.,
                          name=expression(gamma~(xi)),
                          breaks=c(1.111111, 1.555556, 2.000000),
                          labels=c("17.7 (8.3)", "1.8 (0.4)", "0 (-0.5)"))
    )
  
  print(gg_parameter)
  print(gg)
  invisible(mean_data)
}

# Plot function for fixed v interventions for Rothenhaeusler Example ----------
plot_rot_fivi <- function(sim_data, xi_big = 10000) {
  
  sim_data_glare <- sim_data[sim_data$xi != xi_big, ]
  sim_data_big <- sim_data[sim_data$xi == xi_big, ]
  
  mean_data_glare <- sim_data_glare %>%
    group_by(xi) %>%
    summarise(mean_b_dev = mean(b_dev), mean_se_dev = mean(se_dev),
              mean_logLike_pert_dev = mean(logLik_pert_dev),
              mean_b_pea = mean(b_pea), mean_se_pea = mean(se_pea),
              mean_logLike_pert_pea = mean(logLik_pert_pea))
  
  mean_data_big <- sim_data_big %>%
    group_by(xi) %>%
    summarise(mean_b_dev = mean(b_dev), mean_se_dev = mean(se_dev),
              mean_logLike_pert_dev = mean(logLik_pert_dev),
              mean_b_pea = mean(b_pea), mean_se_pea = mean(se_pea),
              mean_logLike_pert_pea = mean(logLik_pert_pea))
  
  gg_data <- data.frame(x = c(mean_data_glare$mean_b_dev,
                              mean_data_glare$mean_b_pea),
                        y = c(mean_data_glare$mean_logLike_pert_dev,
                              mean_data_glare$mean_logLike_pert_pea),
                        xi = c(mean_data_glare$xi, mean_data_glare$xi),
                        Method = c(rep("DEV", length(mean_data_glare$mean_b_dev)),
                                   rep("PEA", length(mean_data_glare$mean_b_pea))))
  
  gg <- ggplot(gg_data, aes(y = y, x = xi, color = Method, group = Method)) +
    
    geom_line() +
    geom_hline(yintercept = mean_data_big$mean_logLike_pert_dev, linetype = "dashed") +
    
    ylab("Averaged 90th percentile of log-likelihood") +
    
    scale_x_continuous(
      
      # Features of the first axis
      name = expression(xi),
      
      # Add a second axis and specify its features
      sec.axis = sec_axis(trans= ~.,
                          name= expression("b"["AR"]~("b"["PEA"])),
                          breaks=c(0, 5, 10),
                          labels=c("1.71 (1.71)", "1.17 (1.03)","1.09 (1.01)"))
      
    )
  
  print(gg)
  invisible(list(mean_data_glare, mean_data_big))
}

# Plot function for fixed xi --------------------------------------------------
plot_fixi <- function(sim_data) {
  
  gg_data <- sim_data
  gg_data$xi_val <- as.character(sim_data$xi) 
  
  gg1 <- ggplot(gg_data, aes(y = logLik_pert, x = v, color = xi_val, group = xi_val)) +
    
    geom_line() +
    
    labs(color = expression(xi)) +
    ylab("90th percentile of log-likelihood") +
    xlab("v")
  
  gg_data2 <- gg_data[gg_data$v <= 3 & gg_data$v >= -3, ]
  
  gg2 <- ggplot(gg_data2, aes(y = logLik_pert, x = v, color = xi_val, group = xi_val)) +
    
    geom_line() +
    
    labs(color = expression(xi)) +
    ylab("90th percentile of log-likelihood") +
    xlab("v")
  
  print(gg1)
  print(gg2)
  invisible(gg_data)
}

# Plot function for fixed xi for deviance res ---------------------------------
plot_deviance <- function(sim_data) {
  
  gg_data <- sim_data
  gg_data$xi_val <- as.character(sim_data$xi) 
  
  gg1 <- ggplot(gg_data, aes(y = deviance, x = v, color = xi_val, group = xi_val)) +
    
    geom_line() +
    
    labs(color = expression(xi)) +
    ylab("Root of averaged squared deviance residuals") +
    xlab("v")
  
  gg_data2 <- gg_data[gg_data$v <= 3 & gg_data$v >= -3, ]
  
  gg2 <- ggplot(gg_data2, aes(y = deviance, x = v, color = xi_val, group = xi_val)) +
    
    geom_line() +
    
    labs(color = expression(xi)) +
    ylab("Root of averaged squared deviance residuals") +
    xlab("v")
  
  print(gg1)
  print(gg2)
  invisible(gg_data)
}


# Plot function for fixed xi for pearson res ----------------------------------
plot_pearson <- function(sim_data) {
  
  gg_data <- sim_data
  gg_data$xi_val <- as.character(sim_data$xi) 
  
  gg1 <- ggplot(gg_data, aes(y = pearson, x = v, color = xi_val, group = xi_val)) +
    
    geom_line() +
    
    labs(color = expression(xi)) +
    ylab("Root of averaged squared pearson residuals") +
    xlab("v")
  
  gg_data2 <- gg_data[gg_data$v <= 3 & gg_data$v >= -3, ]
  
  gg2 <- ggplot(gg_data2, aes(y = pearson, x = v, color = xi_val, group = xi_val)) +
    
    geom_line() +
    
    labs(color = expression(xi)) +
    ylab("Root of averaged squared pearson residuals") +
    xlab("v")
  
  print(gg1)
  print(gg2)
  invisible(gg_data)
}


# Plot function for fixed v interventions for X -------------------------------
plot_fivi_X <- function(sim_data, xi_big = 10000) {
  
  axis_function <- function(value) {mean_data_glare$xi[which(abs(mean_data_glare$mean_b-value)==min(abs(mean_data_glare$mean_b-value)))]}
  
  sim_data_glare <- sim_data[sim_data$xi != xi_big, ]
  sim_data_big <- sim_data[sim_data$xi == xi_big, ]
  
  mean_data_glare <- sim_data_glare %>%
    group_by(xi) %>%
    summarise(mean_b = mean(b), mean_se = mean(b_se),
              mean_logLike_pert = mean(logLik_pert))
  
  mean_data_big <- sim_data_big %>%
    group_by(xi) %>%
    summarise(mean_b = mean(b), mean_se = mean(b_se),
              mean_logLike_pert = mean(logLik_pert))
  
  gg_data <- mean_data_glare
  
  gg <- ggplot(gg_data, aes(y = mean_logLike_pert, x = xi)) +
    
    geom_line() +
    geom_hline(yintercept = mean_data_big$mean_logLike_pert, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    
    ylab("Averaged 90th percentile of log-likelihood") +
    
    scale_x_continuous(
      
      # Features of the first axis
      name = expression(xi),
      # limits = c(-0.142,-0.020),
      # breaks = c(-0.14 ,-0.111, -0.083, -0.054 ,-0.025),
      
      # Add a second axis and specify its features
      sec.axis = sec_axis(trans= ~.,
                          name= "b",
                          breaks=c(0, 2.5, 5, 7.5, 10),
                          labels= round(
                            mean_data_glare$mean_b[mean_data_glare$xi %in% c(0, 2.5, 5, 7.5, 10)], digits = 3)
                            )
      
      #mean_data_glare$mean_b[mean_data_glare$xi == 0]
    )
  
  print(gg)
  invisible(list(mean_data_glare, mean_data_big))
  
  
  gg <- ggplot(gg_data, aes(y = mean_logLike_pert, x = mean_b)) +
    
    geom_line() +
    geom_hline(yintercept = mean_data_big$mean_logLike_pert, linetype = "dashed") +
    geom_vline(xintercept = gg_data$mean_b[gg_data$xi == 0], linetype = "dashed") +
    
    ylab("Averaged 90th percentile of log-likelihood") + 
    #xlab("b")
    
    scale_x_continuous(
      
      
      
      # Features of the first axis
      name = "b",
      
      # Add a second axis and specify its features
      sec.axis = sec_axis(trans= ~.,
                          name=expression(xi),
                          breaks=c(0.50, 0.55, 0.60, 0.65, 0.70, 0.75),
                          labels=sapply(c(0.50, 0.55, 0.60, 0.65, 0.70, 0.75), axis_function))
      
      #mean_data_glare$xi[which(abs(mean_data_glare$mean_b-0.54375)==min(abs(mean_data_glare$mean_b-0.54375)))]
    ) 
  
  print(gg)
  invisible(list(mean_data_glare, mean_data_big))
}

# Plot function for fixed v interventions for XHY -----------------------------
plot_fivi_XHY <- function(sim_data, xi_big = 10000) {
  
  sim_data_glare <- sim_data[sim_data$xi != xi_big, ]
  sim_data_big <- sim_data[sim_data$xi == xi_big, ]
  
  mean_data_glare <- sim_data_glare %>%
    group_by(xi) %>%
    summarise(mean_b = mean(b), mean_se = mean(b_se),
              mean_logLike_pert = mean(logLik_pert))
  
  mean_data_big <- sim_data_big %>%
    group_by(xi) %>%
    summarise(mean_b = mean(b), mean_se = mean(b_se),
              mean_logLike_pert = mean(logLik_pert))
  
  gg_data <- mean_data_glare
  
  gg <- ggplot(gg_data, aes(y = mean_logLike_pert, x = xi)) +
    
    geom_line() +
    geom_hline(yintercept = mean_data_big$mean_logLike_pert, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    
    ylab("Averaged 90th percentile of log-likelihood") +
    
    scale_x_continuous(
      
      # Features of the first axis
      name = expression(xi),
      # limits = c(-0.142,-0.020),
      # breaks = c(-0.14 ,-0.111, -0.083, -0.054 ,-0.025),
      
      # Add a second axis and specify its features
      sec.axis = sec_axis(trans= ~.,
                          name= "b",
                          breaks = c(0, 2.5, 5, 7.5, 10),
                          labels = round(
                            mean_data_glare$mean_b[mean_data_glare$xi %in% c(0, 2.5, 5, 7.5, 10)], digits = 3))
      
      #mean_data_glare$mean_b[mean_data_glare$xi == 0]
    )
  
  print(gg)
  invisible(list(mean_data_glare, mean_data_big))
  
  
  gg <- ggplot(gg_data, aes(y = mean_logLike_pert, x = mean_b)) +
    
    geom_line() +
    geom_hline(yintercept = mean_data_big$mean_logLike_pert, linetype = "dashed") +
    geom_vline(xintercept = gg_data$mean_b[gg_data$xi == 0], linetype = "dashed") +
    
    ylab("Averaged 90th percentile of log-likelihood") + 
    #xlab("b")
    
    scale_x_continuous(
      
      # Features of the first axis
      name = "b",
      # limits = c(-0.142,-0.020),
      breaks = c(-0.300, -0.225, -0.150, -0.075, 0.000),
      
      # Add a second axis and specify its features
      sec.axis = sec_axis(trans= ~.,
                          name=expression(xi),
                          breaks=c(-0.225, -0.150, -0.075),
                          labels=c(  "2.6", "0.7", "0.2"))
      
      #mean_data_glare$xi[which(abs(mean_data_glare$mean_b-0.54375)==min(abs(mean_data_glare$mean_b-0.54375)))]
    ) 
  
  print(gg)
  invisible(list(mean_data_glare, mean_data_big))
}



# Plot IV function ------------------------------------------------------------

plot_IV <- function(data, causal_parameter = 0.4, xi_big = 10000) {
  
  # Average over simulations
  mean_data <- data %>%
    group_by(xi, nobs) %>%
    summarise(mean_b = mean(b), se_b = sd(b))
  
  # Calculate CI
  upper <- mean_data$mean_b + qt(0.975, 100 - 1) * mean_data$se_b / sqrt(100)
  lower <- mean_data$mean_b - qt(0.975, 100 - 1) * mean_data$se_b / sqrt(100)
  mean_data$upper <- upper
  mean_data$lower <- lower
  
  # xi big without CI
  gg_data_big <- mean_data[mean_data$xi == xi_big, ]
  gg1 <- ggplot(gg_data_big, aes(y = mean_b, x = nobs, group = 1)) +
    
    geom_point(stat='summary', fun=sum) +
    stat_summary(fun=sum, geom="line") +
    
    geom_hline(yintercept = causal_parameter, linetype = "dashed") +
    
    xlab("Number of Observation") +
    ylab("Causal Parameter Estimate")
  
  
  # xi big with CI, nobs all
  gg2 <- ggplot(data = gg_data_big, aes(y = mean_b, x = nobs, group = 1)) +
    geom_point() +
    geom_line() + 
    geom_errorbar(aes(x = nobs, ymin = lower, ymax = upper)) +
    
    geom_hline(yintercept = causal_parameter, linetype = "dashed") +
    
    xlab("Number of Observation") +
    ylab("Causal Parameter Estimate")
  
  # xi big with CI, nobs 100-5000  
  gg_data_big2 <- gg_data_big[gg_data_big$nobs >= 100 & gg_data_big$nobs <= 5000, ]
  gg3 <- ggplot(data = gg_data_big2, aes(y = mean_b, x = nobs, group = 1)) +
    geom_point() +
    geom_line() + 
    geom_errorbar(aes(x = nobs, ymin = lower, ymax = upper)) +
    
    geom_hline(yintercept = causal_parameter, linetype = "dashed") +
    
    xlab("Number of Observation") +
    ylab("Causal Parameter Estimate")
  
  # xi all with CI, nobs all
  gg_data_all <- mean_data
  gg_data_all$xi <- as.factor(gg_data_all$xi)
  gg4 <- ggplot(data = gg_data_all, aes(y = mean_b, x = nobs, group = xi, color = xi)) +
    geom_point() +
    geom_line() + 
    geom_errorbar(aes(x = nobs, ymin = lower, ymax = upper)) +
    
    geom_hline(yintercept = causal_parameter, linetype = "dashed") +
    
    xlab("Number of Observation") +
    ylab("Causal Parameter Estimate") +
    labs(color = expression(xi))
  
  # xi all with CI, nobs 100-5000
  gg_data_all2 <- gg_data_all[gg_data_all$nobs >= 100 & gg_data_all$nobs <= 5000, ]
  gg5 <- ggplot(data = gg_data_all2, aes(y = mean_b, x = nobs, group = xi, color = xi)) +
    geom_point() +
    geom_line() + 
    geom_errorbar(aes(x = nobs, ymin = lower, ymax = upper)) +
    
    geom_hline(yintercept = causal_parameter, linetype = "dashed") +
    
    xlab("Number of Observation") +
    ylab("Causal Parameter Estimate") +
    labs(color = expression(xi))
  
  print(gg1)
  print(gg2)
  print(gg3)
  print(gg4)
  print(gg5)
  
  invisible(mean_data)
  
}

# -----------------------------------------------------------------------------
# Plot function for fixed v interventions with varying alpha quantile ---------

plot_fivi_quantiles <- function(sim_data, q_values, xi_big = 10000) {
  
  xi_len <- length(sim_data[[1]][[2]][, 1])
  q_len <- length(q_values)
  nsim <- length(sim_data)
  
  quantile_array <- array(dim = c(nsim, xi_len, q_len))
  for (r in 1:nsim) {
    for (q in 1:q_len) {
      quant <- q_values[q]
      for (x in 1:xi_len) {
        
        quantile_array[r, x, q] <- 
          as.numeric(quantile(sim_data[[r]][[2]][x, -1], quant))
        
      }
    }
  }
  
  # Plot for several xi
  averaged_quantiles <- apply(quantile_array, c(2,3), mean)
  colnames(averaged_quantiles) <- q_values
  
  xi_values <- as.numeric(sim_data[[1]][[2]][, 1])
  xi_values <- matrix(xi_values, ncol = 1)
  colnames(xi_values) <- "xi"
  
  quant_data_frame <- as.data.frame(cbind(xi_values, averaged_quantiles))
  gg_data <-  melt(quant_data_frame, id.vars = 'xi')
  
  gg_data2 <- gg_data[gg_data$xi != xi_big, ,]
  gg_data2$variable <- as.numeric(as.character(gg_data2$variable))
  #gg_data2$xi <- as.factor(gg_data2$xi)
  
  gg1 <- ggplot(gg_data2, aes(y = value, x = variable, color = xi, group = xi)) +
    geom_line() +
    ylab("Averaged quantiles of log-likelihood") +
    xlab(expression(alpha)) +
    labs(color = expression(xi))
  
  
  # Plot for several alpha
  averaged_quantiles <- apply(quantile_array, c(2,3), mean)
  rownames(averaged_quantiles) <- xi_values
  
  qqq <- matrix(q_values, nrow = 1)
  rownames(qqq) <- "alpha"
  
  quant_data_frame <- as.data.frame(t(rbind(qqq, averaged_quantiles)))
  gg_data <-  melt(quant_data_frame, id.vars = 'alpha')
  
  gg_data2 <- gg_data[gg_data$variable != xi_big, ,]
  gg_data2$variable <- as.numeric(as.character(gg_data2$variable))
  
  gg2 <- ggplot(gg_data2, aes(y = value, x = variable, color = alpha, group = alpha)) +
    geom_line() +
    ylab("Averaged quantiles of log-likelihood") +
    xlab(expression(xi)) +
    labs(color = expression(alpha))
  
  print(gg1)
  print(gg2)
  invisible(sim_data)
  
}














