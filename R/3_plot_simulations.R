library(tidyverse)
library(rsimsum)
library(ggplot2)

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
                                  group = c(rep("DEV", length(mean_data$mean_b_dev[!is.na(mean_data$mean_b_dev)])),
                                            rep("PEA", length(mean_data$mean_b_pea[!is.na(mean_data$mean_b_pea)]))))

  gg_parameter <- ggplot(gg_parameter_data, aes(y = y, x = x, color = group, group = group)) +
    
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
                        group = c(rep("AR", length(mean_data$mean_b_ar)),
                                  rep("DEV", length(mean_data$mean_b_dev[!is.na(mean_data$mean_b_dev)])),
                                  rep("PEA", length(mean_data$mean_b_pea[!is.na(mean_data$mean_b_pea)]))))
  
  # Use this data if only AR vs DEV wanted
  gg_data2 <- data.frame(x = c(mean_data$mean_b_ar,
                               mean_data$mean_b_dev[!is.na(mean_data$mean_b_dev)]),
                         y = c(mean_data$mean_MSE_pert_ar,
                               mean_data$mean_MSE_pert_dev[!is.na(mean_data$mean_b_dev)]),
                         group = c(rep("AR", length(mean_data$mean_b_ar)),
                                   rep("DEV", length(mean_data$mean_b_dev[!is.na(mean_data$mean_b_dev)]))))
  
  gg <- ggplot(gg_data, aes(y = y, x = x, color = group, group = group)) +
    
    geom_line() +
    
    ylab("MSE of perturbed data") +
    
    annotate("point", colour = 1, x = b_OLS_ar, y = MSE_OLS_ar) +
    annotate("text", label = "OLS", x = b_OLS_ar - 0.02, y = MSE_OLS_ar + 0.1) +
    annotate("point", colour = 1, x = b_IV_ar, y = MSE_IV_ar) +
    annotate("text", label = "IV", x = b_IV_ar + 0.02, y = MSE_IV_ar + 0.1) +
    annotate("point", colour = 1, x = b_PA_ar, y = MSE_PA_ar) +
    annotate("text", label = "PA", x = b_PA_ar - 0.02, y = MSE_PA_ar + 0.1) +
    
    scale_x_continuous(
      
      limits = c(0.95,2.05),
      
      # Features of the first axis
      name = "b",
      
      # Add a second axis and specify its features
      sec.axis = sec_axis(trans= ~.,
                          name="gamma (xi)",
                          breaks=c(1.111111, 1.555556, 2.000000),
                          labels=c("17.7 (8.3)", "1.8 (0.4)", "0 (-0.5)"))
      
    ) 
  
  print(gg_parameter)
  print(gg)
  invisible(mean_data)
}

# Plot function for fixed v interventions for Rothenhaeusler ------------------
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
                        group = c(rep("DEV", length(mean_data_glare$mean_b_dev)),
                                  rep("PEA", length(mean_data_glare$mean_b_pea))))

  gg <- ggplot(gg_data, aes(y = y, x = xi, color = group, group = group)) +
    
    geom_line() +
    geom_hline(yintercept = mean_data_big$mean_logLike_pert_dev, linetype = "dashed") +
    
    ylab("Averaged quantile of log-likelihood") +
    
    scale_x_continuous(
      
      # Features of the first axis
      name = "xi",
      
      # Add a second axis and specify its features
      sec.axis = sec_axis(trans= ~.,
                          name="b_dev (b_pea)",
                          breaks=c(0, 5, 10),
                          labels=c("1.71 (1.71)", "1.17 (1.03)","1.09 (1.01)"))
      
    )
  
  print(gg)
  invisible(list(mean_data_glare, mean_data_big))
}

# Plot function for fixed xi interventions for Rothenhaeusler -----------------
plot_rot_fixi <- function(sim_data) {

  gg_data <- sim_data
  gg_data$xi_val <- as.character(sim_data$xi) 
  
  gg <- ggplot(gg_data, aes(y = logLik_pert, x = v, color = xi_val, group = xi_val)) +
    
    geom_line() +
    
    ylab("Averaged quantile of log-likelihood") +
    xlab("|v|")

  print(gg)
  invisible(gg_data)
}









# Plot function for fixed v interventions -------------------------------------
plot_fivi <- function(sim_data, xi_big = 10000) {
  
  sim_data_glare <- sim_data[sim_data$xi != xi_big, ]
  sim_data_big <- sim_data[sim_data$xi == xi_big, ]
  
  mean_data_glare <- sim_data_glare %>%
    group_by(xi) %>%
    summarise(mean_b = mean(b), mean_se = mean(se),
              mean_logLike_pert = mean(logLik_pert))
  
  mean_data_big <- sim_data_big %>%
    group_by(xi) %>%
    summarise(mean_b = mean(b), mean_se = mean(se),
              mean_logLike_pert = mean(logLik_pert))
  
  # parameter x-axis
  plot(mean_data_glare$mean_b, mean_data_glare$mean_logLike_pert,
       ylim = c(min(mean_data_glare$mean_logLike_pert,
                    mean_data_big$mean_logLike_pert),
                max(mean_data_glare$mean_logLike_pert,
                    mean_data_big$mean_logLike_pert)),
       col ="1",
       xlab = "Averaged paramerter b",
       ylab = "Averaged quantile of log-likelihood")
  abline(v = mean_data_glare$mean_b[mean_data_glare$xi == 0], lwd=1, lty = 2, col = 2)
  abline(h = mean_data_big$mean_logLike_pert,
         lwd=1, lty = 3, col = 1)
  
  # xi x-axis
  plot(mean_data_glare$xi, mean_data_glare$mean_logLike_pert,
       ylim = c(min(mean_data_glare$mean_logLike_pert,
                    mean_data_big$mean_logLike_pert),
                max(mean_data_glare$mean_logLike_pert,
                    mean_data_big$mean_logLike_pert)),
       type = "l", col ="1",
       xlab = "xi", ylab = "Averaged quantile of log-likelihood")
  abline(v = 0, lwd=1, lty = 2, col = 2)
  abline(h = mean_data_big$mean_logLike_pert,
         lwd=1, lty = 3, col = 1)
  
  invisible(list(mean_data_glare, mean_data_big))
}

# Plot function for fixed xi interventions ------------------------------------
plot_fixi <- function(sim_data) {
  mean_data <- sim_data %>%
    group_by(v) %>%
    summarise(mean_b_glm = mean(b_glm), mean_b_se_glm = mean(se_glm),
              mean_logLike_glm_pert = mean(logLik_glm_pert),
              mean_b_glare = mean(b_glare), mean_se_glare = mean(se_glare),
              mean_logLike_glare_pert = mean(logLik_glare_pert),
              mean_b_big = mean(b_big), mean_se_big = mean(se_big),
              mean_logLike_big_pert = mean(logLik_big_pert))
  
  plot(mean_data$v, mean_data$mean_logLike_glm_pert,
       type = "l", col ="2", ylim = c(min(mean_data$mean_logLike_glm_pert,
                                          mean_data$mean_logLike_glare_pert,
                                          mean_data$mean_logLike_big_pert),
                                      max(mean_data$mean_logLike_glm_pert,
                                          mean_data$mean_logLike_glare_pert,
                                          mean_data$mean_logLike_big_pert)),
       xlab = "| v |", ylab = "Averaged quantile of log-likelihood")
  lines(mean_data$v, mean_data$mean_logLike_glare_pert, col = "3")
  lines(mean_data$v, mean_data$mean_logLike_big_pert, col = "4")
  
  invisible(mean_data)
}
