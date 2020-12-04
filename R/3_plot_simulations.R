library(tidyverse)
library(rsimsum)
library(ggplot2)

# Plot function for anchor regression comparison ------------------------------
plot_rot <- function(sim_data, xi_big = 10000) {
  
  mean_data <- sim_data %>%
    group_by(xi, gamma) %>%
    summarise(mean_b_ar = mean(b_ar), mean_se_ar = mean(se_ar),
              mean_MSE_pert_ar = mean(MSE_pert_ar),
              mean_b_glare = mean(b_glare), mean_se_glare = mean(se_glare),
              mean_MSE_pert_glare = mean(MSE_pert_glare))
  
  # parameter x-axis
  plot(mean_data$mean_b_ar, mean_data$mean_b_glare)
  plot(log(mean_data$mean_b_ar), log(mean_data$mean_b_glare))
  
  
  # Plot ex 2 from Rothenhaeusler
  b_OLS_ar <- mean_data$mean_b_ar[mean_data$gamma == 1]
  b_IV_ar <- mean_data$mean_b_ar[mean_data$xi == xi_big]
  b_PA_ar <- mean_data$mean_b_ar[mean_data$gamma == 0]
  
  MSE_OLS_ar <- mean_data$mean_MSE_pert_ar[mean_data$gamma == 1]
  MSE_IV_ar <- mean_data$mean_MSE_pert_ar[mean_data$xi == xi_big]
  MSE_PA_ar <- mean_data$mean_MSE_pert_ar[mean_data$gamma == 0]
  
  
  gg_data <- data.frame(x = c(mean_data$mean_b_ar, mean_data$mean_b_glare),
                        y = c(mean_data$mean_MSE_pert_ar, mean_data$mean_MSE_pert_glare),
                        group = c(rep("AR", length(mean_data$mean_b_ar)),
                                      rep("GLARE", length(mean_data$mean_b_ar))))
  
  
  
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
      
      limits = c(0.99,2.01),
      
      # Features of the first axis
      name = "b",
      
      # Add a second axis and specify its features
      sec.axis = sec_axis(trans= ~.,
                          name="gamma (xi)",
                          breaks=c(1.111111, 1.555556, 2.000000),
                          labels=c("17.7 (8.3)", "1.8 (0.4)", "0 (-0.5)"))
      
    ) 
  
  # 
  # gg_ar <- ggplot(mean_data, aes(y = mean_MSE_pert_ar)) +
  #   
  #   geom_line(aes(x = mean_b_ar)) +
  #   
  #   annotate("point", colour = 2, x = b_OLS_ar, y = MSE_OLS_ar) +
  #   annotate("text", label = "OLS", x = b_OLS_ar - 0.02, y = MSE_OLS_ar + 0.1) +
  #   annotate("point", colour = 3, x = b_IV_ar, y = MSE_IV_ar) +
  #   annotate("text", label = "IV", x = b_IV_ar + 0.02, y = MSE_IV_ar + 0.1) +
  #   annotate("point", colour = 4, x = b_PA_ar, y = MSE_PA_ar) +
  #   annotate("text", label = "PA", x = b_PA_ar - 0.02, y = MSE_PA_ar + 0.1) +
  #   
  #   scale_x_continuous(
  #     
  #     limits = c(0.99,2.01),
  #     
  #     # Features of the first axis
  #     name = "b",
  #     
  #     # Add a second axis and specify its features
  #     sec.axis = sec_axis(trans= ~.,
  #                         name="gamma",
  #                         breaks=c(1.000000, 1.111111, 1.222222, 1.333333, 1.444444, 1.555556, 1.666667, 1.777778, 1.888889, 2.000000),
  #                         labels=c("inf",17.7,7.9, 4.5, 2.9, 1.8, 1.1, 0.7, 0.3, 0.0))
  #     
  #   )
  # 
  # # Plot ex2 analogon for glare
  # b_OLS_glare <- mean_data$mean_b_glare[mean_data$gamma == 1]
  # b_IV_glare <- mean_data$mean_b_glare[mean_data$xi == 10000]
  # b_PA_glare <- mean_data$mean_b_glare[mean_data$gamma == 0]
  # 
  # MSE_OLS_glare <- mean_data$mean_MSE_pert_glare[mean_data$gamma == 1]
  # MSE_IV_glare <- mean_data$mean_MSE_pert_glare[mean_data$xi == 10000]
  # MSE_PA_glare <- mean_data$mean_MSE_pert_glare[mean_data$gamma == 0]
  # 
  # gg_glare <- ggplot(mean_data, aes(y = mean_MSE_pert_glare)) +
  #   
  #   geom_line(aes(x = mean_b_glare)) +
  #   
  #   annotate("point", colour = 2, x = b_OLS_glare, y = MSE_OLS_glare) +
  #   annotate("text", label = "OLS", x = b_OLS_glare - 0.02, y = MSE_OLS_glare + 0.1) +
  #   annotate("point", colour = 3, x = b_IV_glare, y = MSE_IV_glare) +
  #   annotate("text", label = "IV", x = b_IV_glare + 0.02, y = MSE_IV_glare + 0.1) +
  #   annotate("point", colour = 4, x = b_PA_glare, y = MSE_PA_glare) +
  #   annotate("text", label = "PA", x = b_PA_glare - 0.02, y = MSE_PA_glare + 0.1) +
  #   
  #   scale_x_continuous(
  #     
  #     limits = c(0.99,2.01),
  #     
  #     # Features of the first axis
  #     name = "b",
  #     
  #     # Add a second axis and specify its features
  #     sec.axis = sec_axis(trans= ~.,
  #                         name="xi",
  #                         breaks=c(1.000000, 1.111111, 1.222222, 1.333333, 1.444444, 1.555556, 1.666667, 1.777778, 1.888889, 2.000000),
  #                         labels=c("inf", 8.3, 3.4, 1.8, 0.9, 0.4, 0.1, -0.2, -0.4, -0.5))
  #     
  #   )
  # 
  # print(gg_ar)
  # print(gg_glare)
  
  print(gg)
  invisible(mean_data)
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
