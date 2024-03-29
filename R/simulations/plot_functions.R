# Define plot functions
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
library(tidyverse)
library(rsimsum)
library(ggplot2)
library(reshape2)
library(ggpubr) # for theme and ggarrange
theme_set(theme_bw())

# -----------------------------------------------------------------------------
### Plot Function Definitions

# Plot function for Anchor Regression Comparison ------------------------------
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
  gg_parameter_data <-
    data.frame(x = c(mean_data$mean_b_ar[!is.na(mean_data$mean_b_dev)],
                     mean_data$mean_b_ar[!is.na(mean_data$mean_b_pea)]),
               y = c(mean_data$mean_b_dev[!is.na(mean_data$mean_b_dev)],
                     mean_data$mean_b_pea[!is.na(mean_data$mean_b_pea)]),
               Method = c(rep("DEV",length(mean_data$mean_b_dev[!is.na(mean_data$mean_b_dev)])),
                          rep("PEA", length(mean_data$mean_b_pea[!is.na(mean_data$mean_b_pea)]))))
  
  gg_parameter <- ggplot(gg_parameter_data, aes(y = y, x = x, color = Method,
                                                group = Method)) +
    
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype ="dashed") +
    geom_vline(xintercept = mean_data$mean_b_ar[mean_data$xi == 0],
               linetype = "dashed", col = "gray") +
    geom_vline(xintercept = mean_data$mean_b_ar[which.max(mean_data$xi)],
               linetype = "dashed", col = "gray") +
    
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
  
  gg_data <-
    data.frame(x = c(mean_data$mean_b_ar,
                     mean_data$mean_b_dev[!is.na(mean_data$mean_b_dev)],
                     mean_data$mean_b_pea[!is.na(mean_data$mean_b_pea)]),
               y = c(mean_data$mean_MSE_pert_ar,
                     mean_data$mean_MSE_pert_dev[!is.na(mean_data$mean_b_dev)],
                     mean_data$mean_MSE_pert_pea[!is.na(mean_data$mean_b_pea)]),
               Method = c(rep("AR", length(mean_data$mean_b_ar)),
                          rep("DEV", length(mean_data$mean_b_dev[!is.na(mean_data$mean_b_dev)])),
                          rep("PEA", length(mean_data$mean_b_pea[!is.na(mean_data$mean_b_pea)]))))
  
  # Use this data if only AR vs DEV wanted
  gg_data2 <-
    data.frame(x = c(mean_data$mean_b_ar,
                     mean_data$mean_b_dev[!is.na(mean_data$mean_b_dev)]),
               y = c(mean_data$mean_MSE_pert_ar,
                     mean_data$mean_MSE_pert_dev[!is.na(mean_data$mean_b_dev)]),
               Method = c(rep("AR", length(mean_data$mean_b_ar)),
                          rep("GLARE", length(mean_data$mean_b_dev[!is.na(mean_data$mean_b_dev)]))))
  
  gg <- ggplot(gg_data2, aes(y = y, x = x, color = Method, group = Method)) +
    
    geom_line(aes(linetype = Method)) +
    
    ylab(expression(MSE~of~perturbed~data)) +
    labs(linetype = NULL, color = NULL) +
    
    annotate("point", colour = 1, x = b_OLS_ar, y = MSE_OLS_ar) +
    annotate("text", label = "OLS", x = b_OLS_ar - 0.06, y = MSE_OLS_ar + 0.1) +
    annotate("point", colour = 1, x = b_IV_ar, y = MSE_IV_ar) +
    annotate("text", label = "IV", x = b_IV_ar + 0.00, y = MSE_IV_ar + 0.13) +
    annotate("point", colour = 1, x = b_PA_ar, y = MSE_PA_ar) +
    annotate("text", label = "PA", x = b_PA_ar - 0.06, y = MSE_PA_ar + 0.0) +
    
    scale_x_continuous(
      
      # Features of the first axis
      name = expression(hat(b)),
      
      # Add a second axis and specify its features
      sec.axis = sec_axis(trans= ~.,
                          name=expression(gamma~"["*xi*"]"),
                          breaks=c(1.111111, 1.555556, 2.000000),
                          labels=c("17.7 [8.3]", "1.8 [0.4]", "0 [-0.5]"))
    ) + theme(legend.position = c(0.2, 0.8))
  
  print(gg_parameter)
  print(gg)
  invisible(mean_data)

  # Plots for LaTeX
  ggsave(filename = "ex1sim1.pdf", plot = gg, height = 4, width = 6)
}

# Plot function for fixed v Simulation for Rothenhaeusler Example ----------
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
  
  gg_data <-
    data.frame(x = c(mean_data_glare$mean_b_dev,
                     mean_data_glare$mean_b_pea),
               y = c(mean_data_glare$mean_logLike_pert_dev,
                     mean_data_glare$mean_logLike_pert_pea),
               xi = c(mean_data_glare$xi, mean_data_glare$xi),
               Method = c(rep("DEV", length(mean_data_glare$mean_b_dev)),
                          rep("PEA", length(mean_data_glare$mean_b_pea))))
  
  gg <- ggplot(gg_data, aes(y = y, x = xi, color = Method, group = Method)) +
    
    geom_line() +
    geom_hline(yintercept = mean_data_big$mean_logLike_pert_dev,
               linetype = "dashed") +
    
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

# Plot function for fixed v interventions for X -------------------------------
plot_fivi_X <- function(sim_data, xi_big = 10000) {
  
  axis_function <- function(value) {
    mean_data_glare$xi[which(abs(mean_data_glare$mean_b - value) == 
                               min(abs(mean_data_glare$mean_b - value)))]}
  
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
  
  gg1 <- ggplot(gg_data, aes(y = mean_logLike_pert, x = xi)) +
    
    geom_line() +
    geom_hline(yintercept = mean_data_big$mean_logLike_pert,
               linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    
    ylab("0.9-quantile of -logLik") 

  # Plot output
  print(gg1)
  invisible(list(mean_data_glare, mean_data_big))
  
  # Plots for LaTeX -----------------------------------------------------------
  # # b-axis
  # gg_data_b <- rbind(gg_data, mean_data_big)
  # 
  # gg2 <- ggplot(gg_data_b, aes(y = mean_loglike_pert, x = mean_b)) +
  #   
  #   geom_line() +
  #   geom_hline(yintercept = mean_data_big$mean_loglike_pert,
  #              linetype = "dashed") +
  #   geom_vline(xintercept = gg_data$mean_b[gg_data$xi == 0],
  #              linetype = "dashed") +
  #   
  #   ylab("0.9-quantile of -loglik")
  # #xlab("b")
  # 
  # scale_x_continuous(
  #   
  #   
  #   
  #   # features of the first axis
  #   name = expression(hat(b)),
  #   
  #   # add a second axis and specify its features
  #   sec.axis = sec_axis(trans= ~.,
  #                       name=expression(xi),
  #                       breaks=c(0.50, 0.55, 0.60, 0.65, 0.70, 0.75),
  #                       labels=sapply(c(0.50, 0.55, 0.60, 0.65, 0.70, 0.75),
  #                                     axis_function)))
  # 
  # gg2
  # 
  # bb <- 0.55
  # gg_data_b$xi[which(abs(gg_data_b$mean_b-bb)==min(abs(gg_data_b$mean_b-bb)))]
  # ggsave(filename = "ex4sim1_b.pdf", plot = gg2, height = 4, width = 6)
  # 
  # # b-axis
  # gg_data_b <- rbind(gg_data, mean_data_big)
  # 
  # gg2 <- ggplot(gg_data_b, aes(y = mean_logLike_pert, x = mean_b)) +
  #   
  #   geom_line() +
  #   geom_hline(yintercept = mean_data_big$mean_logLike_pert,
  #              linetype = "dashed") +
  #   geom_vline(xintercept = gg_data$mean_b[gg_data$xi == 0], 
  #              linetype = "dashed") +
  #   
  #   ylab("0.9-quantile of -logLik")
  # #xlab("b")
  # 
  # scale_x_continuous(
  #   
  #   
  #   
  #   # Features of the first axis
  #   name = expression(hat(b)),
  #   
  #   # Add a second axis and specify its features
  #   sec.axis = sec_axis(trans= ~.,
  #                       name=expression(xi),
  #                       breaks=c(0.50, 0.55, 0.60, 0.65, 0.70, 0.75),
  #                       labels=sapply(c(0.50, 0.55, 0.60, 0.65, 0.70, 0.75),
  #                                     axis_function)))
  # 
  # gg2
  # 
  # bb <- 0.55
  # gg_data_b$xi[which(abs(gg_data_b$mean_b - bb)==min(abs(gg_data_b$mean_b-bb)))]
  # ggsave(filename = "ex4sim1_b.pdf", plot = gg2, height = 4, width = 6)
  # 
  # 
  # 
  # 
  # # Example 2 Plots for LaTeX
  # ggsave(filename = "sim1.pdf", plot = gg1, height = 4, width = 6)
  # 
  # gg3 <- ggplot(gg_data[gg_data$xi<=50, ], aes(y = mean_logLike_pert, x = xi)) +
  #   
  #   geom_line() +
  #   geom_hline(yintercept = mean_data_big$mean_logLike_pert,
  #              linetype = "dashed") +
  #   geom_vline(xintercept = 0, linetype = "dashed") +
  #   
  #   ylab("0.9-quantile of -logLik") +
  #   
  #   scale_x_continuous(
  #     
  #     # Features of the first axis
  #     name = expression(xi),
  #     # limits = c(-0.142,-0.020),
  #     # breaks = c(-0.14 ,-0.111, -0.083, -0.054 ,-0.025),
  #     
  #     # Add a second axis and specify its features
  #     sec.axis = sec_axis(trans= ~.,
  #                         name= expression(hat(b)),
  #                         breaks=c(0, 10, 20, 30, 40, 50),
  #                         labels= round(
  #                           mean_data_glare$mean_b[mean_data_glare$xi %in% 
  #                                                    c(0, 10, 20, 30, 40, 50)],
  #                           digits = 3)
  #     )
  #     
  #     #mean_data_glare$mean_b[mean_data_glare$xi == 0]
  #   ) +
  #   #annotate("text", label="b = 0.4", colour = 1, x = 30, y = 8.5)+
  #   annotate(geom = 'text', x = 30, y = 8.5,
  #            label = "atop(b == 0.4,hat(b)[xi == 10000] %~~% 0.445)",
  #            parse = TRUE)
  # gg3
  # # Plots for LaTeX
  # ggsave(filename = "ex2sim1.pdf", plot = gg3, height = 4, width = 6)
  # 
  # 
  # 
  # # Example 3
  # gg4 <- ggplot(gg_data, aes(y = mean_logLike_pert, x = xi)) +
  #   
  #   geom_line() +
  #   geom_hline(yintercept = mean_data_big$mean_logLike_pert,
  #              linetype = "dashed") +
  #   geom_vline(xintercept = 0, linetype = "dashed") +
  #   
  #   ylab("0.9-quantile of -logLik") +
  #   
  #   scale_x_continuous(
  #     
  #     # Features of the first axis
  #     name = expression(xi),
  #     # limits = c(-0.142,-0.020),
  #     # breaks = c(-0.14 ,-0.111, -0.083, -0.054 ,-0.025),
  #     
  #     # Add a second axis and specify its features
  #     sec.axis = sec_axis(trans= ~.,
  #                         name= expression(hat(b)),
  #                         breaks=c(0, 10, 20, 30, 40, 50),
  #                         labels= round(
  #                           mean_data_glare$mean_b[mean_data_glare$xi %in%
  #                                                    c(0, 10, 20, 30, 40, 50)],
  #                           digits = 3))) +
  #   annotate(geom = 'text', x = 30, y = 1.9,
  #            label = "atop(b == 0.4,hat(b)[xi == 10000] %~~% 0.464)",
  #            parse = TRUE)
  # 
  # # Plots for LaTeX
  # ggsave(filename = "ex3sim1.pdf", plot = gg4, height = 4, width = 6)
  # 
  # # Example 4
  # gg4 <- ggplot(gg_data, aes(y = mean_logLike_pert, x = xi)) +
  #   
  #   geom_line() +
  #   geom_hline(yintercept = mean_data_big$mean_logLike_pert,
  #              linetype = "dashed") +
  #   geom_vline(xintercept = 0, linetype = "dashed") +
  #   
  #   ylab("0.9-quantile of -logLik") +
  #   
  #   scale_x_continuous(
  #     
  #     # Features of the first axis
  #     name = expression(xi),
  #     # limits = c(-0.142,-0.020),
  #     # breaks = c(-0.14 ,-0.111, -0.083, -0.054 ,-0.025),
  #     
  #     # Add a second axis and specify its features
  #     sec.axis = sec_axis(trans= ~.,
  #                         name= expression(hat(b)),
  #                         breaks=c(0, 10, 20, 30, 40, 50),
  #                         labels= round(
  #                           mean_data_glare$mean_b[mean_data_glare$xi %in%
  #                                                    c(0, 10, 20, 30, 40, 50)],
  #                           digits = 3))) +
  #   annotate(geom = 'text', x = 30, y = 5,
  #            label = "atop(b == 0.4,hat(b)[xi == 10000] %~~% 0.195)",
  #            parse = TRUE)
  # 
  # # Plots for LaTeX
  # ggsave(filename = "ex4sim1.pdf", plot = gg4, height = 4, width = 6)
  # 
  # # Example 5
  # gg5 <- ggplot(gg_data, aes(y = mean_logLike_pert, x = xi)) +
  #   
  #   geom_line() +
  #   geom_hline(yintercept = mean_data_big$mean_logLike_pert,
  #              linetype = "dashed") +
  #   geom_vline(xintercept = 0, linetype = "dashed") +
  #   
  #   ylab("0.9-quantile of -logLik") +
  #   
  #   scale_x_continuous(
  #     
  #     # Features of the first axis
  #     name = expression(xi),
  #     
  #     # Add a second axis and specify its features
  #     sec.axis = sec_axis(trans= ~.,
  #                         name= expression(hat(b)),
  #                         breaks=c(0, 2.5, 5, 7.5, 10),
  #                         labels= round(
  #                           mean_data_glare$mean_b[mean_data_glare$xi %in%
  #                                                    c(0, 2.5, 5, 7.5, 10)],
  #                           digits = 3))) +
  #   annotate("text", label="b = 0.5", colour = 1, x = 7.5, y = 1.2)
  # 
  # # Plots for LaTeX
  # ggsave(filename = "ex5sim1.pdf", plot = gg5, height = 4, width = 6)
  # 
  # # Example 6
  # gg5 <- ggplot(gg_data, aes(y = mean_logLike_pert, x = xi)) +
  #   
  #   geom_line() +
  #   geom_hline(yintercept = mean_data_big$mean_logLike_pert,
  #              linetype = "dashed") +
  #   geom_vline(xintercept = 0, linetype = "dashed") +
  #   
  #   ylab("0.9-quantile of -logLik") +
  #   
  #   scale_x_continuous(
  #     
  #     # Features of the first axis
  #     name = expression(xi),
  #     
  #     # Add a second axis and specify its features
  #     sec.axis = sec_axis(trans= ~.,
  #                         name= expression(hat(b)),
  #                         breaks=c(0, 2.5, 5, 7.5, 10),
  #                         labels= round(
  #                           mean_data_glare$mean_b[mean_data_glare$xi %in%
  #                                                    c(0, 2.5, 5, 7.5, 10)],
  #                           digits = 3))) +
  #   annotate("text", label="b = 0.5", colour = 1, x = 7.5, y = 0.76)
  # 
  # # Plots for LaTeX
  # ggsave(filename = "ex6sim1.pdf", plot = gg5, height = 4, width = 6)

}

# Plot function for fixed xi --------------------------------------------------
plot_fixi <- function(sim_data) {
  
  viri_opt <- "C" # set color option for viridis
  
  mean_data <- sim_data %>%
    group_by(xi, v) %>%
    summarise(mean_logLik_pert = mean(logLik_pert),
              mean_dev = mean(deviance),
              mean_pea = mean(pearson))
  
  gg_data <- mean_data
  gg_data$xi_val <- as.factor(mean_data$xi)
  
  # # Example 5 and 6:
  # gg_data$xi_val <- factor(mean_data$xi, levels = c("glm", "0", "1", "3", "5",
  #                                                  "8", "10", "20", "50", "10000"))

  # Main Plot: alpha-quantile of -logLik
  gg1 <- ggplot(gg_data, aes(y = mean_logLik_pert, x = v, color = xi_val, group = xi_val)) +
    geom_line() +
    scale_color_viridis_d(option = viri_opt, end = 0.95)+
   labs(color = expression(xi)) +
    ylab("0.9-quantile of -logLik") +
    xlab("")  #xlab("v") 
  
  # Zoomed
  gg_data2 <- gg_data[gg_data$v <= 3 & gg_data$v >= -1, ]
  gg2 <- ggplot(gg_data2, aes(y = mean_logLik_pert, x = v, color = xi_val, group = xi_val)) +
    geom_line() +
    scale_color_viridis_d(option = viri_opt, end = 0.95)+
     labs(color = expression(xi)) +
    ylab("") + # ylab("0.9-quantile of -logLik") +
    xlab("") # xlab("v") 

  # Deviance Plot
  gg3 <- ggplot(gg_data, aes(y = mean_dev, x = v, color = xi_val, group = xi_val)) +
    geom_line() +
    scale_color_viridis_d(option = viri_opt, end = 0.95)+
    labs(color = expression(xi)) +
    ylab(expression(RQSE[dev]^0.9)) +
    xlab("") #xlab("v")
  
  gg4 <- ggplot(gg_data2, aes(y = mean_dev, x = v, color = xi_val, group = xi_val)) +
    geom_line() +
    scale_color_viridis_d(option = viri_opt, end = 0.95)+
    labs(color = expression(xi)) +
    ylab("") + # ylab(expression(RQSE[dev]^0.9)) +
    xlab("")# xlab("v")
  
  # Pearson Plot
  gg5 <- ggplot(gg_data, aes(y = mean_pea, x = v, color = xi_val, group = xi_val)) +
    geom_line() +
    scale_color_viridis_d(option = viri_opt, end = 0.95)+
    labs(color = expression(xi)) +
    ylab(expression(RQSE[pea]^0.9)) +
    xlab("v")
  
  gg6 <- ggplot(gg_data2, aes(y = mean_pea, x = v, color = xi_val, group = xi_val)) +
    geom_line() +
    scale_color_viridis_d(option = viri_opt, end = 0.95)+
    labs(color = expression(xi)) +
    ylab("") + # ylab(expression(RQSE[pea]^0.9)) +
    xlab("v")
  
  # Plots for output
  print(gg1)
  # print(gg2)
  # print(gg3)
  # print(gg4)
  # print(gg5)
  # print(gg6)
  invisible(gg_data)
  
  # Plots for LaTeX -----------------------------------------------------------
  arr <- ggarrange(gg1, gg2,
                   gg3, gg4,
                   gg5, gg6,
                   labels = c("A", "B", "C","D", "E", "F"),
                   ncol = 2, nrow = 3, common.legend = TRUE, legend="top",
                   font.label = list(face = "plain"))
  
  ggsave(filename = "sim2.pdf", plot = arr, height = 7.5, width = 5)
 
  # # plot for ex1
  # gg_data_rot <- gg_data[gg_data$v <= 10 & gg_data$v >= 0, ]
  # 
  # gg_rot <- ggplot(gg_data_rot, aes(y = mean_logLik_pert, x = v, color = xi_val, group = xi_val)) +
  # 
  #   geom_line() +
  #   scale_color_viridis_d(option = viri_opt, end = 0.95)+
  # 
  #   labs(color = expression(xi)) +
  #   ylab("0.9-quantile of -logLik") +
  #   xlab("|v|") +
  #   theme(legend.position = c(0.15, 0.75))
  # 
  # gg_rot
  # 
  # ggsave(filename = "ex1sim3.pdf", plot = gg_rot, height = 4, width = 6)

  # # Plots for RMSE LaTeX
  # arr2 <- ggarrange(gg3, gg4,
  #                  gg5, gg6,
  #                  labels = c("A", "B", "C","D"),
  #                  ncol = 2, nrow = 2, common.legend = TRUE, legend="right",
  #                  font.label = list(face = "plain"))

  # # Plots for LaTeX Ex5
  # arr <- ggarrange(gg1, 
  #                  gg3,
  #                  gg5,
  #                  labels = c("A", "B", "C"),
  #                  ncol = 1, nrow = 3, common.legend = TRUE, legend="right",
  #                  font.label = list(face = "plain"))
  # 
  # ggsave(filename = "ex5sim2.pdf", plot = arr, height = 7, width = 5)
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
    xlab(expression(n[obs])) +
    ylab(expression(hat(beta)))
  
  
  # xi big with CI, nobs all
  gg2 <- ggplot(data = gg_data_big, aes(y = mean_b, x = nobs, group = 1)) +
    geom_point() +
    geom_line() + 
    geom_errorbar(aes(x = nobs, ymin = lower, ymax = upper)) +
    geom_hline(yintercept = causal_parameter, linetype = "dashed") +
    xlab(expression(n[obs])) +
    ylab(expression(hat(beta)))
  
  # xi big with CI, nobs 100-5000  
  gg_data_big2 <- gg_data_big[gg_data_big$nobs >= 100 & gg_data_big$nobs <= 5000, ]
  gg3 <- ggplot(data = gg_data_big2, aes(y = mean_b, x = nobs, group = 1)) +
    geom_point() +
    geom_line() + 
    geom_errorbar(aes(x = nobs, ymin = lower, ymax = upper)) +
    geom_hline(yintercept = causal_parameter, linetype = "dashed") +
    xlab(expression(n[obs])) +
    ylab(expression(hat(beta)))
  
  # xi all with CI, nobs all
  gg_data_all <- mean_data
  gg_data_all$xi <- as.factor(gg_data_all$xi)
  gg4 <- ggplot(data = gg_data_all, aes(y = mean_b, x = nobs, group = xi, color = xi)) +
    geom_point() +
    geom_line() + 
    geom_errorbar(aes(x = nobs, ymin = lower, ymax = upper)) +
    scale_color_viridis_d(option = "C", end = 0.95)+
    geom_hline(yintercept = causal_parameter, linetype = "dashed") +
    xlab(expression(n[obs])) +
    ylab(expression(hat(beta))) +
    labs(color = expression(xi))
  
  # xi all with CI, nobs 100-5000
  gg_data_all2 <- gg_data_all[gg_data_all$nobs >= 100 & gg_data_all$nobs <= 5000, ]
  gg5 <- ggplot(data = gg_data_all2, aes(y = mean_b, x = nobs, group = xi, color = xi)) +
    geom_point(size=0.5) +
    geom_errorbar(aes(x = nobs, ymin = lower, ymax = upper), width=0.1, size=0.1) +
    geom_hline(yintercept = causal_parameter, linetype = "dashed") +
    scale_color_viridis_d(option = "C", end = 0.95)+
    xlab(expression(n[obs])) +
    ylab(expression(hat(beta))) +
    labs(color = expression(xi))
  
  # Plot output
  # print(gg1)
  # print(gg2)
  # print(gg3)
  # print(gg4)
  print(gg5)
  
  invisible(mean_data)
  
  # Plots for LaTeX
  pdf("sim3.pdf", height=5, width=5)
  print(gg5)
  dev.off()
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
  
  gg1 <- ggplot(gg_data, aes(y = mean_logLike_pert, x = xi)) +
    
    geom_line() +
    geom_hline(yintercept = mean_data_big$mean_logLike_pert,
               linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    
    ylab("0.9-quantile of -logLik") +
    
    scale_x_continuous(
      
      # Features of the first axis
      name = expression(xi),
      # limits = c(-0.142,-0.020),
      # breaks = c(-0.14 ,-0.111, -0.083, -0.054 ,-0.025),
      
      # Add a second axis and specify its features
      sec.axis = sec_axis(trans= ~.,
                          name= expression(hat(b)),
                          breaks = c(0, 10, 20, 30, 40, 50),
                          labels = round(
                            mean_data_glare$mean_b[mean_data_glare$xi %in%
                                                     c(0, 10, 20, 30, 40, 50)],
                            digits = 3)))
  
  # Plot output
  print(gg1)
  invisible(list(mean_data_glare, mean_data_big))
  ggsave(filename = "sim1.pdf", plot = gg1, height = 4, width = 6)
  
  # LaTeX plots ---------------------------------------------------------------

  # gg_data_b <- rbind(gg_data, mean_data_big)
  # 
  # gg2 <- ggplot(gg_data_b, aes(y = mean_logLike_pert, x = mean_b)) +
  # 
  #   geom_line() +
  #   geom_hline(yintercept = mean_data_big$mean_logLike_pert,
  #              linetype = "dashed") +
  #   geom_vline(xintercept = gg_data$mean_b[gg_data$xi == 0],
  #              linetype = "dashed") +
  # 
  #   ylab("0.9-quantile of -logLik") +
  #   #xlab("b")
  # 
  #   scale_x_continuous(
  # 
  #     # Features of the first axis
  #     name = expression(hat(b)),
  #     # limits = c(-0.142,-0.020),
  #     breaks = c(0.2, 0.3, 0.4, 0.5),
  # 
  #     # Add a second axis and specify its features
  #     sec.axis = sec_axis(trans= ~.,
  #                         name=expression(xi),
  #                         breaks=c(0.2, 0.3, 0.4, 0.5, 0.562),
  #                         labels=c("50", "4.1", "1.3", "0.3", "0"))
  #   )
  # 
  # gg2
  # 
  # #bb <- 0.562
  # #gg_data_b$xi[which(abs(gg_data_b$mean_b-bb)==min(abs(gg_data_b$mean_b-bb)))]
  # ggsave(filename = "ex4sim1_b.pdf", plot = gg2, height = 4, width = 6)
  # 
  # 
  # gg_data_b <- rbind(gg_data, mean_data_big)
  # 
  # gg2 <- ggplot(gg_data_b, aes(y = mean_logLike_pert, x = mean_b)) +
  # 
  #   geom_line() +
  #   geom_hline(yintercept = mean_data_big$mean_logLike_pert,
  #              linetype = "dashed") +
  #   geom_vline(xintercept = gg_data$mean_b[gg_data$xi == 0],
  #              linetype = "dashed") +
  # 
  #   ylab("0.9-quantile of -logLik") +
  #   #xlab("b")
  # 
  #   scale_x_continuous(
  # 
  #     # Features of the first axis
  #     name = expression(hat(b)),
  #     # limits = c(-0.142,-0.020),
  #     breaks = c(0.2, 0.3, 0.4, 0.5),
  # 
  #     # Add a second axis and specify its features
  #     sec.axis = sec_axis(trans= ~.,
  #                         name=expression(xi),
  #                         breaks=c(0.2, 0.3, 0.4, 0.5, 0.562),
  #                         labels=c("50", "4.1", "1.3", "0.3", "0")))
  # gg2
  # 
  # #bb <- 0.562
  # #gg_data_b$xi[which(abs(gg_data_b$mean_b-bb)==min(abs(gg_data_b$mean_b-bb)))]
  # ggsave(filename = "ex4sim1_b.pdf", plot = gg2, height = 4, width = 6)
  # 
  # 
  # # Plots for LaTeX
  # ggsave(filename = "ex4sim1.pdf", plot = gg1, height = 4, width = 6)
  # 
  # gg4 <- ggplot(gg_data, aes(y = mean_logLike_pert, x = xi)) +
  #   
  #   geom_line() +
  #   geom_hline(yintercept = mean_data_big$mean_logLike_pert,
  #              linetype = "dashed") +
  #   geom_vline(xintercept = 0, linetype = "dashed") +
  #   
  #   ylab("0.9-quantile of -logLik") +
  #   
  #   scale_x_continuous(
  #     
  #     # Features of the first axis
  #     name = expression(xi),
  #     # limits = c(-0.142,-0.020),
  #     # breaks = c(-0.14 ,-0.111, -0.083, -0.054 ,-0.025),
  #     
  #     # Add a second axis and specify its features
  #     sec.axis = sec_axis(trans= ~.,
  #                         name= expression(hat(b)),
  #                         breaks = c(0, 10, 20, 30, 40, 50),
  #                         labels = round(
  #                           mean_data_glare$mean_b[mean_data_glare$xi %in%
  #                                                    c(0, 10, 20, 30, 40, 50)],
  #                           digits = 3))
  #     
  #     #mean_data_glare$mean_b[mean_data_glare$xi == 0]
  #   ) +
  #   annotate("text", label="b = 0.4", colour = 1, x = 30, y = 4.75)
  # 
  # # Plots for LaTeX
  # ggsave(filename = "ex4sim1.pdf", plot = gg4, height = 4, width = 6)
  
}

# Plot function for fixed v interventions with varying alpha quantile ---------
plot_quant <- function(sim_data, q_values, xi_big = 10000,
                       glm_quantile_vector = NULL) {
  
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
    ylab(expression(alpha-quantiles~of~-logLik)) +
    xlab(expression(alpha)) +
    labs(color = expression(xi)) +
    scale_color_viridis_c(option = "C", end = 0.95)
  
  # Zoom in
  gg_data3 <- gg_data2[gg_data2$variable >= 0.9, ]
  
  gg3 <- ggplot(gg_data3, aes(y = value, x = variable, color = xi, group = xi)) +
    geom_line() +
    ylab(expression(alpha-quantiles~of~-logLik)) +
    xlab(expression(alpha)) +
    labs(color = expression(xi)) +
    scale_color_viridis_c(option = "C", end = 0.95)
  
  # Plot for several alpha
  averaged_quantiles <- apply(quantile_array, c(2,3), mean)
  rownames(averaged_quantiles) <- xi_values
  
  qqq <- matrix(q_values, nrow = 1)
  rownames(qqq) <- "alpha"
  
  quant_data_frame <- as.data.frame(t(rbind(qqq, averaged_quantiles)))
  gg_data3 <-  melt(quant_data_frame, id.vars = 'alpha')
  
  gg_data4 <- gg_data3[gg_data3$variable != xi_big, ,]
  gg_data4$variable <- as.numeric(as.character(gg_data4$variable))
  
  gg2 <- ggplot(gg_data4, aes(y = value, x = variable, color = alpha,
                              group = alpha)) +
    geom_line() +
    ylab(expression(alpha-quantiles~of~-logLik)) +
    xlab(expression(xi)) +
    labs(color = expression(alpha)) +
    scale_color_viridis_c(option = "C", end = 0.95)

  # Zoom in
  gg_data5 <- gg_data4[gg_data4$alpha >= 0.8, ]
  
  gg4 <- ggplot(gg_data5, aes(y = value, x = variable, color = alpha,
                              group = alpha)) +
    geom_line() +
    ylab(expression(alpha-quantiles~of~-logLik)) +
    xlab(expression(xi)) +
    labs(color = expression(alpha)) +
    scale_color_viridis_c(option = "C", end = 0.95)
  
  # Plot output
  print(gg2)
  print(gg1)
  # print(gg3)
  # print(gg4)
  invisible(sim_data)
  
  ggsave(filename = "sim4.pdf", plot = gg1, height = 4, width = 6)
  
  # # Plots for LaTeX ---------------------------------------------------------
  # arr0 <- ggarrange(gg1 + rremove("legend"), gg2 + rremove("legend"),
  #                   gg3 + rremove("legend"), gg4 + rremove("legend"), 
  #                   labels = c("A", "B", "C", "D"),
  #                   ncol = 2, nrow = 2)
  # 
  # arr1 <- ggarrange(gg1, gg2,
  #                   labels = c("A", "B"),
  #                   ncol = 1, nrow = 2)
  # 
  # arr2 <- ggarrange(gg3, gg4,
  #                   labels = c("A", "B"),
  #                   ncol = 1, nrow = 2)
  # 
  # g1 <- ggplotGrob(gg3 +
  #                    theme(legend.position = "none") +
  #                    labs(x=element_blank(), y=element_blank()))
  # alpax <- gg1+
  #   annotation_custom(
  #     grob = g1,
  #     xmin = -0.05,
  #     xmax = 0.9,
  #     ymin = 5,
  #     ymax = 15
  #   ) 
  # alpax
  # 
  # arr <- ggarrange(alpax, gg2,
  #                  labels = c("A", "B"),
  #                  ncol = 1, nrow = 2,
  #                  font.label = list(face = "plain"))
  # 
  # # Plots for LaTeX
  # ggsave(filename = "ex4sim4.pdf", plot = arr, height = 4, width = 6)
  # 
  # ggsave(filename = "ex4sim4.1.pdf", plot = alpax, height = 4, width = 6)
  # ggsave(filename = "ex4sim4.2.pdf", plot = gg2, height = 4, width = 6)
  # 
  # print(alpax)
  # print(gg2)
  # 
  # # Example 2
  # g1 <- ggplotGrob(gg3 +
  #                    theme(legend.position = "none") +
  #                    labs(x=element_blank(), y=element_blank()))
  # alpax <- gg1+
  #   annotation_custom(
  #     grob = g1,
  #     xmin = -0.05,
  #     xmax = 0.88,
  #     ymin = 15,
  #     ymax = 47
  #   ) 
  # alpax
  # 
  # arr <- ggarrange(alpax, gg2,
  #                  labels = c("A", "B"),
  #                  ncol = 1, nrow = 2,
  #                  font.label = list(face = "plain"))
  # 
  # # Plots for LaTeX
  # ggsave(filename = "ex2sim4.pdf", plot = arr, height = 4, width = 6)
  # 
  # ggsave(filename = "ex2sim4.1.pdf", plot = alpax, height = 4, width = 6)
  # ggsave(filename = "ex2sim4.2.pdf", plot = gg2, height = 4, width = 6)
  # 
  # print(alpax)
  # print(gg2)
  # 
# 
#   # Example 3
#   g1 <- ggplotGrob(gg3 +
#                      theme(legend.position = "none") +
#                      labs(x=element_blank(), y=element_blank()))
#   alpax <- gg1+
#     annotation_custom(
#       grob = g1,
#       xmin = -0.05,
#       xmax = 0.8,
#       ymin = 1.8,
#       ymax = 3.7
#     )
#   alpax
# 
#   arr <- ggarrange(alpax, gg2,
#                    labels = c("A", "B"),
#                    ncol = 1, nrow = 2,
#                    font.label = list(face = "plain"))
# 
#   # Plots for LaTeX
#   ggsave(filename = "ex3sim4.pdf", plot = arr, height = 4, width = 6)
# 
#   ggsave(filename = "ex3sim4.1.pdf", plot = alpax, height = 4, width = 6)
#   ggsave(filename = "ex3sim4.2.pdf", plot = gg2, height = 4, width = 6)
# 
#   print(alpax)
#   print(gg2)

# 
  # # Example 4
  # g1 <- ggplotGrob(gg3 +
  #                    theme(legend.position = "none") +
  #                    labs(x=element_blank(), y=element_blank()))
  # alpax <- gg1+
  #   annotation_custom(
  #     grob = g1,
  #     xmin = -0.05,
  #     xmax = 0.9,
  #     ymin = 7.5,
  #     ymax = 16.5
  #   )
  # alpax
  # 
  # arr <- ggarrange(alpax, gg2,
  #                  labels = c("A", "B"),
  #                  ncol = 1, nrow = 2,
  #                  font.label = list(face = "plain"))
  # 
  # # Plots for LaTeX
  # ggsave(filename = "ex4sim4.pdf", plot = arr, height = 4, width = 6)
  # 
  # ggsave(filename = "ex4sim4.1.pdf", plot = alpax, height = 4, width = 6)
  # ggsave(filename = "ex4sim4.2.pdf", plot = gg2, height = 4, width = 6)
  # 
  # print(alpax)
  # print(gg2)

  # # Example 5
  # g1 <- ggplotGrob(gg3 +
  #                    theme(legend.position = "none") +
  #                    labs(x=element_blank(), y=element_blank()))
  # alpax <- gg1+
  #   annotation_custom(
  # grob = g1,
  # xmin = 0.45,
  # xmax = 1.05,
  # ymin = -1.82,
  # ymax = -0.5
  #   )
  # alpax
  # 
  # arr <- ggarrange(alpax, gg2,
  #                  labels = c("A", "B"),
  #                  ncol = 1, nrow = 2,
  #                  font.label = list(face = "plain"))
  # 
  # # Plots for LaTeX
  # ggsave(filename = "ex5sim4.pdf", plot = arr, height = 4, width = 6)
  # 
  # ggsave(filename = "ex5sim4.1.pdf", plot = alpax, height = 4, width = 6)
  # ggsave(filename = "ex5sim4.2.pdf", plot = gg2, height = 4, width = 6)
  # 
  # print(alpax)
  # print(gg2)

  # for example 6 LN involved
  if(!is.null(glm_quantile_vector)) {
    
    glm_frame <- data.frame(xi = -1, variable = q_values, value = glm_quantile_vector)
    gg_glm_data <- rbind(gg_data2, glm_frame)
    
    gg_glm <- ggplot(gg_glm_data, aes(y = value, x = variable, color = xi, group = xi)) +
      geom_line() +
      ylab(expression(alpha-quantiles~of~-logLik)) +
      xlab(expression(alpha)) +
      labs(color = expression(xi)) +
      scale_color_viridis_c(option = "C", end = 0.95)
    gg_glm
    
    ggsave(filename = "ex6sim4_glm.pdf", plot = gg_glm, height = 4, width = 6)
  }
  
}














