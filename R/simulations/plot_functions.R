# Define plot functions
#
# Date: 14.12.20     Author: Maic Rakitta
###############################################################################

library(extrafont)
library(tidyverse)
library(rsimsum)
library(ggplot2)
library(gridExtra)
library(reshape2)

# # execute once to add fonts:
# font_import(pattern = "lmroman*")
# font_import(pattern = "lmroman10-regular-webfont.ttf") 
# 
# loadfonts()
# loadfonts(device = "postscript")
# 
# loadfonts(device = "win")
# 
# font_import()
 
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
    
    #theme(text = element_text(size=10, family="Latin Modern Roman 10")) +
    
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
  #grid.arrange(gg, gg_parameter, nrow=2)
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
    ylab("Averaged 90th percentile of log-likelihood") +
    xlab("v")
  
  gg2 <- ggplot(gg_data, aes(y = logLik_pert, x = v, color = xi_val, group = xi_val)) +
    
    geom_line() +
    
    ylab("Averaged 90th percentile of log-likelihood") +
    xlab("v") + scale_y_continuous(trans = scales::pseudo_log_trans(base = 10))

  print(gg1)
  invisible(gg_data)
}

# Plot function for fixed v interventions for X -------------------------------
plot_fivi_X <- function(sim_data, xi_big = 10000) {
  
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
                          labels=c("0.574", "0.490", "0.469", "0.459", "0.454"))
      
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
                          breaks=c(0.4538262, 0.48125, 0.5125, 0.54375, 0.5735984),
                          labels=c("10", "3.3", "1.3", "0.5", "0"))
      
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
    summarise(mean_b = mean(b), mean_se = mean(se),
              mean_logLike_pert = mean(logLik_pert))
  
  mean_data_big <- sim_data_big %>%
    group_by(xi) %>%
    summarise(mean_b = mean(b), mean_se = mean(se),
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
                          labels=c("-0.025", "-0.122", "-0.134", "-0.139", "-0.142"))
      
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
      breaks =c(-0.14 ,-0.111, -0.083, -0.054 ,-0.025),
      
      # Add a second axis and specify its features
      sec.axis = sec_axis(trans= ~.,
                          name=expression(xi),
                          breaks=c(-0.14 ,-0.111, -0.083, -0.054 ,-0.025),
                          labels=c("8.2", "1.6", "0.6", "0.2", "0"))
      
      #mean_data_glare$xi[which(abs(mean_data_glare$mean_b-0.54375)==min(abs(mean_data_glare$mean_b-0.54375)))]
    ) 
  
  print(gg)
  invisible(list(mean_data_glare, mean_data_big))
}



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