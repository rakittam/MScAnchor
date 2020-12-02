library(tidyverse)
library(rsimsum)

# AR vs GLARE -----------------------------------------------------------------

sim_data_rot
head(sim_data_rot)

mean_data <- sim_data_rot %>%
  group_by(xi) %>%
  summarise(mean_b_ar = mean(b_ar),
            mean_b_glare = mean(b_glare),
            mean_logLike_pert = mean(logLik_pert))
mean_data

plot(mean_data$mean_b_ar, mean_data$mean_b_glare)
plot(log(mean_data$mean_b_ar), log(mean_data$mean_b_glare))
plot(log(log(mean_data$mean_b_ar)), log(log(mean_data$mean_b_glare)))

plot(mean_data$xi, mean_data$mean_b_glare)
plot(mean_data$xi, mean_data$mean_logLike_pert)
index_temp <- which.max(mean_data$mean_logLike_pert)

mean_data[index_temp, ] # parameter and xi at optim
2 * mean_data$xi[index_temp] + 1 # gamma at optim
  
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

# Investigate data ------------------------------------------------------------

# Intervention on X
head(sim_data_bin_X)
summary(sim_data_bin_X)
plot_fixi(sim_data_bin_X)

# Intervention on H
plot_fixi(sim_data_bin_H)

# Intervention on Y
plot_fixi(sim_data_bin_Y)

# Investigate data ------------------------------------------------------------

# Intervention on X, H & Y
head(sim_data_bin_XHY)
head(sim_data_bin_XHY_add)

plot_fivi(sim_data_bin_XHY, sim_data_bin_XHY_add)
