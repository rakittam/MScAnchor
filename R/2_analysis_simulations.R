library(tidyverse)
library(rsimsum)

# AR vs GLARE -----------------------------------------------------------------

sim_data_rot
head(sim_data_rot)

mean_data <- sim_data_rot %>%
  group_by(xi) %>%
  summarise(mean_b_ar = mean(b_ar),
            mean_b_glare = mean(b_glare),
            mean_logLike_pert = mean(loglikelihood_pert))
mean_data

plot(mean_data$mean_b_ar, mean_data$mean_b_glare)
plot(mean_data$xi, mean_data$mean_b_glare)
plot(mean_data$xi, mean_data$mean_logLike_pert)
index_temp <- which.max(mean_data$mean_logLike_pert)

mean_data[index_temp, ] # parameter and xi at optim
2 * mean_data$xi[index_temp] + 1 # gamma at optim
  
# Plot function for 1-dim interventions ---------------------------------------
plot_1dim <- function(sim_data) {
  mean_likeli_data <- sim_data %>%
    group_by(v) %>%
    summarise(mean_b_glm = mean(b_glm),
              mean_logLike_glm_pert = mean(loglikelihood_glm_pert),
              mean_b_glare = mean(b_glare),
              mean_logLike_glare_pert = mean(loglikelihood_glare_pert),
              mean_b_big = mean(b_big),
              mean_logLike_big_pert = mean(loglikelihood_big_pert))
  
  plot(mean_likeli_data$v, mean_likeli_data$mean_logLike_glm_pert,
       type = "l", col ="2", ylim = c(min(mean_likeli_data$mean_logLike_glm_pert,
                                          mean_likeli_data$mean_logLike_glare_pert),
                                      max(mean_likeli_data$mean_logLike_glm_pert,
                                          mean_likeli_data$mean_logLike_glare_pert)))
  lines(mean_likeli_data$v, mean_likeli_data$mean_logLike_glare_pert, col = "3")
  lines(mean_likeli_data$v, mean_likeli_data$mean_logLike_big_pert, col = "4")
  
  invisible(mean_likeli_data)
}

# Investigate data ------------------------------------------------------------

# Intervention on X
head(sim_data_bin_X)
summary(sim_data_bin_X)
plot_1dim(sim_data_bin_X)

# Intervention on H
plot_1dim(sim_data_bin_H)

# Intervention on Y
plot_1dim(sim_data_bin_Y)


# Plot function greater then 1-dim interventions ------------------------------
plot_multidim <- function(sim_data, sim_data_add) {
  mean_likeli_data <- sim_data %>%
    group_by(xi) %>%
    summarise(mean_b_glare = mean(b_glare),
              mean_logLike_glare_pert = mean(loglikelihood_glare_pert))
  
  mean_likeli_data_add <- sim_data_add %>%
    group_by() %>%
    summarise(mean_b_glm = mean(b_glm),
              mean_logLike_glm_pert = mean(loglikelihood_glm_pert),
              mean_b_big = mean(b_big),
              mean_logLike_big_pert = mean(loglikelihood_big_pert))
  
  plot(mean_likeli_data$xi, mean_likeli_data$mean_logLike_glare_pert,
       type = "l", col ="1")

  points(0, mean_likeli_data_add$mean_logLike_glm_pert, col = "2")
  abline(h = mean_likeli_data_add$mean_logLike_big_pert, col = "3")
  
  invisible(list(mean_likeli_data, mean_likeli_data_add))
}

# Investigate data ------------------------------------------------------------

# Intervention on X, H & Y
head(sim_data_bin_XHY)
head(sim_data_bin_XHY_add)

plot_multidim(sim_data_bin_XHY, sim_data_bin_XHY_add)
