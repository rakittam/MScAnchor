library(tidyverse)
library(rsimsum)

sim_data_rot
head(sim_data_rot)

mean_data <- sim_data_rot %>%
  group_by(xi) %>%
  summarise(mean_b_ar = mean(b_ar),
            mean_b_glare = mean(b_glare),
            mean_logLike_pert = mean(loglikelihood_pert))
mean_data

plot(mean_data$mean_b_ar, mean_data$mean_b_glare)
plot(mean_data$mean_b_ar, mean_data$mean_b_glare)
plot(mean_data$xi, mean_data$mean_b_glare)
plot(mean_data$xi, mean_data$mean_logLike_pert)
which.max(mean_data$mean_logLike_pert)

mean_data[17, ] # parameter and xi at optim
2 * mean_data$xi[17] + 1 # gamma at optim
  





head(sim_data_bin_IV)
summary(sim_data_bin_IV)

mean_likeli_data <- sim_data_bin_IV %>%
  group_by(v) %>%
  summarise(mean_b_glm = mean(b_glm),
            mean_b_glare = mean(b_glare),
            mean_logLike_glm_pert = mean(loglikelihood_glm_pert),
            mean_logLike_glare_pert = mean(loglikelihood_glare_pert))
mean_likeli_data

plot(mean_likeli_data$v, mean_likeli_data$mean_logLike_glm_pert,
     type = "l", col ="2", ylim = c(min(mean_likeli_data$mean_logLike_glm_pert,
                                        mean_likeli_data$mean_logLike_glare_pert),
                                    max(mean_likeli_data$mean_logLike_glm_pert,
                                        mean_likeli_data$mean_logLike_glare_pert)))
lines(mean_likeli_data$v, mean_likeli_data$mean_logLike_glare_pert, col = "3")




