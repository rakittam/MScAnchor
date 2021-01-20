# Run simulation analysis and plots
#
# Date: 19.01.21     Author: Maic Rakitta
###############################################################################

# -----------------------------------------------------------------------------
### Example Data Analysis

# ex1: Rothenhaeusler example -------------------------------------------------

# sim1: Rothenhaeusler Comparison ---

load(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/simulation_study/ex1/2021-01-15/sim1.Rdata")

data_rot_states <- sim_data_rot$states
data_rot <- sim_data_rot$sim_data

head(data_rot)
summary(data_rot)

plot_rot(data_rot)


# sim2: Rothenhaeusler with fixed v ---

load(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/simulation_study/ex1/2021-01-15/sim2.Rdata")

data_rot_X_fivi_states <- sim_data_rot_X_fivi$states
data_rot_X_fivi <- sim_data_rot_X_fivi$sim_data

head(data_rot_X_fivi)
summary(data_rot_X_fivi)

plot_rot_fivi(data_rot_X_fivi)

# sim3: Rothenhaeusler with fixed xi ---

load(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/simulation_study/ex1/2021-01-15/sim3.Rdata")

# Compute data
data_rot_X_fixi <- data.frame(matrix(nrow = 0, ncol = 8))
colnames(data_rot_X_fixi) <- c("rep", "v", "xi",
                               "b", "se", "logLik_pert", "deviance", "pearson")
for (i in 1:nsim) {
  
  data_temp <- sim_data_rot_X_fixi[[i]][["data"]][["data"]]
  
  data_rot_X_fixi <- rbind(data_rot_X_fixi, data_temp)
}

head(data_rot_X_fixi)
summary(data_rot_X_fixi)

# Plot Quantile of log-likelihood and with residuals
plot_fixi(data_rot_X_fixi)

# ex2: Poisson with IV example ------------------------------------------------

# sim1: Poisson IV with fixed v ---

load(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/simulation_study/ex2/2021-01-18/sim1.Rdata")

data_poi_X_states_fivi <- sim_data_poi_X_fivi$states
data_poi_X_fivi <- sim_data_poi_X_fivi$sim_data

head(data_poi_X_fivi)
summary(data_poi_X_fivi)

plot_fivi_X(data_poi_X_fivi)

# sim2: Poisson IV with fixed xi ---

load(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/simulation_study/ex2/2021-01-18/sim2.Rdata")


# Compute data
data_poi_X_fixi <- data.frame(matrix(nrow = 0, ncol = 8))
colnames(data_poi_X_fixi) <- c("rep", "v", "xi",
                               "b", "se", "logLik_pert", "deviance", "pearson")
for (i in 1:nsim) {
  
  data_temp <- sim_data_poi_X_fixi[[i]][["data"]][["data"]]
  
  data_poi_X_fixi <- rbind(data_poi_X_fixi, data_temp)
}

head(data_poi_X_fixi)
summary(data_poi_X_fixi)

# Plot Quantile of log-likelihood and with residuals
plot_fixi(data_poi_X_fixi)

# Investigate different quantiles

alpha <- 0.9

quantile_data <- data.frame(matrix(nrow = 0, ncol = 6))
colnames(quantile_data) <- c("rep", "v", "xi", "logLik_pert", "deviance", "pearson")

for (i in 1:nsim) {
  
  data_temp <- sim_data_poi_X_fixi[[i]][["data"]][["indv"]]
  
  vlen <- length(data_temp)
  
  for (s in 1:vlen) {
    
    # Log-Likelihood
    logLik <- apply(data_temp[[s]][[1]], 1, function(x) quantile(x, probs = alpha))

    # Deviance Residuals
    deviance <- apply(data_temp[[s]][[2]], 1, function(x) quantile(x, probs = alpha))
    
    # Pearson Residuals
    pearson <- apply(data_temp[[s]][[3]], 1, function(x) quantile(x, probs = alpha))

    temp <- data.frame(rep = i, v = data_temp[[s]][[5]], xi = data_temp[[s]][[4]],
                       logLik_pert = logLik, deviance = deviance, pearson = pearson)
    
    quantile_data <- rbind(quantile_data, temp)
  }
  
  data_temp$logLik
  
  data_poi_X_fixi <- rbind(data_poi_X_fixi, data_temp)
}

plot_fixi(quantile_data)

# sim3: Poisson IV identifiability ---

load(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/simulation_study/ex2/2021-01-18/sim3.Rdata")

IV_b_data <- IV_b_poi$sim_data
plot_IV(IV_b_data, causal_parameter = 0.4)


# ex3: Binomial with IV example -----------------------------------------------

# sim1: Binomial IV with fixed v ---

load(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/simulation_study/ex3/2021-01-18/sim1.Rdata")

data_bin_X_states_fivi <- sim_data_bin_X_fivi$states
data_bin_X_fivi <- sim_data_bin_X_fivi$sim_data

head(data_bin_X_fivi)
summary(data_bin_X_fivi)

plot_fivi_X(data_bin_X_fivi)

# sim2: Binomial IV with fixed xi ---

load(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/simulation_study/ex3/2021-01-18/sim2.Rdata")

# Compute data
data_bin_X_fixi <- data.frame(matrix(nrow = 0, ncol = 8))
colnames(data_bin_X_fixi) <- c("rep", "v", "xi",
                               "b", "se", "logLik_pert", "deviance", "pearson")
for (i in 1:nsim) {
  
  data_temp <- sim_data_bin_X_fixi[[i]][["data"]][["data"]]
  
  data_bin_X_fixi <- rbind(data_bin_X_fixi, data_temp)
}

head(data_bin_X_fixi)
summary(data_bin_X_fixi)

# Plot Quantile of log-likelihood and with residuals
plot_fixi(data_bin_X_fixi)

# sim3: Binomial IV identifiability ---

load(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/simulation_study/ex3/2021-01-18/sim3.Rdata")

IV_b_data <- IV_b_bin$sim_data
plot_IV(IV_b_data, causal_parameter = 1)

# ex4: Poisson with Anchor on X, H and Y --------------------------------------

# sim1: Poisson with Anchor on X, H and Y fixed v ---

load(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/simulation_study/ex4/2021-01-15/sim1.Rdata")

data_poi_XHY_fivi_states <- sim_data_poi_XHY_fivi$states
data_poi_XHY_fivi <- sim_data_poi_XHY_fivi$sim_data

head(data_poi_XHY_fivi)
summary(data_poi_XHY_fivi)

plot_fivi_XHY(data_poi_XHY_fivi)


