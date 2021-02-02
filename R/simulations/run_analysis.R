# Run simulation analysis and plots
#
# Date: 19.01.21     Author: Maic Rakitta
###############################################################################

library(ggplot2)
library(cowplot)
library(reshape2)
library(ggpubr)
library(simstudy)

# -----------------------------------------------------------------------------
### Example Data Analysis

# ex1: Rothenhaeusler Example -------------------------------------------------

# sim1: Rothenhaeusler Comparison ---

load(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/simulation_study/ex1/2021-01-27/sim1.Rdata")

data_rot_states <- sim_data_rot$states
data_rot <- sim_data_rot$sim_data

head(data_rot)
summary(data_rot)

plot_rot(data_rot)


# sim2: Rothenhaeusler with fixed v ---

load(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/simulation_study/ex1/2021-01-27/sim2.Rdata")

data_rot_X_fivi_states <- sim_data_rot_X_fivi$states
data_rot_X_fivi <- sim_data_rot_X_fivi$sim_data

head(data_rot_X_fivi)
summary(data_rot_X_fivi)

plot_rot_fivi(data_rot_X_fivi)

# sim3: Rothenhaeusler with fixed xi ---

load(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/simulation_study/ex1/2021-01-27/sim3.Rdata")

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

# ex2: Poisson with IV Example ------------------------------------------------

# Define variables for unperturbed and perturbed data set
def_poi_X <- defData(varname = "A", dist = "normal", 
                     formula = 0, variance = 0.25)
def_poi_X <- defData(def_poi_X, varname = "H", dist = "normal",
                     formula = 0, variance = 0.25)
def_poi_X <- defData(def_poi_X, varname = "X", dist = "normal", 
                     formula = "2 * H + A", variance = 0.25)
def_poi_X <- defData(def_poi_X, varname = "Y", dist = "poisson", link = "log", 
                     formula = "0.4 * X + 1 * H", variance = 1)

def_poi_X_pert <- defData(varname = "A", dist = "normal", 
                          formula = 0, variance = 0.25)
def_poi_X_pert <- defData(def_poi_X_pert, varname = "v", 
                          formula = 3) # set perturbation strength
def_poi_X_pert <- defData(def_poi_X_pert, varname = "H", dist = "normal",
                          formula = 0, variance = 0.25)
def_poi_X_pert <- defData(def_poi_X_pert, varname = "X", dist = "normal", 
                          formula = "2 * H + v", variance = 0.25)
def_poi_X_pert <- defData(def_poi_X_pert, varname = "Y", dist = "poisson", link = "log",
                          formula = "0.4 * X + 1 * H", variance = 1)

# Parameters of linear predictor
b <- 0.4
h <- 1
a <- 0
a_pert <- 3

# sim0: General ---

set.seed(8151)
# Distribution histogram
dd <- genData(1000, def_poi_X)
dd_pert <- genData(1000, def_poi_X_pert)

eta <- b*dd$X + h * dd$H + a * dd$A
eta_pert <- b*dd_pert$X + h * dd_pert$H + a_pert * dd_pert$A

gg_hist <- data.frame(eta = c(eta, eta_pert),
                      mu = c(exp(eta), exp(eta_pert)),
                      Y = c(dd$Y,dd_pert$Y),
                      data = c(rep("unpert", length(eta)),
                               rep("pert", length(eta))))
# eta plot
phist <- gghistogram(
  gg_hist, x = "eta", 
  add = "mean", rug = FALSE,
  fill = "data", palette = c("#FF6666", "#00AFBB"),
  bins = 30
) + xlab(expression(eta)) + rremove("legend")

pdensity <- ggdensity(
  gg_hist, x = "eta", 
  color= "data", palette = c("#FF6666", "#00AFBB"),
  alpha = 0,
) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), position = "right")  +
  theme_half_open(11, rel_small = 1) +
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("legend")

aligned_plots <- align_plots(phist, pdensity, align="hv", axis="tblr")
gg_eta <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

# Y plot
gg_Y <- ggplot(gg_hist, aes(Y, fill = data)) + 
  geom_bar(position = 'identity', alpha = .4) +
  scale_fill_manual("data", values = c("pert" = "#FF6666", "unpert" = "#00AFBB"))

arr <- ggarrange(gg_eta, gg_Y,
                 labels = c("A", "B"),
                 font.label = list(face = "plain"),
                 ncol = 1, nrow = 2, common.legend = TRUE)
arr

# ggsave(filename = "ex2hist.pdf", plot = arr, height = 4, width = 6)

# sim1: Poisson IV with fixed v ---

load(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/simulation_study/ex2/2021-01-27/sim1.Rdata")

data_poi_X_states_fivi <- sim_data_poi_X_fivi$states
data_poi_X_fivi <- sim_data_poi_X_fivi$sim_data

head(data_poi_X_fivi)
summary(data_poi_X_fivi)

plot_fivi_X(data_poi_X_fivi)

# sim2: Poisson IV with fixed xi ---

load(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/simulation_study/ex2/2021-01-27/sim2.Rdata")


# Compute data
data_poi_X_fixi <- data.frame(matrix(nrow = 0, ncol = 8))
colnames(data_poi_X_fixi) <- c("rep", "v", "xi",
                               "b", "se", "logLik_pert", "deviance", "pearson")

nsim <- 100

for (i in 1:nsim) {
  
  data_temp <- sim_data_poi_X_fixi[[i]][["data"]][["data"]]
  
  data_poi_X_fixi <- rbind(data_poi_X_fixi, data_temp)
}

head(data_poi_X_fixi)
summary(data_poi_X_fixi)

# Plot Quantile of log-likelihood and with residuals
poi_X_fixi_chosen <-
  data_poi_X_fixi[data_poi_X_fixi$xi %in% c(0, 1, 3, 5, 10000), ]

plot_fixi(poi_X_fixi_chosen)

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
    # deviance <- apply(data_temp[[s]][[2]], 1, function(x) quantile(x, probs = alpha))
    deviance <- sqrt(apply(data_temp[[s]][[2]]^2, 1, function(x) quantile(x, probs = alpha)))
    
    # Pearson Residuals
    # pearson <- apply(data_temp[[s]][[3]], 1, function(x) quantile(x, probs = alpha))
    pearson <- sqrt(apply(data_temp[[s]][[3]]^2, 1, function(x) quantile(x, probs = alpha)))
    
    temp <- data.frame(rep = i, v = data_temp[[s]][[5]], xi = data_temp[[s]][[4]],
                       logLik_pert = logLik, deviance = deviance, pearson = pearson)
    
    quantile_data <- rbind(quantile_data, temp)
  }
  
  data_temp$logLik
  
  data_poi_X_fixi <- rbind(data_poi_X_fixi, data_temp)
}

quantile_data_chosen <-
  quantile_data[quantile_data$xi %in% c(0, 1, 3, 5, 10000), ]

plot_fixi(quantile_data_chosen)

# sim3: Poisson IV identifiability ---

load(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/simulation_study/ex2/2021-01-27/sim3.Rdata")

IV_b_data <- IV_b_poi$sim_data
plot_IV(IV_b_data, causal_parameter = 0.4)

# sim4: Poisson IV several quantiles ---

load(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/simulation_study/ex2/2021-01-27/sim4.Rdata")

data_poi_X_quant <- sim_data_poi_X_quant$sim_data

head(data_poi_X_quant)
str(data_poi_X_quant)

q_values <- seq(0, 1, by = 0.01)
plot_quant(data_poi_X_quant, q_values, xi_big = 10000)



# ex3: Binomial with IV Example -----------------------------------------------

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

# Define variables for unperturbed and perturbed data set
def_poi_XHY <- defData(varname = "A", dist = "normal", variance = 0.25,
                       formula = 0)
def_poi_XHY <- defData(def_poi_XHY, varname = "H", dist = "normal",
                       variance = 0.25,
                       formula = "0.6 + A",)
def_poi_XHY <- defData(def_poi_XHY, varname = "X", dist = "normal",
                       variance = 0.25,
                       formula = "H + A")
def_poi_XHY <- defData(def_poi_XHY, varname = "Y",
                       dist = "poisson", link = "log", 
                       formula = "0.8 * X - H - A")

def_poi_XHY_pert <- defData(varname = "A", dist = "normal", 
                            variance = 0.25,
                            formula = 0)
def_poi_XHY_pert <- defData(def_poi_XHY_pert, varname = "H", dist = "normal",
                            variance = 0.25,
                            formula = "0.6 + 0.2 * A") # set perturbation
def_poi_XHY_pert <- defData(def_poi_XHY_pert, varname = "X", dist = "normal", 
                            variance = 0.25,
                            formula = "H - 1 * A") # set perturbation 
def_poi_XHY_pert <- defData(def_poi_XHY_pert, varname = "Y",
                            dist = "poisson", link = "log", 
                            formula = "0.8 * X - H + 2 * A") # set perturbation 

# Parameters of linear predictor
b <- 0.8
h <- -1
a <- -1

a_pert <- 2

# sim0: General ---

set.seed(21255)
# Distribution histogram
dd <- genData(1000, def_poi_XHY)
dd_pert <- genData(1000, def_poi_XHY_pert)

eta <- b*dd$X + h * dd$H + a * dd$A
eta_pert <- b*dd_pert$X + h * dd_pert$H + a_pert * dd_pert$A

gg_hist <- data.frame(eta = c(eta, eta_pert),
                      mu = c(exp(eta), exp(eta_pert)),
                      Y = c(dd$Y,dd_pert$Y),
                      data = c(rep("unpert", length(eta)),
                               rep("pert", length(eta))))
# eta plot
phist <- gghistogram(
  gg_hist, x = "eta", 
  add = "mean", rug = FALSE,
  fill = "data", palette = c("#FF6666", "#00AFBB"),
  bins = 30
) + xlab(expression(eta)) + rremove("legend")

pdensity <- ggdensity(
  gg_hist, x = "eta", 
  color= "data", palette = c("#FF6666", "#00AFBB"),
  alpha = 0,
) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), position = "right")  +
  theme_half_open(11, rel_small = 1) +
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("legend")

aligned_plots <- align_plots(phist, pdensity, align="hv", axis="tblr")
gg_eta <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

# Y plot
gg_Y <- ggplot(gg_hist, aes(Y, fill = data)) + 
  geom_bar(position = 'identity', alpha = .4) +
  scale_fill_manual("data", values = c("pert" = "#FF6666", "unpert" = "#00AFBB"))

arr <- ggarrange(gg_eta, gg_Y,
                 labels = c("A", "B"),
                 font.label = list(face = "plain"),
                 ncol = 1, nrow = 2, common.legend = TRUE)
arr

#ggsave(filename = "ex4hist.pdf", plot = arr, height = 4, width = 6)


# sim1: Poisson with Anchor on X, H and Y fixed v ---

load(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/simulation_study/ex4/2021-01-27/sim1.Rdata")

data_poi_XHY_fivi_states <- sim_data_poi_XHY_fivi$states
data_poi_XHY_fivi <- sim_data_poi_XHY_fivi$sim_data

head(data_poi_XHY_fivi)
summary(data_poi_XHY_fivi)

plot_fivi_XHY(data_poi_XHY_fivi)

# sim4: Poisson IV several quantiles ---

load("C:/Users/maicr/Desktop/Github/MScAnchor/data sets/simulation_study/ex4/2021-01-27/sim4.Rdata")

data_poi_XHY_quant <- sim_data_poi_XHY_quant$sim_data

head(data_poi_XHY_quant)
str(data_poi_XHY_quant)

q_values <- seq(0, 1, by = 0.01)
plot_quant(data_poi_XHY_quant, q_values, xi_big = 10000)


# ex5: Label Noise Example ----------------------------------------------------

# sim1: Label Noise with fixed v ---

load(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/simulation_study/ex5/2021-01-29/sim1.Rdata")

data_LN_states_fivi <- sim_data_LN_fivi$states
data_LN_fivi <- sim_data_LN_fivi$sim_data

head(data_LN_fivi)
summary(data_LN_fivi)

plot_fivi_X(data_LN_fivi)

# sim2: Label Noise with fixed xi ---

load(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/simulation_study/ex5/2021-01-29/sim2.Rdata")

# Compute data
data_LN_fixi <- data.frame(matrix(nrow = 0, ncol = 8))
colnames(data_LN_fixi) <- c("rep", "v", "xi",
                               "b", "se", "logLik_pert", "deviance", "pearson")
for (i in 1:nsim) {
  
  data_temp <- sim_data_LN_fixi[[i]][["data"]][["data"]]
  
  data_LN_fixi <- rbind(data_LN_fixi, data_temp)
}

head(data_LN_fixi)
summary(data_LN_fixi)

# Plot Quantile of log-likelihood and with residuals
plot_fixi(data_LN_fixi)

# sim4: Label Noise several quantiles ---

load(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/simulation_study/ex5/2021-01-29/sim4.Rdata")

data_LN_quant <- sim_data_LN_quant$sim_data

head(data_LN_quant)
str(data_LN_quant)

q_values <- seq(0, 1, by = 0.01)
plot_quant(data_LN_quant, q_values, xi_big = 10000)

