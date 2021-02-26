# Run simulation analysis and plots
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
library(ggplot2)
library(cowplot)
library(reshape2)
library(ggpubr)
library(simstudy)

# -----------------------------------------------------------------------------
### Example Data Analysis

# ex1: Rothenhaeusler Example -------------------------------------------------

# sim1: Rothenhaeusler Comparison ---

load("./data sets/simulation_study/ex1/2021-01-27/sim1.Rdata")

data_rot_states <- sim_data_rot$states
data_rot <- sim_data_rot$sim_data

head(data_rot)
summary(data_rot)

plot_rot(data_rot)


# sim2: Rothenhaeusler with fixed v ---

load("./data sets/simulation_study/ex1/2021-01-27/sim2.Rdata")

data_rot_X_fivi_states <- sim_data_rot_X_fivi$states
data_rot_X_fivi <- sim_data_rot_X_fivi$sim_data

head(data_rot_X_fivi)
summary(data_rot_X_fivi)

plot_rot_fivi(data_rot_X_fivi)

# sim3: Rothenhaeusler with fixed xi ---

load("./data sets/simulation_study/ex1/2021-01-27/sim3.Rdata")

# Compute data
data_rot_X_fixi <- data.frame(matrix(nrow = 0, ncol = 8))
colnames(data_rot_X_fixi) <- c("rep", "v", "xi",
                               "b", "se", "logLik_pert", "deviance", "pearson")
for (i in 1:100) {
  
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
def_poi_X_pert <- defData(def_poi_X_pert, varname = "Y", dist = "poisson",
                          link = "log",
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

eta <- b * dd$X + h * dd$H + a * dd$A
eta_pert <- b * dd_pert$X + h * dd_pert$H + a_pert * dd_pert$A

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
  bins = 30) + xlab(expression(eta)) + rremove("legend")

pdensity <- ggdensity(
  gg_hist, x = "eta", 
  color= "data", palette = c("#FF6666", "#00AFBB"),
  alpha = 0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)),
                     position = "right")  +
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

load("./data sets/simulation_study/ex2/2021-02-22/sim1.Rdata")

data_poi_X_states_fivi <- sim_data_poi_X_fivi$states
data_poi_X_fivi <- sim_data_poi_X_fivi$sim_data

head(data_poi_X_fivi)
summary(data_poi_X_fivi)

plot_fivi_X(data_poi_X_fivi)

# sim2: Poisson IV with fixed xi ---

load("./data sets/simulation_study/ex2/2021-02-18/sim2.Rdata")

nsim <- 100

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

# Choose hyperparameter of interest
poi_X_fixi_chosen <-
  data_poi_X_fixi[data_poi_X_fixi$xi %in% c(0, 1, 3, 5, 10000), ]

# Investigate different quantiles for deviance and Pearson residuals
alpha <- 0.9 # set alpha-quantile of interest for plots

quantile_data <- data.frame(matrix(nrow = 0, ncol = 6))
colnames(quantile_data) <- c("rep", "v", "xi",
                             "logLik_pert", "deviance", "pearson")

for (i in 1:nsim) {
  
  data_temp <- sim_data_poi_X_fixi[[i]][["data"]][["indv"]]
  
  vlen <- length(data_temp)
  
  for (s in 1:vlen) {
    
    # Log-Likelihood
    logLik <- apply(data_temp[[s]][[1]], 1,
                    function(x) quantile(x, probs = alpha))
    
    # Deviance Residuals
    deviance <- sqrt(apply(data_temp[[s]][[2]]^2, 1,
                           function(x) quantile(x, probs = alpha)))
    
    # Pearson Residuals
    pearson <- sqrt(apply(data_temp[[s]][[3]]^2, 1,
                          function(x) quantile(x, probs = alpha)))
    
    temp <- data.frame(rep = i, v = data_temp[[s]][[5]],
                       xi = data_temp[[s]][[4]], logLik_pert = logLik,
                       deviance = deviance, pearson = pearson)
    
    quantile_data <- rbind(quantile_data, temp)
  }
  
  data_temp$logLik
  data_poi_X_fixi <- rbind(data_poi_X_fixi, data_temp)
}

quantile_data_chosen <-
  quantile_data[quantile_data$xi %in% c(0, 1, 3, 5, 10000), ]

plot_fixi(quantile_data_chosen)

# sim3: Poisson IV identifiability ---

load("./data sets/simulation_study/ex2/2021-01-27/sim3.Rdata")

IV_b_data <- IV_b_poi$sim_data
plot_IV(IV_b_data, causal_parameter = 0.4)

# sim4: Poisson IV several quantiles ---

load("./data sets/simulation_study/ex2/2021-02-18/sim4.Rdata")

data_poi_X_quant <- sim_data_poi_X_quant$sim_data

head(data_poi_X_quant)
str(data_poi_X_quant)

q_values <- seq(0, 1, by = 0.001)
plot_quant(data_poi_X_quant, q_values, xi_big = 10000)

# ex3: Binomial with IV Example -----------------------------------------------

# Define variables for unperturbed and perturbed data set
def_bin_X <- defData(varname = "A", dist = "normal", 
                     formula = 0, variance = 0.25)
def_bin_X <- defData(def_bin_X, varname = "H", dist = "normal",
                     formula = 0, variance = 0.25)
def_bin_X <- defData(def_bin_X, varname = "X", dist = "normal", 
                     formula = "2 * H + A", variance = 0.25)
def_bin_X <- defData(def_bin_X, varname = "Y", dist = "binomial", link = "logit", 
                     formula = "0.4 * X + 1 * H", variance = 1)

def_bin_X_pert <- defData(varname = "A", dist = "normal", 
                          formula = 0, variance = 0.25)
def_bin_X_pert <- defData(def_bin_X_pert, varname = "v", 
                          formula = 3) # set perturbation strength
def_bin_X_pert <- defData(def_bin_X_pert, varname = "H", dist = "normal",
                          formula = 0, variance = 0.25)
def_bin_X_pert <- defData(def_bin_X_pert, varname = "X", dist = "normal", 
                          formula = "2 * H + v", variance = 0.25)
def_bin_X_pert <- defData(def_bin_X_pert, varname = "Y", dist = "binomial",
                          link = "logit",
                          formula = "0.4 * X + 1 * H", variance = 1)

# Parameters
b <- 0.4
h <- 1
v <- 3

# sim 0: Histograms and bar plots
set.seed(96853)
dd <- genData(1000, def_poi_X)
dd_pert <- genData(1000, def_poi_X_pert)

eta <- b * dd$X + h * dd$H
eta_pert <- b * dd_pert$X + h * dd_pert$H

p <- exp(eta) / (1 + exp(eta))
p_pert <- exp(eta_pert) / (1 + exp(eta_pert))

gg_hist <- data.frame(eta = c(eta, eta_pert),
                      p = c(p, p_pert),
                      Y = c(dd$Y,dd_pert$Y),
                      data = c(rep("unpert", length(eta)),
                               rep("pert", length(eta))))
# eta plot
phist <- gghistogram(
  gg_hist, x = "eta", 
  add = "mean", rug = FALSE,
  fill = "data", palette = c("#FF6666", "#00AFBB"),
  bins = 30) + xlab(expression(eta)) + rremove("legend")

pdensity <- ggdensity(
  gg_hist, x = "eta", 
  color= "data", palette = c("#FF6666", "#00AFBB"),
  alpha = 0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)),
                     position = "right")  +
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

# ggsave(filename = "ex3hist.pdf", plot = arr, height = 4, width = 6)


# sim1: Binomial IV with fixed v ---

load("./data sets/simulation_study/ex3/2021-02-19/sim1.Rdata")

data_bin_X_states_fivi <- sim_data_bin_X_fivi$states
data_bin_X_fivi <- sim_data_bin_X_fivi$sim_data

head(data_bin_X_fivi)
summary(data_bin_X_fivi)

plot_fivi_X(data_bin_X_fivi)

# sim2: Binomial IV with fixed xi ---

load("./data sets/simulation_study/ex3/2021-02-19/sim2.Rdata")

# Compute data
data_bin_X_fixi <- data.frame(matrix(nrow = 0, ncol = 8))
colnames(data_bin_X_fixi) <- c("rep", "v", "xi",
                               "b", "se", "logLik_pert", "deviance", "pearson")
for (i in 1:100) {
  
  data_temp <- sim_data_bin_X_fixi[[i]][["data"]][["data"]]
  
  data_bin_X_fixi <- rbind(data_bin_X_fixi, data_temp)
}

head(data_bin_X_fixi)
summary(data_bin_X_fixi)

# Investigate different quantiles of deviance and Pearson residuals

alpha <- 0.9 # set quantile

quantile_data <- data.frame(matrix(nrow = 0, ncol = 6))
colnames(quantile_data) <- c("rep", "v", "xi",
                             "logLik_pert", "deviance", "pearson")

for (i in 1:nsim) {
  
  data_temp <- sim_data_bin_X_fixi[[i]][["data"]][["indv"]]
  
  vlen <- length(data_temp)
  
  for (s in 1:vlen) {
    
    # Log-Likelihood
    logLik <- apply(data_temp[[s]][[1]], 1,
                    function(x) quantile(x, probs = alpha))
    
    # Deviance Residuals
    deviance <- sqrt(apply(data_temp[[s]][[2]]^2, 1,
                           function(x) quantile(x, probs = alpha)))
    
    # Pearson Residuals
    pearson <- sqrt(apply(data_temp[[s]][[3]]^2, 1,
                          function(x) quantile(x, probs = alpha)))
    
    temp <- data.frame(rep = i, v = data_temp[[s]][[5]],
                       xi = data_temp[[s]][[4]], logLik_pert = logLik,
                       deviance = deviance, pearson = pearson)
    
    quantile_data <- rbind(quantile_data, temp)
  }
  
  data_temp$logLik
  data_bin_X_fixi <- rbind(data_bin_X_fixi, data_temp)
}

quantile_data_chosen <-
  quantile_data[quantile_data$xi %in% c(0, 1, 3, 5, 10000), ]

plot_fixi(quantile_data_chosen)

# sim3: Binomial IV identifiability ---

load("./data sets/simulation_study/ex3/2021-02-19/sim3.Rdata")

IV_b_data <- IV_b_bin$sim_data
plot_IV(IV_b_data, causal_parameter = 1)

# sim4: Binomial IV several quantiles ---

load("./data sets/simulation_study/ex3/2021-02-19/sim4.Rdata")

data_bin_X_quant <- sim_data_bin_X_quant$sim_data

head(data_bin_X_quant)
str(data_bin_X_quant)

q_values <- seq(0, 1, by = 0.001)
plot_quant(data_bin_X_quant, q_values, xi_big = 10000)

# ex4: Poisson with Anchor on X, H and Y --------------------------------------

# Define variables for unperturbed and perturbed data set
def_poi_XHY <- defData(varname = "A", dist = "normal", 
                       formula = 0, variance = 0.25)
def_poi_XHY <- defData(def_poi_XHY, varname = "H", dist = "normal",
                       formula = "0.5 * A", variance = 0.25)
def_poi_XHY <- defData(def_poi_XHY, varname = "X", dist = "normal", 
                       formula = "2 * H + A", variance = 0.25)
def_poi_XHY <- defData(def_poi_XHY, varname = "Y", dist = "poisson",
                       link = "log", 
                       formula = "0.4 * X + 1 * H - A", variance = 1)

def_poi_XHY_pert <- defData(varname = "A", dist = "normal", 
                            formula = 0, variance = 0.25)
def_poi_XHY_pert <- defData(def_poi_XHY_pert, varname = "v", 
                            formula = 2) # set perturbation strength
def_poi_XHY_pert <- defData(def_poi_XHY_pert, varname = "H", dist = "normal",
                            formula = "0.5", variance = 0.25)
def_poi_XHY_pert <- defData(def_poi_XHY_pert, varname = "X", dist = "normal", 
                            formula = "2 * H + v", variance = 0.25)
def_poi_XHY_pert <- defData(def_poi_XHY_pert, varname = "Y", dist = "poisson",
                            link = "log",
                            formula = "0.4 * X + 1 * H - 1.5", variance = 1)

# Parameters of linear predictor
b <- 0.4
h <- 1
a <- -1

a_pert <- -1.5

# sim0: General ---

set.seed(21255)
# Distribution histogram
dd <- genData(1000, def_poi_XHY)
dd_pert <- genData(1000, def_poi_XHY_pert)

eta <- b * dd$X + h * dd$H + a * dd$A
eta_pert <- b * dd_pert$X + h * dd_pert$H + a_pert * dd_pert$A

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
  bins = 30) + xlab(expression(eta)) + rremove("legend")

pdensity <- ggdensity(
  gg_hist, x = "eta", 
  color= "data", palette = c("#FF6666", "#00AFBB"),
  alpha = 0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)),
                     position = "right")  +
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

load("./data sets/simulation_study/ex4/2021-02-21/sim1.Rdata")

data_poi_XHY_fivi_states <- sim_data_poi_XHY_fivi$states
data_poi_XHY_fivi <- sim_data_poi_XHY_fivi$sim_data

head(data_poi_XHY_fivi)
summary(data_poi_XHY_fivi)

plot_fivi_XHY(data_poi_XHY_fivi)

# sim4: Poisson IV several quantiles ---

load("./data sets/simulation_study/ex4/2021-02-19/sim4.Rdata")

data_poi_XHY_quant <- sim_data_poi_XHY_quant$sim_data

head(data_poi_XHY_quant)
str(data_poi_XHY_quant)

q_values <- seq(0, 1, by = 0.001)
plot_quant(data_poi_XHY_quant, q_values, xi_big = 10000)

# ex5: Label Noise Example ----------------------------------------------------

# Define variables for unperturbed and perturbed data set
def_LN <- defData(varname = "A", dist = "normal", 
                  formula = 0, variance = 0.01)
def_LN <- defData(def_LN, varname = "C", dist = "binomial", link = "logit",
                  formula = "-2 + A", variance = 1)
def_LN <- defData(def_LN, varname = "H", dist = "normal",
                  formula = 0, variance = 0.01)
def_LN <- defData(def_LN, varname = "X", dist = "normal", 
                  formula = "2 * H", variance = 0.01)
def_LN <- defData(def_LN, varname = "Yorg", dist = "binomial", link = "logit", 
                  formula = "-1 + 0.4 * X + 1 * H", variance = 1)
def_LN <- defData(def_LN, varname = "Y", dist = "binary", 
                  formula = "abs(Yorg - C)")

def_LN_pert <- defData(varname = "A", dist = "normal", 
                       formula = 0, variance = 0.01)
def_LN_pert <- defData(def_LN_pert, varname = "v", 
                       formula = 0.8) # set perturbation strength
def_LN_pert <- defData(def_LN_pert, varname = "C", dist = "binomial",
                       link = "logit",
                       formula = "-2 + v", variance = 1)
def_LN_pert <- defData(def_LN_pert, varname = "H", dist = "normal",
                       formula = 0, variance = 0.01)
def_LN_pert <- defData(def_LN_pert, varname = "X", dist = "normal", 
                       formula = "2 * H", variance = 0.01)
def_LN_pert <- defData(def_LN_pert, varname = "Yorg", dist = "binomial",
                       link = "logit", 
                       formula = "-1 + 0.4 * X + 1 * H", variance = 1)
def_LN_pert <- defData(def_LN_pert, varname = "Y", dist = "binary", 
                       formula = "abs(Yorg - C)")

# Parameters of linear predictor
b <- 0.4
h <- 1
c <- -1

t <- -2
v_pert <- 0.8

# sim0: Histograms and bar plots
set.seed(68465)
dd <- genData(1000, def_LN)
dd_pert <- genData(1000, def_LN_pert)

eta <- c + b * dd$X + h * dd$H
eta_pert <- c + b * dd_pert$X + h * dd_pert$H

tau <- t + dd$A
tau_pert <- t + v_pert
exp(tau_pert) / (1 + exp(tau_pert))

gg_hist <- data.frame(tau = tau,
                      mu = exp(tau) / (1 + exp(tau)),
                      C = dd$C,
                      data = rep("unpert", length(tau)))

gg_bar <- data.frame(C = c(dd$C, dd_pert$C),
                     data = c(rep("unpert", length(tau)),
                              rep("pert", length(tau))))
# tau plot
phist <- gghistogram(
  gg_hist, x = "mu", 
  add = "mean", rug = FALSE,
  fill = "data", palette = c("#00AFBB"),
  bins = 30) + xlab(expression(p[C])) + rremove("legend")

pdensity <- ggdensity(
  gg_hist, x = "mu", 
  color= "data", palette = c( "#00AFBB"),
  alpha = 0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)),
                     position = "right")  +
  theme_half_open(11, rel_small = 1) +
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("legend") +
  annotate("text", label= expression(p[C]^pert), colour = 1, x = 0.143, y = 25) +
  annotate("text", label= "= 0.231", colour = 1, x = 0.1503, y = 24.8)

aligned_plots <- align_plots(phist, pdensity, align="hv", axis="tblr")
gg_tau <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

# C plot
gg_C <- ggplot(gg_bar, aes(C, fill = data)) + 
  geom_bar(position = 'identity', alpha = .4) +
  scale_fill_manual("data", values = c("pert" = "#FF6666", "unpert" = "#00AFBB"))

arr <- ggarrange(gg_tau, gg_C,
                 labels = c("A", "B"),
                 font.label = list(face = "plain"),
                 ncol = 1, nrow = 2, common.legend = TRUE)
arr

#ggsave(filename = "exLNhist.pdf", plot = arr, height = 4, width = 6)

# sim1: Label Noise with fixed v ---

load("./data sets/simulation_study/ex5/2021-02-18/sim1.Rdata")

data_LN_states_fivi <- sim_data_LN_fivi$states
data_LN_fivi <- sim_data_LN_fivi$sim_data

head(data_LN_fivi)
summary(data_LN_fivi)

plot_fivi_X(data_LN_fivi)

# sim2: Label Noise with fixed xi ---

load("./data sets/simulation_study/ex5/2021-02-03/sim2.Rdata")

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

# Plot with glm(Y~X+A-1)
data_LN_fixi_glm <- data.frame(matrix(nrow = 0, ncol = 8))
colnames(data_LN_fixi_glm) <- c("rep", "v", "xi",
                                "b", "se", "logLik_pert", "deviance", "pearson")
for (i in 1:nsim) {
  
  data_temp <- sim_data_LN_fixi_glm[[i]][["data"]][["data"]]
  
  data_LN_fixi_glm <- rbind(data_LN_fixi_glm, data_temp)
}

head(data_LN_fixi_glm)
summary(data_LN_fixi_glm)

data_LN_fixi_glm$xi <- "glm"

data_LN_fixi_all <- rbind(data_LN_fixi_glm, data_LN_fixi)

plot_fixi(data_LN_fixi_all)

# sim4: Label Noise several quantiles ---

load("./data sets/simulation_study/ex5/2021-02-18/sim4.Rdata")

data_LN_quant <- sim_data_LN_quant$sim_data

head(data_LN_quant)
str(data_LN_quant)

q_values <- seq(0, 1, by = 0.01)
plot_quant(data_LN_quant, q_values, xi_big = 10000)

# ex6: Label Noise Example ----------------------------------------------------

# Define variables for unperturbed and perturbed data set
def_LN <- defData(varname = "A", dist = "normal", 
                  formula = 0, variance = 0.01)
def_LN <- defData(def_LN, varname = "B", dist = "binomial", link = "logit",
                  formula = "-2 + A", variance = 1)
def_LN <- defData(def_LN, varname = "H", dist = "normal",
                  formula = "A", variance = 0.01)
def_LN <- defData(def_LN, varname = "X", dist = "normal", 
                  formula = "2 * H + A", variance = 0.01)
def_LN <- defData(def_LN, varname = "Yorg", dist = "binomial", link = "logit", 
                  formula = "-1 + 0.4 * X + 1 * H", variance = 1)
def_LN <- defData(def_LN, varname = "Y", dist = "binary", 
                  formula = "abs(Yorg - B)")

def_LN_pert <- defData(varname = "A", dist = "normal", 
                       formula = 0, variance = 0.01)
def_LN_pert <- defData(def_LN_pert, varname = "v", 
                       formula = 0.6) # set perturbation strength
def_LN_pert <- defData(def_LN_pert, varname = "B", dist = "binomial",
                       link = "logit",
                       formula = "-2 + v", variance = 1)
def_LN_pert <- defData(def_LN_pert, varname = "vH", 
                       formula = 0.4) # set perturbation strength
def_LN_pert <- defData(def_LN_pert, varname = "H", dist = "normal",
                       formula = "vH", variance = 0.01)
def_LN_pert <- defData(def_LN_pert, varname = "vX", 
                       formula = -0.5) # set perturbation strength
def_LN_pert <- defData(def_LN_pert, varname = "X", dist = "normal", 
                       formula = "2 * H + vX", variance = 0.01)
def_LN_pert <- defData(def_LN_pert, varname = "Yorg", dist = "binomial",
                       link = "logit", 
                       formula = "-1 + 0.4 * X + 1 * H", variance = 1)
def_LN_pert <- defData(def_LN_pert, varname = "Y", dist = "binary", 
                       formula = "abs(Yorg - B)")

# Parameters of linear predictor
b <- 0.4
h <- 1
c <- -1

t <- -2
v_pert <- 0.6

# sim0: Histograms and bar plots
dd <- genData(1000, def_LN)
dd_pert <- genData(1000, def_LN_pert)

eta <- c + b * dd$X + h * dd$H
eta_pert <- c + b * dd_pert$X + h * dd_pert$H

tau <- t + dd$A
tau_pert <- t + v_pert

# sim0: histograms and bar plots
gg_hist <- data.frame(tau = tau,
                      mu = exp(tau)/(1+exp(tau)),
                      data = rep("unpert", length(tau)))

gg_bar <- data.frame(B = c(dd$B, dd_pert$B),
                     data = c(rep("unpert", length(tau)),
                              rep("pert", length(tau))))
# tau plot
phist <- gghistogram(
  gg_hist, x = "mu", 
  add = "mean", rug = FALSE,
  fill = "data", palette = c("#00AFBB"),
  bins = 30
) + xlab(expression(p[C])) + rremove("legend")

pdensity <- ggdensity(
  gg_hist, x = "mu", 
  color= "data", palette = c( "#00AFBB"),
  alpha = 0,
) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)),
                     position = "right")  +
  theme_half_open(11, rel_small = 1) +
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("legend") +
  annotate(geom = 'text', x = 0.145, y = 25,
             label = "p[C]^pert %~~% 0.198",
             parse = TRUE)

aligned_plots <- align_plots(phist, pdensity, align="hv", axis="tblr")
gg_tau <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

# B plot
gg_C <- ggplot(gg_bar, aes(B, fill = data)) + 
  geom_bar(position = 'identity', alpha = .4) +
  scale_fill_manual("data", values = c("pert" = "#FF6666", "unpert" = "#00AFBB"))

arr <- ggarrange(gg_tau, gg_C,
                 labels = c("A", "B"),
                 font.label = list(face = "plain"),
                 ncol = 1, nrow = 2, common.legend = TRUE)
arr

#ggsave(filename = "exLNhist.pdf", plot = arr, height = 4, width = 6)

# sim1: Label Noise with fixed v ---

load("./data sets/simulation_study/ex6/2021-02-17/sim1.Rdata")

data_LN_states_fivi <- sim_data_LN_fivi$states
data_LN_fivi <- sim_data_LN_fivi$sim_data

head(data_LN_fivi)
summary(data_LN_fivi)

plot_fivi_X(data_LN_fivi)

# sim4: Label Noise several quantiles ---

load("./data sets/simulation_study/ex6/2021-02-17/sim4.Rdata")

data_LN_quant <- sim_data_LN_quant$sim_data

head(data_LN_quant)
str(data_LN_quant)

q_values <- seq(0, 1, by = 0.01)
plot_quant(data_LN_quant, q_values, xi_big = 10000)

# Plot with glm(Y~X+A-1)

load("./data sets/simulation_study/ex6/2021-02-17/sim4_glm.Rdata")

data_LN_quant_glm <- sim_data_LN_quant_glm$sim_data

q_len <- length(q_values)

quantile_array <- matrix(nrow = nsim, ncol = q_len)

for (i in 1:nsim) {
  for (q in 1:q_len) {
    quant <- q_values[q]
    
    quantile_array[i, q] <- 
      as.numeric(quantile(data_LN_quant_glm[[i]][[2]][, 2], quant))
    
  }
}
glm_quantile_vector <- apply(quantile_array, 2, mean)

plot_quant(data_LN_quant, q_values, xi_big = 10000, glm_quantile_vector)

