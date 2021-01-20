
library(glare)
library(simstudy)
library(extrafont)
library(tidyverse)
library(rsimsum)
library(ggplot2)
library(gridExtra)
library(reshape2)

# NSIM = 100 fuer fixi ----------------------
# Function definition for fixed xi --------------------------------------------

onerep_fixi <- function(rep, nobs = 300, data_table, data_pert_table, 
                        formula, A_formula, xi_values, v_values,
                        family, type, quant_value = 0.9) {
  
  v_len <- length(v_values)
  
  data_frame <- data.frame()
  
  for (s in 1:v_len) {
    data_pert_table_temp <- updateDef(data_pert_table,
                                      changevar = "v",
                                      newformula = v_values[s])
    
    # Generate data with nobs observations
    set.seed(rep) # use the same seeds for each v to use same errors
    dd <- genData(nobs, data_table)
    set.seed(rep)
    dd_pert <- genData(nobs, data_pert_table_temp)
    
    # Fit glares
    xi_len <- length(xi_values)
    b <- matrix(nrow = xi_len, ncol = 1)
    b_se <- matrix(nrow = xi_len, ncol = 1)
    logLik_pert_indiv <- matrix(nrow = xi_len, ncol = nobs)
    logLik_pert <- matrix(nrow = xi_len, ncol = 1)
    
    devres <- matrix(nrow = xi_len, ncol = nobs)
    deviance <- matrix(nrow = xi_len, ncol = 1)
    
    peares <- matrix(nrow = xi_len, ncol = nobs)
    pearson <- matrix(nrow = xi_len, ncol = 1)
    
    for (j in 1:xi_len) {
      xi <- xi_values[j]
      
      fit <- glare(formula = formula,
                   A_formula = A_formula,
                   data = dd,
                   xi = xi,
                   family = family,
                   type = type)
      
      b[j] <- as.numeric(coef(fit))
      b_se[j] <- fit$coef_se
      
      logLik_pert_indiv[j, ] <- logLik(fit, newdata = dd_pert, indiv = TRUE)
      logLik_pert[j] <- as.numeric(quantile(logLik_pert_indiv[j, ], quant_value))
      
      if(family == "poisson") {
        
        link_inv <- poisson()$linkinv
        vari <- poisson()$variance
        mu <- link_inv(dd_pert$X * b[j])
        
        r <- mu
        p <- which(dd_pert$Y > 0)
        r[p] <- (dd_pert$Y * log(dd_pert$Y / mu) - (dd_pert$Y - mu))[p]
        
        devres[j, ] <- sign(dd_pert$Y - mu) * sqrt(2 * r)
        
      } else if(family == "binomial") {
        
        link_inv <- binomial()$linkinv
        vari <- binomial()$variance
        mu <- link_inv(dd_pert$X * b[j])
        
        r <- dd_pert$Y * log(dd_pert$Y / (mu)) + (1 - dd_pert$Y) *
          log((1 - dd_pert$Y) / (1 - mu))
        p <- which(dd_pert$Y == 0)
        r[p] <- ((1 - dd_pert$Y) * log((1 - dd_pert$Y) / (1 - mu)))[p]
        q <- which(dd_pert$Y == 1)
        r[q] <- (dd_pert$Y * log(dd_pert$Y / (mu)))[q]
        
        devres[j, ] <- sign(dd_pert$Y / 1 - mu) * sqrt(2 * r)
        
      } else if(family == "gaussian"){
        
        link_inv <- gaussian()$linkinv
        vari <- gaussian()$variance
        mu <- link_inv(dd_pert$X * b[j])
        
        n <- length(dd_pert$Y)
        s <- sqrt(1/n * sum((dd_pert$Y-mu)^2))
        
        devres[j, ] <- 1 / s * (dd_pert$Y - mu)
        
      } else {
        stop("No valid family! Use either 'poisson', 'binomial' or 'gaussian'.")
      }
      
      # Deviance
      deviance[j] <- sqrt(mean(devres[j, ]^2))
      
      # Pearson
      V <- vari(mu)
      peares[j, ] <- (dd_pert$Y - mu)/sqrt(V)
      pearson[j] <-  sqrt(mean(peares[j, ]^2))
    }
    
    # Return
    data_frame_temp <- data.frame(rep = rep, v = v_values[s], xi = xi_values,
                                  b = b,
                                  se = b_se,
                                  logLik_pert = logLik_pert,
                                  deviance = deviance,
                                  pearson = pearson)
    
    data_frame <- rbind(data_frame, data_frame_temp)
  }
  
  list(data = data_frame,
       logLik = logLik_pert_indiv,
       devres = devres,
       peares = peares)
}


# Define simulation function
simulate_fixi <- function(nsim = 100, nobs = 300, xi_values, v_values, 
                          data_table, data_pert_table, 
                          formula, A_formula, family, type,
                          quant_value = 0.9) {
  
  data_list <- list()
  
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  for (r in 1:nsim) {
    
    data_temp <-
      onerep_fixi(rep = r, nobs = nobs,
                  data_table = data_table,
                  data_pert_table = data_pert_table, 
                  formula = formula, A_formula = A_formula,
                  xi_values = xi_values, v_values = v_values,
                  family = family, type = type,
                  quant_value = quant_value)
    
    data_list[[r]] <- list(rep = r,
                           data = data_temp)
    
    Sys.sleep(0.1)
    setTxtProgressBar(pb, r)
  }
  close(pb)
  
  data_list
}
# Function definition for IV parameter identifiability ------------------------

onerep_IV <- function(rep, nobs_values = 300, data_table, 
                      formula, A_formula, xi_values,
                      family, type) {
  
  data <- as.data.frame(matrix(ncol = 5, nrow = 0))
  colnames(data) <- c("rep", "nobs", "xi", "b", "b_se")
  
  for (k in 1:length(nobs_values)) {
    
    nobs <- nobs_values[k]
    
    # Generate data set with nobs observations
    dd <- genData(nobs, data_table)
    # dd_pert <- genData(nobs, data_pert_table)
    
    xi_len <- length(xi_values)
    
    b <- matrix(nrow = xi_len, ncol = 1)
    b_se <- matrix(nrow = xi_len, ncol = 1)
    
    # logLik_pert <- numeric(xi_len)
    
    for (i in 1:xi_len) {
      
      xi <- xi_values[i]
      
      fit_temp <- glare(formula = formula,
                        A_formula = A_formula,
                        data = dd,
                        xi = xi,
                        family = family,
                        type = type)
      
      b[i] <- as.numeric(coef(fit_temp))
      b_se[i] <- fit_temp$coef_se
      
      # logLik_pert_indiv <- logLik(fit_temp, newdata = dd_pert, indiv = TRUE)
      # logLik_pert[i] <-
      #   as.numeric(quantile(logLik_pert_indiv, quant_value))
    }
    
    # Return
    data_temp <- data.frame(rep = rep, nobs = nobs,
                            xi = xi_values,
                            b = b,
                            b_se = b_se)
    
    data <- rbind(data, data_temp)
  }
  
  data
}

# Define simulation function
simulate_IV <- function(nsim, nobs_values = 300, xi_values,
                        data_table, formula, A_formula, family, type) {
  
  xi_len <- length(xi_values)
  
  states <- array(dim = c(nsim, 626))
  
  b_data <- as.data.frame(matrix(ncol = 5, nrow = 0))
  colnames(b_data) <- c("rep", "nobs", "xi", "b", "b_se")
  
  
  # 
  # sim_data <- data.frame(matrix(nrow = nsim * (xi_len + 2), ncol = 5))
  # colnames(sim_data) <- c("rep", "xi",
  #                         "b", "se", "logLik_pert")
  # 
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  for (r in 1:nsim) {
    states[r, ] <- .Random.seed
    
    data_temp <- onerep_IV(rep = r, nobs_values = nobs_values,
                           data_table = data_table, 
                           formula = formula, A_formula = A_formula,
                           xi_values = xi_values,
                           family = family, type = type)
    
    b_data <- rbind(b_data, data_temp)
    
    # sim_data[((xi_len + 2) * (r - 1) + + 1):((xi_len + 2) * r), ] <- 
    #   onerep_fivi(rep = r,
    #               data_table = data_table, data_pert_table = data_pert_table, 
    #               formula = formula, A_formula = A_formula,
    #               xi_values = xi_values, xi_big = xi_big,
    #               family = family, type = type,
    #               quant_value = quant_value)
    
    Sys.sleep(0.1)
    setTxtProgressBar(pb, r)
  }
  close(pb)
  
  list(states = states, sim_data = b_data)
}



# Plot function for fixed xi --------------------------------------------------
plot_fixi2 <- function(sim_data) {
  
  mean_data <- sim_data %>%
    group_by(xi, v) %>%
    summarise(mean_logLik_pert = mean(logLik_pert),
              mean_dev = mean(deviance),
              mean_pea = mean(pearson))
  
  gg_data <- mean_data
  gg_data$xi_val <- as.character(mean_data$xi) 
  
  gg1 <- ggplot(gg_data, aes(y = mean_logLik_pert, x = v, color = xi_val, group = xi_val)) +
    
    geom_line() +
    
    labs(color = expression(xi)) +
    ylab("90th percentile of log-likelihood") +
    xlab("v")
  
  gg_data2 <- gg_data[gg_data$v <= 3 & gg_data$v >= -3, ]
  
  gg2 <- ggplot(gg_data2, aes(y = mean_logLik_pert, x = v, color = xi_val, group = xi_val)) +
    
    geom_line() +
    
    labs(color = expression(xi)) +
    ylab("90th percentile of log-likelihood") +
    xlab("v")
  
  #deviance
  gg3 <- ggplot(gg_data, aes(y = mean_dev, x = v, color = xi_val, group = xi_val)) +
    
    geom_line() +
    
    labs(color = expression(xi)) +
    ylab("Averaged root of averaged squared deviance residuals") +
    xlab("v")
  
  gg4 <- ggplot(gg_data2, aes(y = mean_dev, x = v, color = xi_val, group = xi_val)) +
    
    geom_line() +
    
    labs(color = expression(xi)) +
    ylab("Averaged root of averaged squared deviance residuals") +
    xlab("v")
  
  #pearson
  gg5 <- ggplot(gg_data, aes(y = mean_pea, x = v, color = xi_val, group = xi_val)) +
    
    geom_line() +
    
    labs(color = expression(xi)) +
    ylab("Averaged root of averaged squared pearson residuals") +
    xlab("v")
  
  gg6 <- ggplot(gg_data2, aes(y = mean_pea, x = v, color = xi_val, group = xi_val)) +
    
    geom_line() +
    
    labs(color = expression(xi)) +
    ylab("Averaged root of averaged squared pearson residuals") +
    xlab("v")
  
  print(gg1)
  print(gg2)
  print(gg3)
  print(gg4)
  print(gg5)
  print(gg6)
  invisible(gg_data)
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
    
    xlab("Number of Observation") +
    ylab("Causal Parameter Estimate")
  
  
  # xi big with CI, nobs all
  gg2 <- ggplot(data = gg_data_big, aes(y = mean_b, x = nobs, group = 1)) +
    geom_point() +
    geom_line() + 
    geom_errorbar(aes(x = nobs, ymin = lower, ymax = upper)) +
    
    geom_hline(yintercept = causal_parameter, linetype = "dashed") +
    
    xlab("Number of Observation") +
    ylab("Causal Parameter Estimate")
  
  # xi big with CI, nobs 100-5000  
  gg_data_big2 <- gg_data_big[gg_data_big$nobs >= 100 & gg_data_big$nobs <= 5000, ]
  gg3 <- ggplot(data = gg_data_big2, aes(y = mean_b, x = nobs, group = 1)) +
    geom_point() +
    geom_line() + 
    geom_errorbar(aes(x = nobs, ymin = lower, ymax = upper)) +
    
    geom_hline(yintercept = causal_parameter, linetype = "dashed") +
    
    xlab("Number of Observation") +
    ylab("Causal Parameter Estimate")
  
  # xi all with CI, nobs all
  gg_data_all <- mean_data
  gg_data_all$xi <- as.factor(gg_data_all$xi)
  gg4 <- ggplot(data = gg_data_all, aes(y = mean_b, x = nobs, group = xi, color = xi)) +
    geom_point() +
    geom_line() + 
    geom_errorbar(aes(x = nobs, ymin = lower, ymax = upper)) +
    
    geom_hline(yintercept = causal_parameter, linetype = "dashed") +
    
    xlab("Number of Observation") +
    ylab("Causal Parameter Estimate") +
    labs(color = expression(xi))
  
  # xi all with CI, nobs 100-5000
  gg_data_all2 <- gg_data_all[gg_data_all$nobs >= 100 & gg_data_all$nobs <= 5000, ]
  gg5 <- ggplot(data = gg_data_all2, aes(y = mean_b, x = nobs, group = xi, color = xi)) +
    geom_point() +
    geom_line() + 
    geom_errorbar(aes(x = nobs, ymin = lower, ymax = upper)) +
    
    geom_hline(yintercept = causal_parameter, linetype = "dashed") +
    
    xlab("Number of Observation") +
    ylab("Causal Parameter Estimate") +
    labs(color = expression(xi))
  
  print(gg1)
  print(gg2)
  print(gg3)
  print(gg4)
  print(gg5)
  
  invisible(mean_data)
  
}


# EX2: poisson ---------------------------------------------------------------

Define variables for unperturbed and perturbed data set
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
                          formula = 4) # set perturbation strength
def_poi_X_pert <- defData(def_poi_X_pert, varname = "H", dist = "normal",
                          formula = 0, variance = 0.25)
def_poi_X_pert <- defData(def_poi_X_pert, varname = "X", dist = "normal", 
                          formula = "2 * H + v", variance = 0.25)
def_poi_X_pert <- defData(def_poi_X_pert, varname = "Y", dist = "poisson", link = "log",
                          formula = "0.4 * X + 1 * H", variance = 1)

b <- 0.4
h <- 1

dd <- genData(300, def_poi_X)
dd_pert <- genData(300, def_poi_X_pert)
par(mfrow = c(2,2))
hist(b*dd$X + h * dd$H, breaks = 100)
hist(dd$Y, breaks = 100)
hist(b*dd_pert$X + h * dd_pert$H, breaks = 100)
hist(dd_pert$Y, breaks = 100)

type <- "deviance"

# library(Hmisc)
#  hist.data.frame(dd[ , c("Y", "X", "A", "H")])

# sim2: Poisson IV for fixed xi ---

# Initialize for fixed xi
set.seed(96574)
nsim <- 10

xi_values <- c(0, 1, 50, 10000)
v_values <- seq(-10, 10, by = 1)

data_table <- def_poi_X
data_pert_table <- def_poi_X_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- "poisson"


sim_data_poi_X_fixi <- simulate_fixi(nsim = nsim, nobs = 300,
                                     xi_values = xi_values,
                                     v_values = v_values, 
                                     data_table = data_table,
                                     data_pert_table = data_pert_table, 
                                     formula = formula, A_formula = A_formula,
                                     family = family, type = type,
                                     quant_value = 0.9)


# # Save simulated data list
# path_name <- "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/simulation_study/ex2/"
# dir.create(paste(path_name, Sys.Date(), sep = ""))
# save(sim_data_poi_X_fixi, file = paste(paste(path_name, Sys.Date(), sep = ""), "/sim2_smooth.Rdata", sep =""))



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

plot_fixi2(data_poi_X_fixi)


# sim3: IV investigation for poisson ---

# Initialize for fixed v
set.seed(88435)
# nobs_values <- c(seq(10, 100, by = 10),
#                  seq(100, 1000, by = 100)[-1],
#                  seq(1000, 10000, by = 1000)[-1])

nobs_values <- c(300, 5000)

# xi_values <- seq(0, 100, by = 10)

# Simulate
IV_b_poi <- simulate_IV(nsim = nsim, nobs_values = nobs_values,
                        xi_values = xi_values,
                        data_table = data_table,
                        formula = formula, A_formula = A_formula,
                        family = family, type = type)

IV_b_data <- IV_b_poi$sim_data
plot_IV(IV_b_data, causal_parameter = b)





# EX3: binomial ---------------------------------------------------------------


# sim1: Binomial IV for fixed v ---

# Investigate mean
def_bin_X <- defData( varname = "Y", dist = "binomial", link = "logit", 
                      formula = "0.5", variance = 1)
dd <- genData(5000, def_bin_X)
hist(dd$Y)

# Define variables for unperturbed and perturbed data set
# def_bin_X <- defData(varname = "randA", dist = "binary",
#                       formula = 0.5)
# def_bin_X <- defData(def_bin_X, varname = "A", 
#                       formula = "2 * randA - 1")


def_bin_X <- defData(varname = "A", dist = "normal", 
                     formula = 1, variance = 0.25)
def_bin_X <- defData(def_bin_X, varname = "H", dist = "normal",
                     formula = 2, variance = 0.25)
def_bin_X <- defData(def_bin_X, varname = "X", dist = "normal", 
                     formula = "- 3 * H + 2.5 * A", variance = 0.25)
def_bin_X <- defData(def_bin_X, varname = "Y", dist = "binomial", link = "logit", 
                     formula = "1 * X + 2 * H", variance = 1)

def_bin_X_pert <- defData(varname = "A", dist = "normal", 
                          formula = 1, variance = 0.25)
def_bin_X_pert <- defData(def_bin_X_pert, varname = "H", dist = "normal",
                          formula = 2, variance = 0.25)
def_bin_X_pert <- defData(def_bin_X_pert, varname = "v", 
                          formula = 3.5) # set perturbation strength
def_bin_X_pert <- defData(def_bin_X_pert, varname = "X", dist = "normal", 
                          formula = "- 3 * H + v ", variance = 0.25)
def_bin_X_pert <- defData(def_bin_X_pert, varname = "Y", dist = "binomial",
                          link = "logit", 
                          formula = "1 * X + 2 * H", variance = 1)

b <- 1
h <- 2

dd <- genData(300, def_bin_X)
dd_pert <- genData(300, def_bin_X_pert)
par(mfrow = c(2,3))
hist(b*dd$X + h * dd$H, breaks = 100)
hist(exp(b*dd$X + h * dd$H)/(1+exp(b*dd$X + h * dd$H)), breaks = 100)
hist(dd$Y, breaks = 100)
hist(b*dd_pert$X + h * dd_pert$H, breaks = 100)
hist(exp(b*dd_pert$X + h * dd_pert$H)/(1+exp(b*dd_pert$X + h * dd_pert$H)), breaks = 100)
hist(dd_pert$Y, breaks = 100)

# library(Hmisc)
#  hist.data.frame(dd[ , c("Y", "X", "A", "H")])



# sim2: Binomial IV for fixed xi ---

# Initialize for fixed xi
set.seed(16336)

nsim <- 10

xi_values <- c(0, 1, 10000)
v_values <- seq(-10, 10, by = 1)

data_table <- def_bin_X
data_pert_table <- def_bin_X_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- "binomial"
type <- "deviance"

# Simulate
sim_data_bin_X_fixi <- simulate_fixi(nsim = nsim, nobs = 300,
                                     xi_values = xi_values,
                                     v_values = v_values, 
                                     data_table = data_table,
                                     data_pert_table = data_pert_table, 
                                     formula = formula, A_formula = A_formula,
                                     family = family, type = type,
                                     quant_value = 0.9)

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

plot_fixi2(data_bin_X_fixi)

# sim3: IV investigation for binomial ---

# Initialize for fixed v
set.seed(58435)

# nobs_values <- c(seq(10, 100, by = 10),
#                  seq(100, 1000, by = 100)[-1],
#                  seq(1000, 10000, by = 1000)[-1])

nobs_values <- c(300, 5000)


# Simulate
IV_b_bin <- simulate_IV(nsim = nsim, nobs_values = nobs_values,
                        xi_values = xi_values,
                        data_table = data_table,
                        formula = formula, A_formula = A_formula,
                        family = family, type = type)

IV_b_data_bin <- IV_b_bin$sim_data
plot_IV(IV_b_data_bin, causal_parameter = b)

