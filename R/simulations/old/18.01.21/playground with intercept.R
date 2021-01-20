# Ausprobieren mit intercept --------------------------------------------------

# Function definition for IV parameter identifiability ------------------------

onerep_IV_mult <- function(rep, nobs_values = 300, data_table, 
                           formula, A_formula, xi_values,
                           family, type, dim = 2) {
  
  data <- as.data.frame(matrix(ncol = 3 + (2 * dim), nrow = 0))
  colnames(data) <- c("rep", "nobs", "xi", "intercept", "int_se", "b", "b_se")
  
  for (k in 1:length(nobs_values)) {
    
    nobs <- nobs_values[k]
    
    # Generate data set with nobs observations
    dd <- genData(nobs, data_table)
    # dd_pert <- genData(nobs, data_pert_table)
    
    xi_len <- length(xi_values)
    
    b <- matrix(nrow = xi_len, ncol = dim)
    b_se <- matrix(nrow = xi_len, ncol = dim)
    
    # logLik_pert <- numeric(xi_len)
    
    for (i in 1:xi_len) {
      
      xi <- xi_values[i]
      
      fit_temp <- glare(formula = formula,
                        A_formula = A_formula,
                        data = dd,
                        xi = xi,
                        family = family,
                        type = type)
      
      b[i, ] <- as.numeric(coef(fit_temp))
      b_se[i, ] <- fit_temp$coef_se
      
      # logLik_pert_indiv <- logLik(fit_temp, newdata = dd_pert, indiv = TRUE)
      # logLik_pert[i] <-
      #   as.numeric(quantile(logLik_pert_indiv, quant_value))
    }
    
    # Return
    data_temp <- data.frame(rep = rep, nobs = nobs,
                            xi = xi_values,
                            intercept = b[, 1],
                            int_se = b_se[, 1],
                            b = b[, 2],
                            b_se = b_se[, 2])
    
    data <- rbind(data, data_temp)
  }
  
  data
}

# Define simulation function
simulate_IV_mult <- function(nsim, nobs_values = 300, xi_values,
                             data_table, formula, A_formula, family, type,
                             dim =2) {
  
  xi_len <- length(xi_values)
  
  states <- array(dim = c(nsim, 626))
  
  b_data <- as.data.frame(matrix(ncol = 3 + (2 * dim), nrow = 0))
  colnames(b_data) <- c("rep", "nobs", "xi", "intercept", "int_se", "b", "b_se")
  
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  for (r in 1:nsim) {
    states[r, ] <- .Random.seed
    
    data_temp <- onerep_IV_mult(rep = r, nobs_values = nobs_values,
                                data_table = data_table, 
                                formula = formula, A_formula = A_formula,
                                xi_values = xi_values,
                                family = family, type = type, dim = dim)
    
    b_data <- rbind(b_data, data_temp)
    
    Sys.sleep(0.1)
    setTxtProgressBar(pb, r)
  }
  close(pb)
  
  list(states = states, sim_data = b_data)
}

# Function definition for fixed xi --------------------------------------------

onerep_fixi_mult <- function(nobs = 300, data_table, data_pert_table, 
                             formula, A_formula, xi_values,
                             family, type, quant_value = 0.9, dim = 2) {
  
  # Generate data with nobs observations
  set.seed(6423) # use the same seeds for each v to use same errors
  dd <- genData(nobs, data_table)
  set.seed(6423)
  dd_pert <- genData(nobs, data_pert_table)
  
  X_pert <- model.matrix(object = formula, data = dd_pert)
  
  # Fit glares
  xi_len <- length(xi_values)
  b <- matrix(nrow = xi_len, ncol = dim)
  b_se <- matrix(nrow = xi_len, ncol = dim)
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
    
    b[j, ] <- as.numeric(coef(fit))
    b_se[j, ] <- fit$coef_se
    
    logLik_pert_indiv[j, ] <- logLik(fit, newdata = dd_pert, indiv = TRUE)
    logLik_pert[j] <- as.numeric(quantile(logLik_pert_indiv[j, ], quant_value))
    
    if(family == "poisson") {
      
      link_inv <- poisson()$linkinv
      vari <- poisson()$variance
      mu <- link_inv(X_pert %*% b[j, ])
      
      r <- mu
      p <- which(dd_pert$Y > 0)
      r[p] <- (dd_pert$Y * log(dd_pert$Y / mu) - (dd_pert$Y - mu))[p]
      
      devres[j, ] <- sign(dd_pert$Y - mu) * sqrt(2 * r)
      
    } else if(family == "binomial") {
      
      link_inv <- binomial()$linkinv
      vari <- binomial()$variance
      mu <- link_inv(X_pert %*% b[j, ])
      
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
      mu <- link_inv(X_pert %*% b[j, ])
      
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
  data_frame <- data.frame(xi = xi_values,
                           b = b,
                           se = b_se,
                           logLik_pert = logLik_pert,
                           deviance = deviance,
                           pearson = pearson)
  
  list(data = data_frame,
       logLik = logLik_pert_indiv,
       devres = devres,
       peares = peares)
}


# Define simulation function
simulate_fixi_mult <- function(nobs = 300, xi_values, v_values, 
                               data_table, data_pert_table, 
                               formula, A_formula, family, type,
                               quant_value = 0.9, dim = 2) {
  
  v_len <- length(v_values)
  
  data_list <- list()
  
  pb <- txtProgressBar(min = 0, max = v_len, style = 3)
  for (s in 1:v_len) {
    data_pert_table_temp <- updateDef(data_pert_table,
                                      changevar = "v",
                                      newformula = v_values[s])
    data_temp <-
      onerep_fixi_mult(nobs = nobs,
                       data_table = data_table,
                       data_pert_table = data_pert_table_temp, 
                       formula = formula, A_formula = A_formula,
                       xi_values = xi_values, family = family, type = type,
                       quant_value = quant_value, dim = dim)
    
    data_list[[s]] <- list(v_value = v_values[s],
                           data = data_temp)
    
    Sys.sleep(0.1)
    setTxtProgressBar(pb, s)
  }
  close(pb)
  
  data_list
}

# Plot IV function ------------------------------------------------------------

plot_IV_mult <- function(data,
                         causal_intercept = 1,
                         causal_parameter = 0.4,
                         xi_big = 10000) {
  
  # Average over simulations
  mean_data <- data %>%
    group_by(xi, nobs) %>%
    summarise(mean_intercept = mean(intercept), se_int = sd(intercept),
              mean_b = mean(b), se_b = sd(b))
  
  # Calculate CI
  upper_int <- mean_data$mean_intercept +
    qt(0.975, 100 - 1) * mean_data$se_int / sqrt(100)
  lower_int <- mean_data$mean_intercept -
    qt(0.975, 100 - 1) * mean_data$se_int / sqrt(100)
  mean_data$upper_int <- upper_int
  mean_data$lower_int <- lower_int
  upper_b <- mean_data$mean_b +
    qt(0.975, 100 - 1) * mean_data$se_b / sqrt(100)
  lower_b <- mean_data$mean_b -
    qt(0.975, 100 - 1) * mean_data$se_b / sqrt(100)
  mean_data$upper_b <- upper_b
  mean_data$lower_b <- lower_b
  
  # Intercept
  gg_data1 <- mean_data
  gg_data1$xi <- as.factor(gg_data1$xi)
  gg1 <- ggplot(data = gg_data1, aes(y = mean_intercept, x = nobs, group = xi, color = xi)) +
    geom_point() +
    geom_line() + 
    geom_errorbar(aes(x = nobs, ymin = lower_int, ymax = upper_int)) +
    
    geom_hline(yintercept = causal_intercept, linetype = "dashed") +
    
    xlab("Number of Observation") +
    ylab("Intercept Estimate")
  
  # Causal Parameter
  gg_data2 <- gg_data1
  gg2 <- ggplot(data = gg_data2, aes(y = mean_b, x = nobs, group = xi, color = xi)) +
    geom_point() +
    geom_line() + 
    geom_errorbar(aes(x = nobs, ymin = lower_b, ymax = upper_b)) +
    
    geom_hline(yintercept = causal_parameter, linetype = "dashed") +
    
    xlab("Number of Observation") +
    ylab("Causal Parameter Estimate")
  
  
  # nobs 100-5000  
  # gg_data_big2 <- gg_data_big[gg_data_big$nobs >= 100 & gg_data_big$nobs <= 5000, ]
  
  
  print(gg1)
  print(gg2)
  
  invisible(mean_data)
}

# Plot function for fixed xi --------------------------------------------------
plot_fixi_mult <- function(sim_data) {
  
  gg_data <- sim_data
  gg_data$xi_val <- as.character(sim_data$xi) 
  
  gg1 <- ggplot(gg_data, aes(y = logLik_pert, x = v, color = xi_val, group = xi_val)) +
    
    geom_line() +
    
    labs(color = expression(xi)) +
    ylab("90th percentile of log-likelihood") +
    xlab("v")
  
  gg_data2 <- gg_data[gg_data$v <= 3 & gg_data$v >= -3, ]
  
  gg2 <- ggplot(gg_data2, aes(y = logLik_pert, x = v, color = xi_val, group = xi_val)) +
    
    geom_line() +
    
    labs(color = expression(xi)) +
    ylab("90th percentile of log-likelihood") +
    xlab("v")
  
  print(gg1)
  print(gg2)
  invisible(gg_data)
}

# Plot function for fixed xi for deviance res ---------------------------------
plot_deviance <- function(sim_data) {
  
  gg_data <- sim_data
  gg_data$xi_val <- as.character(sim_data$xi) 
  
  gg1 <- ggplot(gg_data, aes(y = deviance, x = v, color = xi_val, group = xi_val)) +
    
    geom_line() +
    
    labs(color = expression(xi)) +
    ylab("Root of averaged squared deviance residuals") +
    xlab("v")
  
  gg_data2 <- gg_data[gg_data$v <= 3 & gg_data$v >= -3, ]
  
  gg2 <- ggplot(gg_data2, aes(y = deviance, x = v, color = xi_val, group = xi_val)) +
    
    geom_line() +
    
    labs(color = expression(xi)) +
    ylab("Root of averaged squared deviance residuals") +
    xlab("v")
  
  print(gg1)
  print(gg2)
  invisible(gg_data)
}


# Plot function for fixed xi for pearson res ----------------------------------
plot_pearson <- function(sim_data) {
  
  gg_data <- sim_data
  gg_data$xi_val <- as.character(sim_data$xi) 
  
  gg1 <- ggplot(gg_data, aes(y = pearson, x = v, color = xi_val, group = xi_val)) +
    
    geom_line() +
    
    labs(color = expression(xi)) +
    ylab("Root of averaged squared pearson residuals") +
    xlab("v")
  
  gg_data2 <- gg_data[gg_data$v <= 3 & gg_data$v >= -3, ]
  
  gg2 <- ggplot(gg_data2, aes(y = pearson, x = v, color = xi_val, group = xi_val)) +
    
    geom_line() +
    
    labs(color = expression(xi)) +
    ylab("Root of averaged squared pearson residuals") +
    xlab("v")
  
  print(gg1)
  print(gg2)
  invisible(gg_data)
}


# ------------------------------------------------------------------------------
# POISSON EXAMPLE ---------------------------------
# Define variables for unperturbed and perturbed data set
def_poi_X <- defData(varname = "A", dist = "normal", 
                     formula = 0, variance = 0.25)
def_poi_X <- defData(def_poi_X, varname = "H", dist = "normal",
                     formula = 0, variance = 0.25)
def_poi_X <- defData(def_poi_X, varname = "X", dist = "normal", 
                     formula = "0.8 * H + 1 * A", variance = 0.25)
def_poi_X <- defData(def_poi_X, varname = "Y", dist = "poisson", link = "log", 
                     formula = "2 + 0.4 * X + 0.5 * H", variance = 1)

def_poi_X_pert <- defData(varname = "A", dist = "normal", 
                          formula = 0, variance = 0.25)
def_poi_X_pert <- defData(def_poi_X_pert, varname = "v", 
                          formula = 5) # set perturbation strength
def_poi_X_pert <- defData(def_poi_X_pert, varname = "H", dist = "normal",
                          formula = 0, variance = 0.25)
def_poi_X_pert <- defData(def_poi_X_pert, varname = "X", dist = "normal", 
                          formula = "0.8 * H + v ", variance = 0.25)
def_poi_X_pert <- defData(def_poi_X_pert, varname = "Y", dist = "poisson",
                          formula = "2 + 0.4 * X + 0.5 * H", link = "log")

b <- 0.4
h <- 0.5
c <- 2

dd <- genData(300, def_poi_X)
dd_pert <- genData(300, def_poi_X_pert)
par(mfrow = c(2,2))
hist(c+b*dd$X + h * dd$H, breaks = 100)
hist(dd$Y, breaks = 100)
hist(c+b*dd_pert$X + h * dd_pert$H, breaks = 100)
hist(dd_pert$Y, breaks = 100)


# Initialize for fixed xi
set.seed(7562)

xi_values <- c(0, 1, 50, 10000)
v_values <- seq(0, 20, by = 1)
nobs_values <- c(300,5000)

data_table <- def_poi_X
data_pert_table <- def_poi_X_pert 
formula <- Y ~ X
A_formula <- ~ A - 1
family <- "poisson"
type <- "deviance"

# Simulate

IV_temp <- simulate_IV_mult(nsim = nsim, nobs_values = nobs_values,
                            xi_values = xi_values,
                            data_table = data_table,
                            formula = formula, A_formula = A_formula,
                            family = family, type = type, dim = 2)

IV_temp_data <- IV_temp$sim_data
plot_IV_mult(IV_temp_data, causal_intercept = c, causal_parameter = b)





# sim2: Binomial IV for fixed xi ---

# Initialize for fixed xi
set.seed(13336)

xi_values <- c(0, 1, 10000)
v_values <- seq(-10, 10, by = 1)

# Simulate
sim_data_poi_X_fixi_mult <- simulate_fixi_mult(nobs = 300,
                                               xi_values = xi_values,
                                               v_values = v_values, 
                                               data_table = data_table,
                                               data_pert_table = data_pert_table, 
                                               formula = formula, A_formula = A_formula,
                                               family = family, type = type,
                                               quant_value = 0.9, dim = 2)

# Compute data
v_len <- length(sim_data_poi_X_fixi_mult)

data_poi_X_fixi_mult <- data.frame(matrix(nrow = 0, ncol = 9))
for (i in 1:v_len) {
  
  v <- sim_data_poi_X_fixi_mult[[i]][["v_value"]]
  
  data_temp <- sim_data_poi_X_fixi_mult[[i]][["data"]][["data"]]
  
  data_poi_X_fixi_mult <- rbind(data_poi_X_fixi_mult, data.frame(v = v, data_temp))
}

head(data_poi_X_fixi_mult)
summary(data_poi_X_fixi_mult)

plot_fixi_mult(data_poi_X_fixi_mult)

# Plot pearson and deviance
plot_deviance(data_poi_X_fixi_mult)
plot_pearson(data_poi_X_fixi_mult)


# # Save simulated data list
# path_name <- "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/simulation_study/mult_try/"
# dir.create(paste(path_name, Sys.Date(), sep = ""))
# save(sim_data_bin_X_fixi_mult, file = paste(paste(path_name, Sys.Date(), sep = ""), "/sim2.Rdata", sep =""))






# BINOMIAL EXAMPLE -------------------------------------------------------------

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
                     formula = 0.5, variance = 0.25)
def_bin_X <- defData(def_bin_X, varname = "H", dist = "normal",
                     formula = 0.5, variance = 0.25)
def_bin_X <- defData(def_bin_X, varname = "X", dist = "normal", 
                     formula = "0.5 * H + 0.7 * A", variance = 0.25)
def_bin_X <- defData(def_bin_X, varname = "Y", dist = "binomial", link = "logit", 
                     formula = "0.2 + 0.6 * X + 0.4 * H", variance = 1)

def_bin_X_pert <- defData(varname = "A", dist = "normal", 
                          formula = 0.5, variance = 0.25)
def_bin_X_pert <- defData(def_bin_X_pert, varname = "H", dist = "normal",
                          formula = 0.5, variance = 0.25)
def_bin_X_pert <- defData(def_bin_X_pert, varname = "v", 
                          formula = 2) # set perturbation strength
def_bin_X_pert <- defData(def_bin_X_pert, varname = "X", dist = "normal", 
                          formula = "0.5 * H + v ", variance = 0.25)
def_bin_X_pert <- defData(def_bin_X_pert, varname = "Y", dist = "binomial",
                          link = "logit", 
                          formula = "0.2 + 0.6 * X + 0.4 * H", variance = 1)

b <- 0.6
h <- 0.4
c <- 0.3

dd <- genData(300, def_bin_X)
dd_pert <- genData(300, def_bin_X_pert)
par(mfrow = c(2,3))
hist(c+b*dd$X + h * dd$H, breaks = 100)
hist(exp(c+b*dd$X + h * dd$H)/(1+exp(c+b*dd$X + h * dd$H)), breaks = 100)
hist(dd$Y, breaks = 100)
hist(c+b*dd_pert$X + h * dd_pert$H, breaks = 100)
hist(exp(c+b*dd_pert$X + h * dd_pert$H)/(1+exp(c+b*dd_pert$X + h * dd_pert$H)), breaks = 100)
hist(dd_pert$Y, breaks = 100)

# library(Hmisc)
#  hist.data.frame(dd[ , c("Y", "X", "A", "H")])

# Initialize for fixed v
set.seed(987)
nsim <- 10

# sim2: Binomial IV for fixed xi ---

# Initialize for fixed xi
set.seed(16336)

xi_values <- c(0, 1, 10000)
v_values <- seq(-10, 10, by = 0.1)

data_table <- def_bin_X
data_pert_table <- def_bin_X_pert
formula <- Y ~ X
A_formula <- ~ A - 1
family <- "binomial"
type <- "deviance"

# Simulate
sim_data_bin_X_fixi_mult <- simulate_fixi_mult(nobs = 300,
                                               xi_values = xi_values,
                                               v_values = v_values, 
                                               data_table = data_table,
                                               data_pert_table = data_pert_table, 
                                               formula = formula, A_formula = A_formula,
                                               family = family, type = type,
                                               quant_value = 0.9)

# Compute data
v_len <- length(sim_data_bin_X_fixi_mult)

data_bin_X_fixi_mult <- data.frame(matrix(nrow = 0, ncol = 9))
for (i in 1:v_len) {
  
  v <- sim_data_bin_X_fixi_mult[[i]][["v_value"]]
  
  data_temp <- sim_data_bin_X_fixi_mult[[i]][["data"]][["data"]]
  
  data_bin_X_fixi_mult <- rbind(data_bin_X_fixi_mult, data.frame(v = v, data_temp))
}

head(data_bin_X_fixi_mult)
summary(data_bin_X_fixi_mult)

plot_fixi_mult(data_bin_X_fixi_mult)

# Plot pearson and deviance
plot_deviance(data_bin_X_fixi_mult)
plot_pearson(data_bin_X_fixi_mult)

# sim3: IV investigation for binomial ---

causal_parameter <- b

# Initialize for fixed v
set.seed(58435)
nsim <- 100
# nobs_values <- c(seq(10, 100, by = 10),
#                  seq(100, 1000, by = 100)[-1],
#                  seq(1000, 10000, by = 1000)[-1])

nobs_values <- c(300, 5000)

# xi_values <- seq(0, 100, by = 10)
xi_values <- c(0, 1, 50, 10000)

data_table <- def_bin_X
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- "binomial"
type <- "deviance"

# Simulate
IV_b_bin_mult <- simulate_IV_mult(nsim = nsim, nobs_values = nobs_values,
                        xi_values = xi_values,
                        data_table = data_table,
                        formula = formula, A_formula = A_formula,
                        family = family, type = type)

IV_temp_data <- IV_b_bin_mult$sim_data
plot_IV_mult(IV_temp_data, causal_intercept = c, causal_parameter = b)






# neues theme fuer ggplots
install.packages("ggpubr")
theme_pubr
+theme_set(ggpubr::theme_pubr)


# packages for plot together: patchwork or cowplot
# ggarrange in einem der packages