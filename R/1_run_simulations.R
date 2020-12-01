# Run simulation study to produce estimates and states for analysis
#
# Date: 24.11.20     Author: Maic Rakitta
###############################################################################

library(glare)
library(simstudy)

# Rothenhaeusler Comparison ---------------------------------------------------

# Define variables for unperturbed and perturbed data set
def_rot <- defData(varname = "randA", 
                   formula = 0.5, dist = "binary")
def_rot <- defData(def_rot, varname = "H", dist = "normal",
                   formula = 0, variance = 1)
def_rot <- defData(def_rot, varname = "A", 
                   formula = "2 * randA - 1")
def_rot <- defData(def_rot, varname = "X", dist = "normal", 
                   formula = "A + H", variance = 1)
def_rot <- defData(def_rot, varname = "Y", dist = "normal", 
                   formula = "X + 2 * H", variance = 1)

def_rot_pert <- defData(varname = "randA", 
                        formula = 0.5, dist = "binary")
def_rot_pert <- defData(def_rot_pert, varname = "H", dist = "normal",
                        formula = 0, variance = 1)
def_rot_pert <- defData(def_rot_pert, varname = "A", 
                        formula = "2 * randA - 1")
def_rot_pert <- defData(def_rot_pert, varname = "X", dist = "normal", 
                        formula = "1.8 + H", variance = 1)
def_rot_pert <- defData(def_rot_pert, varname = "Y", dist = "normal", 
                        formula = "X + 2 * H", variance = 1)

# Define fitting function for one repetition

#for testing:
# rep <- 1
# nobs <- 300
# data_table <- def_rot
# data_pert_table <- def_rot_pert
# formula <- Y~X-1
# A_formula <- ~A-1
# xi_values <- seq(-0.5, 4, by = 0.1)
# family <- gaussian
# type <- "deviance"

onerep_rot <- function(rep, nobs = 300, data_table, data_pert_table, 
                       formula, A_formula, xi_values = seq(-1, 2, by = 0.1),
                       family, type) {
  
  # Generate data set with nobs observations
  dd <- genData(nobs, data_table)
  dd_pert <- genData(nobs, data_pert_table)
  
  # Fit
  b_ar <- matrix(nrow = length(xi_values), ncol = 1)
  b_ar_se <- matrix(nrow = length(xi_values), ncol = 1)
  
  b_glare <- matrix(nrow = length(xi_values), ncol = 1)
  b_glare_se <- matrix(nrow = length(xi_values), ncol = 1)
  
  loglikelihood_pert <- numeric(length(xi_values))
  
  gamma_values <- 2 * xi_values + 1
  
  for (i in 1:length(xi_values)) {
    
    xi <- xi_values[i]
    gamma <-gamma_values[i]
    
    fit_ar <- anchor_regression(formula = formula,
                                A_formula = A_formula,
                                data = dd,
                                gamma = gamma)
    
    fit_glare <- glare(formula = formula,
                       A_formula = A_formula,
                       data = dd,
                       xi = xi,
                       family = family,
                       type = type)
    
    b_ar[i] <- as.numeric(coef(fit_ar))
    b_ar_se[i] <- coef(summary(fit_ar))[, 2]
    
    b_glare[i] <- as.numeric(coef(fit_glare))
    b_glare_se[i] <- fit_glare$coef_se
    
    loglikelihood_pert[i] <- logLik(fit_glare, newdata = dd_pert)
  }
  
  # Return
  data.frame(rep = rep,
             xi_values = xi_values,
             gamma_values = gamma_values,
             b_ar = b_ar,
             b_ar_se = b_ar_se,
             b_glare = b_glare,
             b_glare_se = b_glare_se,
             loglikelihood_pert = loglikelihood_pert)
}

# Make nsim simulation runs
set.seed(3516)
nsim <- 10

xi_values <- seq(-0.5, 4, by = 0.1)

nxi <- length(xi_values)

sim_data_rot <- data.frame(matrix(ncol = 8, nrow = nsim * nxi))
colnames(sim_data_rot) <- c("rep", "xi", "gamma",
                            "b_ar", "se_ar", "b_glare", "se_glare",
                            "loglikelihood_pert")
states_rot <- matrix(ncol = 626, nrow = nsim * nxi)

for (r in 1:nsim) {
  states_rot[r, ] <- .Random.seed
  sim_data_rot[(nxi * (r - 1) + 1):(nxi * r), ] <-
    onerep_rot(rep = r, nobs = 300, def_rot, def_rot_pert, 
               formula = Y~X-1, A_formula = ~A-1,
               xi_values = xi_values,
               family = gaussian, type = "deviance")
}

head(sim_data_rot)
summary(sim_data_rot)


# Function definition for one-dim interventions -------------------------------

onerep <- function(rep, nobs = 300, data_table, data_pert_table, 
                   formula, A_formula, xi,
                   family, type) {
  
  # Generate data set with nobs observations
  dd <- genData(nobs, data_table)
  dd_pert <- genData(nobs, data_pert_table)
  
  # Fit GLM
  fit_glm <- glare(formula = formula,
                   A_formula = A_formula,
                   data = dd,
                   xi = 0,
                   family = family,
                   type = type)
  b_glm <- as.numeric(coef(fit_glm))
  b_glm_se <- fit_glm$coef_se
  loglikelihood_glm_pert <- logLik(fit_glm, newdata = dd_pert)
  
  # Fit glare
  fit_glare <- glare(formula = formula,
                     A_formula = A_formula,
                     data = dd,
                     xi = xi,
                     family = family,
                     type = type)
  b_glare <- as.numeric(coef(fit_glare))
  b_glare_se <- fit_glare$coef_se
  loglikelihood_glare_pert <- logLik(fit_glare, newdata = dd_pert)
  
  # Fit glare with xi big
  fit_big <- glare(formula = formula,
                   A_formula = A_formula,
                   data = dd,
                   xi = 10000,
                   family = family,
                   type = type)
  b_big <- as.numeric(coef(fit_big))
  b_big_se <- fit_big$coef_se
  loglikelihood_big_pert <- logLik(fit_big, newdata = dd_pert)
  
  # Return
  data.frame(rep = rep,
             xi = xi,
             b_glm = b_glm,
             b_glm_se = b_glm_se,
             loglikelihood_glm_pert = loglikelihood_glm_pert,
             b_glare = b_glare,
             b_glare_se = b_glare_se,
             loglikelihood_glare_pert = loglikelihood_glare_pert,
             b_big = b_big,
             b_big_se = b_big_se,
             loglikelihood_big_pert = loglikelihood_big_pert)
}


# Define simulation function
simulate_function <- function(nsim, nobs = 300, xi, v_values, 
                              data_table, data_pert_table, 
                              formula, A_formula, family, type) {
  
  v_len <- length(v_values)
  
  states <- array(dim = c(nsim, v_len, 626))
  sim_data <- data.frame(matrix(nrow = nsim * v_len, ncol = 12))
  colnames(sim_data) <- c("v", "rep", "xi",
                          "b_glm", "se_glm", "loglikelihood_glm_pert",
                          "b_glare", "se_glare", "loglikelihood_glare_pert",
                          "b_big", "se_big", "loglikelihood_big_pert")
  
  for (r in 1:nsim) {
    for (s in 1:v_len) {
      states[r, s, ] <- .Random.seed
      
      data_pert_table_temp <- updateDef(data_pert_table,
                                        changevar = "v", newformula = v_values[s])
      
      sim_data[(v_len * (r - 1) + s), 1] <- v_values[s]
      sim_data[(v_len * (r - 1) + s), 2:12] <-
        onerep(rep = r, nobs = nobs, data_table, data_pert_table_temp, 
               formula = formula, A_formula = A_formula,
               xi = xi,
               family = family, type = type)
    }
  }
  
  list(states = states, sim_data = sim_data)
}

# Anchor on X (IV setting) ----------------------------------------------------

# Define variables for unperturbed and perturbed data set
def_bin_X <- defData(varname = "A", dist = "normal", 
                     formula = 0, variance = 1)
def_bin_X <- defData(def_bin_X, varname = "H", dist = "normal",
                     formula = 0, variance = 1)
def_bin_X <- defData(def_bin_X, varname = "X", dist = "normal", 
                     formula = "2 * H + A", variance = 1)
def_bin_X <- defData(def_bin_X, varname = "m",
                     formula = 5)
def_bin_X <- defData(def_bin_X, varname = "Y", dist = "binomial", link = "logit", 
                     formula = "0.8* X + 2 * H", variance = 5)

def_bin_X_pert <- defData(varname = "A", dist = "normal", 
                          formula = 0, variance = 1)
def_bin_X_pert <- defData(def_bin_X_pert, varname = "v", 
                          formula = 1.5) # set perturbation strength
def_bin_X_pert <- defData(def_bin_X_pert, varname = "H", dist = "normal",
                          formula = 0, variance = 1)
def_bin_X_pert <- defData(def_bin_X_pert, varname = "X", dist = "normal", 
                          formula = "2 * H + v * A", variance = 1)
def_bin_X_pert <- defData(def_bin_X_pert, varname = "m",
                          formula = 5)
def_bin_X_pert <- defData(def_bin_X_pert, varname = "Y", dist = "binomial", link = "logit", 
                          formula = "0.8 * X + 2 * H", variance = 5)

dd <- genData(300, def_bin_X)
dd_pert <- genData(300, def_bin_X_pert)
hist(dd$Y)

# Initialize
set.seed(3516)
nsim <- 10

xi <- 2
v_values <- seq(0, 10, by = 0.5)

data_table <- def_bin_X
data_pert_table <- def_bin_X_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- binomial
type <- "pearson"

# Simulate
simulated_data_X <- simulate_function(nsim, nobs = 300, xi, v_values, 
                                      data_table, data_pert_table, 
                                      formula, A_formula, family, type)

sim_data_bin_X_states <- simulated_data_X$states
sim_data_bin_X <- simulated_data_X$sim_data

head(sim_data_bin_X)
summary(sim_data_bin_X)

# Anchor on H -----------------------------------------------------------------

# Define variables for unperturbed and perturbed data set
def_bin_H <- defData(varname = "A", dist = "normal", 
                     formula = 0, variance = 1)
def_bin_H <- defData(def_bin_H, varname = "H", dist = "normal",
                     formula = "A", variance = 1)
def_bin_H <- defData(def_bin_H, varname = "X", dist = "normal", 
                     formula = "H", variance = 1)
def_bin_H <- defData(def_bin_H, varname = "m",
                     formula = 5)
def_bin_H <- defData(def_bin_H, varname = "Y", dist = "binomial", link = "logit", 
                     formula = "3 * X + H", variance = 5)

def_bin_H_pert <- defData(varname = "A", dist = "normal", 
                          formula = 0, variance = 1)
def_bin_H_pert <- defData(def_bin_H_pert, varname = "v", 
                          formula = 1.5) # set perturbation strength
def_bin_H_pert <- defData(def_bin_H_pert, varname = "H", dist = "normal",
                          formula = "v * A", variance = 1)
def_bin_H_pert <- defData(def_bin_H_pert, varname = "X", dist = "normal", 
                          formula = "H", variance = 1)
def_bin_H_pert <- defData(def_bin_H_pert, varname = "m",
                          formula = 5)
def_bin_H_pert <- defData(def_bin_H_pert, varname = "Y", dist = "binomial", link = "logit", 
                          formula = "3 * X + H", variance = 5)

# Initialize
set.seed(3516)
nsim <- 10

xi <- 5
v_values <- seq(0, 10, by = 1)

data_table <- def_bin_H
data_pert_table <- def_bin_H_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- binomial
type <- "pearson"

# Simulate
simulated_data_H <- simulate_function(nsim, nobs = 300, xi, v_values, 
                                      data_table, data_pert_table, 
                                      formula, A_formula, family, type)

sim_data_bin_H_states <- simulated_data_H$states
sim_data_bin_H <- simulated_data_H$sim_data

head(sim_data_bin_H)
summary(sim_data_bin_H)

# Anchor on Y -----------------------------------------------------------------

# Define variables for unperturbed and perturbed data set
def_bin_Y <- defData(varname = "A", dist = "normal", 
                     formula = 0, variance = 1)
def_bin_Y <- defData(def_bin_Y, varname = "H", dist = "normal",
                     formula = 0, variance = 1)
def_bin_Y <- defData(def_bin_Y, varname = "X", dist = "normal", 
                     formula = "H", variance = 1)
def_bin_Y <- defData(def_bin_Y, varname = "m",
                     formula = 5)
def_bin_Y <- defData(def_bin_Y, varname = "Y", dist = "binomial", link = "logit", 
                     formula = "3 * X + H + A", variance = 5)

def_bin_Y_pert <- defData(varname = "A", dist = "normal", 
                          formula = 0, variance = 1)
def_bin_Y_pert <- defData(def_bin_Y_pert, varname = "v", 
                          formula = 1.5) # set perturbation strength
def_bin_Y_pert <- defData(def_bin_Y_pert, varname = "H", dist = "normal",
                          formula = 0, variance = 1)
def_bin_Y_pert <- defData(def_bin_Y_pert, varname = "X", dist = "normal", 
                          formula = "H", variance = 1)
def_bin_Y_pert <- defData(def_bin_Y_pert, varname = "m",
                          formula = 5)
def_bin_Y_pert <- defData(def_bin_Y_pert, varname = "Y", dist = "binomial", link = "logit", 
                          formula = "3 * X + H + v *  A", variance = 5)

# Initialize
set.seed(3516)
nsim <- 10

xi <- 5
v_values <- seq(0, 10, by = 1)

data_table <- def_bin_Y
data_pert_table <- def_bin_Y_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- binomial
type <- "pearson"

# Simulate
simulated_data_Y <- simulate_function(nsim, nobs = 300, xi, v_values, 
                                      data_table, data_pert_table, 
                                      formula, A_formula, family, type)

sim_data_bin_Y_states <- simulated_data_Y$states
sim_data_bin_Y <- simulated_data_Y$sim_data

head(sim_data_bin_Y)
summary(sim_data_bin_Y)


















# Function definition for greater then one-dim interventions ------------------

onerep_multdim <- function(rep, nobs = 300, data_table, data_pert_table, 
                           formula, A_formula, xi_values,
                           family, type) {
  
  # Generate data set with nobs observations
  dd <- genData(nobs, data_table)
  dd_pert <- genData(nobs, data_pert_table)
  
  xi_len <- length(xi_values)
  
  b_glare <- matrix(nrow = xi_len, ncol = 1)
  b_glare_se <- matrix(nrow = xi_len, ncol = 1)
  
  loglikelihood_glare_pert <- numeric(xi_len)
  
  for (i in 1:xi_len) {
    
    xi <- xi_values[i]
    
    fit_temp <- glare(formula = formula,
                      A_formula = A_formula,
                      data = dd,
                      xi = xi,
                      family = family,
                      type = type)
    b_glare[i] <- as.numeric(coef(fit_temp))
    b_glare_se[i] <- fit_temp$coef_se
    loglikelihood_glare_pert[i] <- logLik(fit_temp, newdata = dd_pert)
  }
  
  # Fit GLM
  fit_glm <- glare(formula = formula,
                   A_formula = A_formula,
                   data = dd,
                   xi = 0,
                   family = family,
                   type = type)
  b_glm <- as.numeric(coef(fit_glm))
  b_glm_se <- fit_glm$coef_se
  loglikelihood_glm_pert <- logLik(fit_glm, newdata = dd_pert)
  
  # Fit glare with xi big
  xi_big <- 10000
  fit_big <- glare(formula = formula,
                   A_formula = A_formula,
                   data = dd,
                   xi = xi_big,
                   family = family,
                   type = type)
  b_big <- as.numeric(coef(fit_big))
  b_big_se <- fit_big$coef_se
  loglikelihood_big_pert <- logLik(fit_big, newdata = dd_pert)
  
  # Return
  sim_data_one <- data.frame(rep = rep,
                             xi_values = xi_values,
                             b_glare = b_glare,
                             b_glare_se = b_glare_se,
                             loglikelihood_glare_pert = loglikelihood_glare_pert)
  sim_data_add <- data.frame(rep = rep,
                             xi_glm = 0,
                             b_glm = b_glm,
                             b_glm_se = b_glm_se,
                             loglikelihood_glm_pert = loglikelihood_glm_pert,
                             xi_big = xi_big,
                             b_big = b_big,
                             b_big_se = b_big_se,
                             loglikelihood_big_pert = loglikelihood_big_pert)
  
  list(sim_data_one = sim_data_one, sim_data_add = sim_data_add)
}


# Define simulation function
simulate_multdim_function <- function(nsim, nobs = 300, xi_values, 
                                      data_table, data_pert_table, 
                                      formula, A_formula, family, type) {
  
  xi_len <- length(xi_values)
  
  states <- array(dim = c(nsim, 626))
  
  sim_data <- data.frame(matrix(nrow = nsim * xi_len, ncol = 5))
  colnames(sim_data) <- c("rep", "xi",
                          "b_glare", "se_glare", "loglikelihood_glare_pert")
  sim_data_add <- data.frame(matrix(nrow = nsim, ncol = 9))
  colnames(sim_data_add) <- c("rep", "xi_glm",
                              "b_glm", "se_glm", "loglikelihood_glm_pert",
                              "xi_big",
                              "b_big", "se_big", "loglikelihood_big_pert")
  
  for (r in 1:nsim) {
    
    states[r, ] <- .Random.seed
    sim_data_temp <- 
      onerep_multdim(rep = r,
                     data_table = data_table, data_pert_table = data_pert_table, 
                     formula = formula, A_formula = A_formula,
                     xi_values = xi_values,
                     family = family, type = type)
    sim_data[(xi_len * (r - 1) + + 1):(xi_len * r),] <- sim_data_temp$sim_data_one
    sim_data_add[r,] <- sim_data_temp$sim_data_add
  }
  
  list(states = states, sim_data = sim_data, sim_data_add = sim_data_add)
}




# Anchor on X, H and Y --------------------------------------------------------

# Define variables for unperturbed and perturbed data set
def_bin_XHY <- defData(varname = "A", dist = "normal", 
                       formula = 0, variance = 1)
def_bin_XHY <- defData(def_bin_XHY, varname = "H", dist = "normal",
                       formula = "A", variance = 1)
def_bin_XHY <- defData(def_bin_XHY, varname = "X", dist = "normal", 
                       formula = "H + A", variance = 1)
def_bin_XHY <- defData(def_bin_XHY, varname = "m",
                       formula = 5)
def_bin_XHY <- defData(def_bin_XHY, varname = "Y", dist = "binomial", link = "logit", 
                       formula = "3 * X + H + A", variance = 5)
# HOW TO INTERVERNE?
def_bin_XHY_pert <- defData(varname = "A", dist = "normal", 
                            formula = 0, variance = 1)
def_bin_XHY_pert <- defData(def_bin_XHY_pert, varname = "H", dist = "normal",
                            formula = "0 * A", variance = 1) # set perturbation
def_bin_XHY_pert <- defData(def_bin_XHY_pert, varname = "X", dist = "normal", 
                            formula = "H - 1 * A", variance = 1) # set perturbation 
def_bin_XHY_pert <- defData(def_bin_XHY_pert, varname = "m",
                            formula = 5)
def_bin_XHY_pert <- defData(def_bin_XHY_pert, varname = "Y",
                            dist = "binomial", link = "logit", 
                            formula = "3 * X + H + 2 * A",
                            variance = 5) # set perturbation

# Initialize
set.seed(3516)
nsim <- 10

xi_values <- seq(0, 5, by = 0.1)

data_table <- def_bin_XHY
data_pert_table <- def_bin_XHY_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- binomial
type <- "pearson"

# Simulate
simulated_data_XHY <- simulate_multdim_function(nsim, xi_values = xi_values, 
                                                data_table = data_table,
                                                data_pert_table = data_pert_table, 
                                                formula = formula,
                                                A_formula = A_formula,
                                                family = family, type = type)

sim_data_bin_XHY_states <- simulated_data_XHY$states
sim_data_bin_XHY <- simulated_data_XHY$sim_data
sim_data_bin_XHY_add <- simulated_data_XHY$sim_data_add

head(sim_data_bin_XHY)
head(sim_data_bin_XHY_add)
summary(sim_data_bin_XHY)






# -----------------------------------------------------------------------------

#path_name <- "C:/Users/maicr/Desktop/Github/glare/data sets/simulation 1/data1"
#write.table(dd1, file=paste(path_name,Sys.Date(),sep = "_"))