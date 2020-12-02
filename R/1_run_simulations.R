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
  
  logLik_pert <- numeric(length(xi_values))
  
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
    
    logLik_pert[i] <- logLik(fit_glare, newdata = dd_pert)
  }
  
  # Return
  data.frame(rep = rep,
             xi_values = xi_values,
             gamma_values = gamma_values,
             b_ar = b_ar,
             b_ar_se = b_ar_se,
             b_glare = b_glare,
             b_glare_se = b_glare_se,
             logLik_pert = logLik_pert)
}

# Make nsim simulation runs
set.seed(3516)
nsim <- 10

xi_values <- seq(-0.5, 4, by = 0.1)

nxi <- length(xi_values)

sim_data_rot <- data.frame(matrix(ncol = 8, nrow = nsim * nxi))
colnames(sim_data_rot) <- c("rep", "xi", "gamma",
                            "b_ar", "se_ar", "b_glare", "se_glare",
                            "logLik_pert")
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


# Function definition for fixed v simulation ----------------------------------

onerep_fivi <- function(rep, nobs = 300, data_table, data_pert_table, 
                        formula, A_formula, xi_values, xi_big = 10000,
                        family, type, quant_value = 0.9) {
  
  # Generate data set with nobs observations
  dd <- genData(nobs, data_table)
  dd_pert <- genData(nobs, data_pert_table)
  
  xi_len <- length(xi_values)
  
  b_glare <- matrix(nrow = xi_len, ncol = 1)
  b_glare_se <- matrix(nrow = xi_len, ncol = 1)
  
  logLik_glare_pert <- numeric(xi_len)
  
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
    
    logLik_glare_pert_indiv <- logLik(fit_temp, newdata = dd_pert, indiv = TRUE)
    logLik_glare_pert[i] <-
      as.numeric(quantile(logLik_glare_pert_indiv, quant_value))
    
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
  
  logLik_glm_pert_indiv <- logLik(fit_glm, newdata = dd_pert, indiv = TRUE)
  logLik_glm_pert <- as.numeric(quantile(logLik_glm_pert_indiv, quant_value))
  
  
  # Fit glare with xi big
  fit_big <- glare(formula = formula,
                   A_formula = A_formula,
                   data = dd,
                   xi = xi_big,
                   family = family,
                   type = type)
  
  b_big <- as.numeric(coef(fit_big))
  b_big_se <- fit_big$coef_se
  
  logLik_big_pert_indiv <- logLik(fit_big, newdata = dd_pert, indiv = TRUE)
  logLik_big_pert <- as.numeric(quantile(logLik_big_pert_indiv, quant_value))
  
  # Return
  data.frame(rep = rep,
             xi = c(0, xi_big, xi_values),
             b = c(b_glm, b_big, b_glare),
             b_se = c(b_glm_se, b_big_se, b_glare_se),
             logLik_pert = c(logLik_glm_pert,
                             logLik_big_pert,
                             logLik_glare_pert))
  
  # sim_data_one <- data.frame(rep = rep,
  #                            xi_values = xi_values,
  #                            b_glare = b_glare,
  #                            b_glare_se = b_glare_se,
  #                            logLik_glare_pert = logLik_glare_pert)
  # sim_data_add <- data.frame(rep = rep,
  #                            xi_glm = 0,
  #                            b_glm = b_glm,
  #                            b_glm_se = b_glm_se,
  #                            logLik_glm_pert = logLik_glm_pert,
  #                            xi_big = xi_big,
  #                            b_big = b_big,
  #                            b_big_se = b_big_se,
  #                            logLik_big_pert = logLik_big_pert)
  # 
  # list(sim_data_one = sim_data_one, sim_data_add = sim_data_add)
}

# Define simulation function
simulate_fivi <- function(nsim, nobs = 300, xi_values, xi_big = 10000,
                          data_table, data_pert_table, 
                          formula, A_formula, family, type,
                          quant_value = 0.9) {
  
  xi_len <- length(xi_values)
  
  states <- array(dim = c(nsim, 626))
  
  sim_data <- data.frame(matrix(nrow = nsim * (xi_len + 2), ncol = 5))
  colnames(sim_data) <- c("rep", "xi",
                          "b", "se", "logLik_pert")
  
  for (r in 1:nsim) {
    
    states[r, ] <- .Random.seed
    sim_data[((xi_len + 2) * (r - 1) + + 1):((xi_len + 2) * r), ] <- 
      onerep_fivi(rep = r,
                  data_table = data_table, data_pert_table = data_pert_table, 
                  formula = formula, A_formula = A_formula,
                  xi_values = xi_values, xi_big = xi_big,
                  family = family, type = type,
                  quant_value = quant_value)
  }
  
  list(states = states, sim_data = sim_data)
}

# Function definition for fixed xi simulation ---------------------------------

onerep_fixi <- function(rep, nobs = 300, data_table, data_pert_table, 
                        formula, A_formula, xi, xi_big = 10000,
                        family, type, quant_value = 0.9) {
  
  # Generate data with nobs observations
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
  
  logLik_glm_pert_indiv <- logLik(fit_glm, newdata = dd_pert, indiv = TRUE)
  logLik_glm_pert <- as.numeric(quantile(logLik_glm_pert_indiv, quant_value))
  
  # Fit glare
  fit_glare <- glare(formula = formula,
                     A_formula = A_formula,
                     data = dd,
                     xi = xi,
                     family = family,
                     type = type)
  
  b_glare <- as.numeric(coef(fit_glare))
  b_glare_se <- fit_glare$coef_se
  
  logLik_glare_pert_indiv <- logLik(fit_glare, newdata = dd_pert, indiv = TRUE)
  logLik_glare_pert <- as.numeric(quantile(logLik_glare_pert_indiv, quant_value))
  
  # Fit glare with xi big
  fit_big <- glare(formula = formula,
                   A_formula = A_formula,
                   data = dd,
                   xi = xi_big,
                   family = family,
                   type = type)
  
  b_big <- as.numeric(coef(fit_big))
  b_big_se <- fit_big$coef_se
  
  logLik_big_pert_indiv <- logLik(fit_big, newdata = dd_pert, indiv = TRUE)
  logLik_big_pert <- as.numeric(quantile(logLik_big_pert_indiv, quant_value))
  
  # Return
  data.frame(rep = rep,
             xi = xi,
             b_glm = b_glm,
             b_glm_se = b_glm_se,
             logLik_glm_pert = logLik_glm_pert,
             b_glare = b_glare,
             b_glare_se = b_glare_se,
             logLik_glare_pert = logLik_glare_pert,
             b_big = b_big,
             b_big_se = b_big_se,
             logLik_big_pert = logLik_big_pert)
}


# Define simulation function
simulate_fixi <- function(nsim, nobs = 300, xi, xi_big = 10000, v_values, 
                          data_table, data_pert_table, 
                          formula, A_formula, family, type,
                          quant_value = 0.9) {
  
  v_len <- length(v_values)
  
  states <- array(dim = c(nsim, v_len, 626))
  sim_data <- data.frame(matrix(nrow = nsim * v_len, ncol = 12))
  colnames(sim_data) <- c("v", "rep", "xi",
                          "b_glm", "se_glm", "logLik_glm_pert",
                          "b_glare", "se_glare", "logLik_glare_pert",
                          "b_big", "se_big", "logLik_big_pert")
  
  for (r in 1:nsim) {
    for (s in 1:v_len) {
      states[r, s, ] <- .Random.seed
      
      data_pert_table_temp <- updateDef(data_pert_table,
                                        changevar = "v",
                                        newformula = v_values[s])
      
      sim_data[(v_len * (r - 1) + s), 1] <- v_values[s]
      sim_data[(v_len * (r - 1) + s), 2:12] <-
        onerep_fixi(rep = r, nobs = nobs,
                    data_table = data_table, data_pert_table = data_pert_table, 
                    formula = formula, A_formula = A_formula,
                    xi = xi, xi_big = xi_big, family = family, type = type,
                    quant_value = quant_value)
    }
  }
  
  list(states = states, sim_data = sim_data)
}

# Anchor on X (IV setting) ----------------------------------------------------

# Define variables for unperturbed and perturbed data set
def_bin_X <- defData(varname = "A", dist = "normal", 
                     formula = 0, variance = 0.25)
def_bin_X <- defData(def_bin_X, varname = "H", dist = "normal",
                     formula = 0, variance = 0.25)
def_bin_X <- defData(def_bin_X, varname = "X", dist = "normal", 
                     formula = "0.7 * H + A", variance = 0.25)
def_bin_X <- defData(def_bin_X, varname = "Y", dist = "binomial", link = "logit", 
                     formula = "-log(3) - 0.3 * X + H", variance = 1)

def_bin_X_pert <- defData(varname = "A", dist = "normal", 
                          formula = 0, variance = 0.25)
def_bin_X_pert <- defData(def_bin_X_pert, varname = "v", 
                          formula = 1.5) # set perturbation strength
def_bin_X_pert <- defData(def_bin_X_pert, varname = "H", dist = "normal",
                          formula = 0, variance = 0.25)
def_bin_X_pert <- defData(def_bin_X_pert, varname = "X", dist = "normal", 
                          formula = "0.7 * H + v * A", variance = 0.25)
def_bin_X_pert <- defData(def_bin_X_pert, varname = "Y", dist = "binomial", link = "logit", 
                          formula = "-log(3) - 0.3 * X + H", variance = 1)

dd <- genData(300, def_bin_X)
dd_pert <- genData(300, def_bin_X_pert)
hist(dd$Y)

# Initialize for fixed v
set.seed(4981)
nsim <- 10

xi_values <- seq(0, 10, by = 0.1)
xi_values <- xi_values[-1]
xi_big <- 10000

data_table <- def_bin_X
data_pert_table <- def_bin_X_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- binomial
type <- "pearson"

# nobs <- 300
# rep <- 1
# data_table <- def_bin_X
# data_pert_table <- def_bin_X_pert

# Simulate
data_X_fivi <- simulate_fivi(nsim, nobs = 300, xi_values, xi_big,
                             data_table, data_pert_table, 
                             formula, A_formula, family, type,
                             quant_value = 0.9) 

data_bin_X_states_fivi <- data_X_fivi$states
data_bin_X_fivi <- data_X_fivi$sim_data

head(data_bin_X_fivi)
summary(data_bin_X_fivi)

plot_fivi(data_bin_X_fivi, xi_big = xi_big)

# Initialize for fixed xi
set.seed(3316)
nsim <- 10

xi <- 2
xi_big <- 10000
v_values <- seq(0, 100, by = 1)

data_table <- def_bin_X
data_pert_table <- def_bin_X_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- binomial
type <- "pearson"

# Simulate
data_X_fixi <- simulate_fixi(nsim = nsim, nobs = 300,
                             xi = xi, xi_big = xi_big, v_values = v_values, 
                             data_table = data_table,
                             data_pert_table = data_pert_table, 
                             formula = formula, A_formula = A_formula,
                             family = family, type = type,
                             quant_value = 0.9)

data_bin_X_states_fixi <- data_X_fixi$states
data_bin_X_fixi <- data_X_fixi$sim_data

head(data_bin_X_fixi)
summary(data_bin_X_fixi)

# red = GLM, green = glare, blue = xi_big
plot_fixi(data_bin_X_fixi)






# Anchor on X, H and Y --------------------------------------------------------

# Define variables for unperturbed and perturbed data set
def_bin_XHY <- defData(varname = "A", dist = "normal", 
                       formula = 0, variance = 0.25)
def_bin_XHY <- defData(def_bin_XHY, varname = "H", dist = "normal",
                       formula = "0.6 + A", variance = 0.25)
def_bin_XHY <- defData(def_bin_XHY, varname = "X", dist = "normal", 
                       formula = "H + A", variance = 0.25)
def_bin_XHY <- defData(def_bin_XHY, varname = "Y", dist = "binomial", link = "logit", 
                       formula = "3 * X + H + A", variance = 1)
# HOW TO INTERVERNE?
def_bin_XHY_pert <- defData(varname = "A", dist = "normal", 
                            formula = 0, variance = 0.25)
def_bin_XHY_pert <- defData(def_bin_XHY_pert, varname = "H", dist = "normal",
                            formula = "0.6 + 0 * A", variance = 0.25) # set perturbation
def_bin_XHY_pert <- defData(def_bin_XHY_pert, varname = "X", dist = "normal", 
                            formula = "H - 1 * A", variance = 0.25) # set perturbation 
def_bin_XHY_pert <- defData(def_bin_XHY_pert, varname = "Y",
                            dist = "binomial", link = "logit", 
                            formula = "3 * X + H + 2 * A", # set perturbation
                            variance = 1) 

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
data_XHY <- simulate_fivi(nsim, xi_values = xi_values, 
                          data_table = data_table,
                          data_pert_table = data_pert_table, 
                          formula = formula,
                          A_formula = A_formula,
                          family = family, type = type)

sim_data_bin_XHY_states <- data_XHY$states
sim_data_bin_XHY <- data_XHY$sim_data
sim_data_bin_XHY_add <- data_XHY$sim_data_add

head(sim_data_bin_XHY)
head(sim_data_bin_XHY_add)
summary(sim_data_bin_XHY)






# -----------------------------------------------------------------------------

#path_name <- "C:/Users/maicr/Desktop/Github/glare/data sets/simulation 1/data1"
#write.table(dd1, file=paste(path_name,Sys.Date(),sep = "_"))