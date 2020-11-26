# Run simulation study to produce estimates and states for analysis
#
# Date: 24.11.20     Author: Maic Rakitta
###############################################################################

library(glare)
library(simstudy)

# Rothenhaeusler Comparison ---------------------------------------------------

# Define variables for unperturbed and perturbed data set
def_rot <- defData(varname = "H", dist = "normal",
                   formula = 0, variance = 1)
def_rot <- defData(def_rot, varname = "randA", 
                   formula = 0.5, dist = "binary")
def_rot <- defData(def_rot, varname = "A", 
                   formula = "2 * randA - 1")
def_rot <- defData(def_rot, varname = "X", dist = "normal", 
                   formula = "A + H", variance = 1)
def_rot <- defData(def_rot, varname = "Y", dist = "normal", 
                   formula = "X + 2 * H", variance = 1)

def_rot_pert <- defData(varname = "H", dist = "normal",
                        formula = 0, variance = 1)
def_rot_pert <- defData(def_rot_pert, varname = "randA", 
                        formula = 0.5, dist = "binary")
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

# Run all nsim reps
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








# Anchor on X (IV setting) ----------------------------------------------------

# Define variables for unperturbed and perturbed data set
def_bin_IV <- defData(varname = "H", dist = "normal",
                      formula = 0, variance = 1)
def_bin_IV <- defData(def_bin_IV, varname = "A", dist = "normal", 
                      formula = 0, variance = 1)
def_bin_IV <- defData(def_bin_IV, varname = "X", dist = "normal", 
                      formula = "H + A", variance = 1)
def_bin_IV <- defData(def_bin_IV, varname = "m",
                      formula = 5)
def_bin_IV <- defData(def_bin_IV, varname = "Y", dist = "binomial", link = "logit", 
                      formula = "3 * X + H", variance = 5)

def_bin_IV_pert <- defData(varname = "H", dist = "normal",
                           formula = 0, variance = 1)
def_bin_IV_pert <- defData(def_bin_IV_pert, varname = "A", dist = "normal", 
                           formula = 0, variance = 1)
def_bin_IV_pert <- defData(def_bin_IV_pert, varname = "v", 
                           formula = 1.5) # set perturbation strength
def_bin_IV_pert <- defData(def_bin_IV_pert, varname = "X", dist = "normal", 
                           formula = "H + v * A", variance = 1)
def_bin_IV_pert <- defData(def_bin_IV_pert, varname = "m",
                           formula = 5)
def_bin_IV_pert <- defData(def_bin_IV_pert, varname = "Y", dist = "binomial", link = "logit", 
                           formula = "3 * X + H", variance = 5)

# Define fitting function for one repetition

# #for testing:
# rep <- 1
# nobs <- 300
# data_table <- def_bin_IV
# data_pert_table <- def_bin_IV_pert
# formula <- Y~X-1
# A_formula <- ~A-1
# xi_values <- 5  or  xi_values <- seq(-0.5, 4, by = 0.1)
# family <- binomial
# type <- "pearson"

onerep_bin_IV <- function(rep, nobs = 300, data_table, data_pert_table, 
                          formula, A_formula, xi_values = seq(-1, 2, by = 0.1),
                          family, type) {
  
  # Generate data set with nobs observations
  dd <- genData(nobs, data_table)
  dd_pert <- genData(nobs, data_pert_table)
  
  # Fit
  b_glm <- matrix(nrow = length(xi_values), ncol = 1)
  b_glm_se <- matrix(nrow = length(xi_values), ncol = 1)
  
  b_glare <- matrix(nrow = length(xi_values), ncol = 1)
  b_glare_se <- matrix(nrow = length(xi_values), ncol = 1)
  
  loglikelihood_glm_pert <- numeric(length(xi_values))
  loglikelihood_glare_pert <- numeric(length(xi_values))
  
  for (i in 1:length(xi_values)) {
    
    xi <- xi_values[i]
    
    fit_glm <- glare(formula = formula,
                     A_formula = A_formula,
                     data = dd,
                     xi = 0,
                     family = family,
                     type = type)
    
    fit_glare <- glare(formula = formula,
                       A_formula = A_formula,
                       data = dd,
                       xi = xi,
                       family = family,
                       type = type)
    
    b_glm[i] <- as.numeric(coef(fit_glm))
    b_glm_se[i] <- fit_glm$coef_se
    
    b_glare[i] <- as.numeric(coef(fit_glare))
    b_glare_se[i] <- fit_glare$coef_se
    
    loglikelihood_glm_pert[i] <- logLik(fit_glm, newdata = dd_pert)
    loglikelihood_glare_pert[i] <- logLik(fit_glare, newdata = dd_pert)
  }
  
  # Return
  data.frame(rep = rep,
             xi_values = xi_values,
             b_glm = b_glm,
             b_glm_se = b_glm_se,
             b_glare = b_glare,
             b_glare_se = b_glare_se,
             loglikelihood_glm_pert = loglikelihood_glm_pert,
             loglikelihood_glare_pert = loglikelihood_glare_pert)
}

# Make nsim simulation runs
set.seed(3516)
nsim <- 10

xi_values <- 5
v_values <- seq(-5, 5, by = 0.5)
v_len <- length(v_values)

sim_data_bin_IV <- data.frame(matrix(ncol = 9, nrow = nsim * v_len))
colnames(sim_data_bin_IV) <- c("v", "rep", "xi",
                               "b_glm", "se_glm", "b_glare", "se_glare",
                               "loglikelihood_glm_pert",
                               "loglikelihood_glare_pert")
states_bin_IV <- array(dim = c(nsim, v_len, 626))
for (r in 1:nsim) {
  for (s in 1:v_len) {
    states_bin_IV[r, s, ] <- .Random.seed
    
    def_bin_IV_pert <- updateDef(def_bin_IV_pert,
                                 changevar = "v", newformula = v_values[s])
    
    sim_data_bin_IV[(v_len * (r - 1) + s), 1] <- v_values[s]
    sim_data_bin_IV[(v_len * (r - 1) + s), 2:9] <-
      onerep_bin_IV(rep = r, nobs = 300, def_bin_IV, def_bin_IV_pert, 
                    formula = Y ~ X - 1, A_formula = ~ A - 1,
                    xi_values = xi_values,
                    family = binomial, type = "pearson")
  }
}

head(sim_data_bin_IV)
summary(sim_data_bin_IV)



# Anchor on ..... ----------------------------------------------------











#path_name <- "C:/Users/maicr/Desktop/Github/glare/data sets/simulation 1/data1"
#write.table(dd1, file=paste(path_name,Sys.Date(),sep = "_"))