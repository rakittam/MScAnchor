# Run simulation study to produce estimates and states for analysis
#
# Date: 24.11.20     Author: Maic Rakitta
###############################################################################

library(glare)
library(simstudy)

# Defining variables of data sets for one repetition --------------------------

# Data 1
def_1 <- defData(varname = "H", dist = "normal",
                 formula = 0, variance = 1)
def_1 <- defData(def_1, varname = "A", dist = "normal", 
                 formula = 0, variance = 1)
def_1 <- defData(def_1, varname = "X", dist = "normal", 
                 formula = "H + 1.5 * A", variance = 1)
def_1 <- defData(def_1, varname = "Y", dist = "normal", 
                 formula = "X + 2 * H - 3 * A", variance = 1)

# Perturbed data 1
def_pert_1 <- defData(varname = "H", dist = "normal",
                      formula = 0, variance = 1)
def_pert_1 <- defData(def_pert_1, varname = "A", dist = "normal", 
                      formula = 0, variance = 1)
def_pert_1 <- defData(def_pert_1, varname = "X", dist = "normal", 
                      formula = "H - 1.5 * A", variance = 1)
def_pert_1 <- defData(def_pert_1, varname = "Y", dist = "normal", 
                      formula = "X + 2 * H - 1 * A", variance = 1)

rep <- 1
nobs <- 300
data_table <- def_1
data_pert_table <- def_pert_1
formula <- Y~X-1
A_formula <- ~A-1
nest <- 1
xi_values = seq(-1, 2, by = 0.1)
family <- gaussian
type <- "deviance"



# Define function for one repetition ------------------------------------------
# rep is an index of a simulation run
onerep <- function(rep, nobs = 300, data_table, data_pert_table, 
                   formula, A_formula, nest,
                   xi_values = seq(-1, 2, by = 0.1),
                   family, type) {
  
  # Generate data set with nobs observations
  dd <- genData(nobs, data_table)
  dd_pert <- genData(nobs, data_pert_table)
  
  # Fit
  b_matrix <- matrix(nrow = length(xi_values), ncol = nest)
  b_se_matrix <- matrix(nrow = length(xi_values), ncol = nest)
  loglikelihood_pert <- numeric(length(xi_values))
  
  for (i in 1:length(xi_values)) {
    
    xi <- xi_values[i]
    fit_temp <- glare(formula = formula,
                      A_formula = A_formula,
                      data = dd,
                      xi = xi,
                      family = family,
                      type = type)
    
    b_matrix[i, ] <- as.numeric(coef(fit_temp))
    b_se_matrix[i, ] <- fit_temp$coef_se
    loglikelihood_pert[i] <- logLik(fit_temp, newdata = dd_pert)
  }
  
  # Return
  data.frame(rep = rep,
             xi_values = xi_values,
             b = b_matrix,
             b_se = b_se_matrix,
             loglikelihood_pert = loglikelihood_pert)
}

# Make nsim simulation runs --------------------------------------------------
set.seed(3516)
nsim <- 10

nxi <- length(xi_values)

sim_data_1 <- data.frame(matrix(ncol = 5, nrow = nsim * nxi))
colnames(sim_data_1) <- c("rep", "xi", "b", "se", "loglikelihood_pert")
states_1 <- matrix(ncol = 626, nrow = nsim * nxi)

# Run all nsim reps
for (r in 1:nsim) {
  states_1[r, ] <- .Random.seed
  sim_data_1[(nxi * (r - 1) + 1):(nxi * r), ] <-
    onerep(rep = r, nobs = 300, data_table, data_pert_table, 
           formula, A_formula, nest,
           xi_values = seq(-1, 2, by = 0.1),
           family, type)
}

head(sim_data_1)








# Data 2
def_2 <- defData(varname = "H", dist = "normal",
                formula = 0, variance = 1)
def_2 <- defData(def_2, varname = "A", dist = "normal", 
                formula = 0, variance = 1)
def_2 <- defData(def_2, varname = "X", dist = "normal", 
                formula = "H + 0.5 * A", variance = 1)
def_2 <- defData(def_2, varname = "m",
                formula = 5)
def_2 <- defData(def_2, varname = "Y", dist = "binomial", link = "logit", 
                formula = "3 * X + H - 2 * A", variance = 5)

formula <- Y~X-1
A_formula <- ~A-1
xi <- 2
family <- binomial
type <- "deviance"
data_table <- def_2
nobs <- 300


#path_name <- "C:/Users/maicr/Desktop/Github/glare/data sets/simulation 1/data1"
#write.table(dd1, file=paste(path_name,Sys.Date(),sep = "_"))