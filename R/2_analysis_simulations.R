# Rothenhaeusler Comparison ---------------------------------------------------

# Define variables for unperturbed and perturbed data set
def_rot <- defData(varname = "randA", dist = "binary",
                   formula = 0.5)
def_rot <- defData(def_rot, varname = "H", dist = "normal", variance = 1,
                   formula = 0)
def_rot <- defData(def_rot, varname = "A", 
                   formula = "2 * randA - 1")
def_rot <- defData(def_rot, varname = "X", dist = "normal", variance = 1,
                   formula = "A + H")
def_rot <- defData(def_rot, varname = "Y", dist = "normal", variance = 1,
                   formula = "X + 2 * H")

def_rot_pert <- defData(varname = "randA", dist = "binary",
                        formula = 0.5)
def_rot_pert <- defData(def_rot_pert, varname = "H", dist = "normal",
                        variance = 1,
                        formula = 0)
def_rot_pert <- defData(def_rot_pert, varname = "A", 
                        formula = "2 * randA - 1")
def_rot_pert <- defData(def_rot_pert, varname = "v", 
                        formula = 1.8)
def_rot_pert <- defData(def_rot_pert, varname = "X", dist = "normal", 
                        variance = 1,
                        formula = "v + H")
def_rot_pert <- defData(def_rot_pert, varname = "Y", dist = "normal", 
                        variance = 1,
                        formula = "X + 2 * H")

dd <- genData(300, def_rot)
dd_pert <- genData(300, def_rot_pert)
hist(dd$Y)


# Initialize for Rothenhaeusler comparison
set.seed(1813)

# Simulate
sim_data_rot <- simulate_rot(nsim = 10, nobs = 1000,
                             xi_values = seq(-0.5, 50, by = 0.1), xi_big = 10000,
                             data_table = def_rot, data_pert_table = def_rot_pert, 
                             formula = Y ~ X - 1, A_formula = ~ A - 1,
                             family = gaussian) 

data_rot_states <- sim_data_rot$states
data_rot <- sim_data_rot$sim_data

head(data_rot)
summary(data_rot)

plot_rot(data_rot)

# Anchor on X normal Rothenhaeusler (IV setting) ------------------------------

def_rot <- defData(varname = "randA", dist = "binary",
                   formula = 0.5)
def_rot <- defData(def_rot, varname = "H", dist = "normal", variance = 1,
                   formula = 0)
def_rot <- defData(def_rot, varname = "A", 
                   formula = "2 * randA - 1")
def_rot <- defData(def_rot, varname = "X", dist = "normal", variance = 1,
                   formula = "A + H")
def_rot <- defData(def_rot, varname = "Y", dist = "normal", variance = 1,
                   formula = "X + 2 * H")

def_rot_pert <- defData(varname = "randA", dist = "binary",
                        formula = 0.5)
def_rot_pert <- defData(def_rot_pert, varname = "H", dist = "normal",
                        variance = 1,
                        formula = 0)
def_rot_pert <- defData(def_rot_pert, varname = "A", 
                        formula = "2 * randA - 1")
def_rot_pert <- defData(def_rot_pert, varname = "v", 
                        formula = 1.8)
def_rot_pert <- defData(def_rot_pert, varname = "X", dist = "normal", 
                        variance = 1,
                        formula = "v + H")
def_rot_pert <- defData(def_rot_pert, varname = "Y", dist = "normal", 
                        variance = 1,
                        formula = "X + 2 * H")

dd <- genData(300, def_rot)
dd_pert <- genData(300, def_rot_pert)
hist(dd$Y)

# Initialize for fixed v
set.seed(8653)
nsim <- 2

xi_values <- seq(0, 10, by = 0.1)
xi_values <- xi_values[-1]
xi_big <- 10000

data_table <- def_rot
data_pert_table <- def_rot_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- gaussian

# Simulate
sim_data_rot_X_fivi <- simulate_rot_fivi(nsim = nsim, nobs = 300,
                                         xi_values = xi_values, xi_big = xi_big,
                                         data_table = data_table,
                                         data_pert_table = data_pert_table, 
                                         formula = formula, A_formula = A_formula,
                                         family = family,
                                         quant_value = 0.9) 

data_rot_X_states_fivi <- sim_data_rot_X_fivi$states
data_rot_X_fivi <- sim_data_rot_X_fivi$sim_data

head(data_rot_X_fivi)
summary(data_rot_X_fivi)

plot_rot_fivi(data_rot_X_fivi, xi_big = xi_big)

# Initialize for fixed xi
set.seed(7863)

xi_values <- c(0.2, 0.5, 0.8)
xi_big <- 10000
v_values <- seq(0, 10, by = 1)

data_table <- def_rot
data_pert_table <- def_rot_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- gaussian
type <- "pearson"

# Simulate
sim_data_rot_X_fixi <- simulate_rot_fixi(nobs = 300,
                                         xi_values = xi_values, xi_big = xi_big,
                                         v_values = v_values, 
                                         data_table = data_table,
                                         data_pert_table = data_pert_table, 
                                         formula = formula, A_formula = A_formula,
                                         family = family, type = type,
                                         quant_value = 0.9)

data_rot_X_fixi <- sim_data_rot_X_fixi$sim_data

head(data_rot_X_fixi)
summary(data_rot_X_fixi)

# red = GLM, green = glare, blue = xi_big
plot_rot_fixi(data_rot_X_fixi)








# Anchor on X, H and Y normal continuous anchor -------------------------------

# Define variables for unperturbed and perturbed data set
def_nor_XHY <- defData(varname = "A", dist = "normal", variance = 0.25,
                       formula = 0)
def_nor_XHY <- defData(def_nor_XHY, varname = "H", dist = "normal",
                       variance = 0.25,
                       formula = "-1 + A",)
def_nor_XHY <- defData(def_nor_XHY, varname = "X", dist = "normal",
                       variance = 0.25,
                       formula = "1.5 * H + 1 * A")
def_nor_XHY <- defData(def_nor_XHY, varname = "Y", dist = "normal", 
                       formula = "1.7 * X - 2 * H - A", variance = 0.25)

def_nor_XHY_pert <- defData(varname = "A", dist = "normal", 
                            variance = 0.25, formula = 0)
def_nor_XHY_pert <- defData(def_nor_XHY_pert, varname = "H", dist = "normal",
                            variance = 0.25,
                            formula = "-1 + 3") # set perturbation
def_nor_XHY_pert <- defData(def_nor_XHY_pert, varname = "X", dist = "normal", 
                            variance = 0.25,
                            formula = "1.5 * H - 2.5") # set perturbation 
def_nor_XHY_pert <- defData(def_nor_XHY_pert, varname = "Y", dist = "normal",
                            formula = "1.7 * X - 2 * H - 4.5", variance = 0.25) # set perturbation 

dd <- genData(300, def_nor_XHY)
dd_pert <- genData(300, def_nor_XHY_pert)
hist(dd$Y)

# Initialize
set.seed(68531)
nsim <- 100

xi_values <- seq(0, 100, by = 1)
xi_values <- xi_values[-1]
xi_big <- 10000

data_table <- def_nor_XHY
data_pert_table <- def_nor_XHY_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- gaussian
type <- "pearson"

# Simulate
sim_data_nor_XHY_fivi <- simulate_fivi(nsim = nsim, nobs = 300,
                                       xi_values = xi_values, xi_big = xi_big,
                                       data_table = data_table,
                                       data_pert_table = data_pert_table, 
                                       formula = formula, A_formula = A_formula,
                                       family = family, type = type,
                                       quant_value = 0.9)

data_nor_XHY_states_fivi <- sim_data_nor_XHY_fivi$states
data_nor_XHY_fivi <- sim_data_nor_XHY_fivi$sim_data

head(data_nor_XHY_fivi)
summary(data_nor_XHY_fivi)

plot_fivi(data_nor_XHY_fivi, xi_big = xi_big)












# Anchor on X poisson (IV setting) -------------------------------------------

# Define variables for unperturbed and perturbed data set
def_poi_X <- defData(varname = "A", dist = "normal", 
                     formula = 0, variance = 0.25)
def_poi_X <- defData(def_poi_X, varname = "H", dist = "normal",
                     formula = 0, variance = 0.25)
def_poi_X <- defData(def_poi_X, varname = "X", dist = "normal", 
                     formula = "1.2 * H + A", variance = 0.25)
def_poi_X <- defData(def_poi_X, varname = "Y", dist = "poisson", link = "log", 
                     formula = "0.4 * X + 0.5 * H", variance = 1)

def_poi_X_pert <- defData(varname = "A", dist = "normal", 
                          formula = 0, variance = 0.25)
def_poi_X_pert <- defData(def_poi_X_pert, varname = "v", 
                          formula = 4) # set perturbation strength
def_poi_X_pert <- defData(def_poi_X_pert, varname = "H", dist = "normal",
                          formula = 0, variance = 0.25)
def_poi_X_pert <- defData(def_poi_X_pert, varname = "X", dist = "normal", 
                          formula = "1.3 * H + v", variance = 0.25)
def_poi_X_pert <- defData(def_poi_X_pert, varname = "Y", dist = "poisson", link = "log",
                          formula = "0.4 * X + 0.5 * H", variance = 1)

dd <- genData(300, def_poi_X)
dd_pert <- genData(300, def_poi_X_pert)
hist(dd$Y)
hist(dd_pert$Y)

# Initialize for fixed v
set.seed(14981)
nsim <- 100

xi_values <- seq(0, 10, by = 0.1)
xi_values <- xi_values[-1]
xi_big <- 10000

data_table <- def_poi_X
data_pert_table <- def_poi_X_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- poisson
type <- "pearson"

# nobs <- 300
# rep <- 1
# data_table <- def_poi_X
# data_pert_table <- def_poi_X_pert

# Simulate
sim_data_poi_X_fivi <- simulate_fivi(nsim = nsim, nobs = 300,
                                     xi_values = xi_values, xi_big = xi_big,
                                     data_table = data_table,
                                     data_pert_table = data_pert_table, 
                                     formula = formula, A_formula = A_formula,
                                     family = family, type = type,
                                     quant_value = 0.9) 

data_poi_X_states_fivi <- sim_data_poi_X_fivi$states
data_poi_X_fivi <- sim_data_poi_X_fivi$sim_data

head(data_poi_X_fivi)
summary(data_poi_X_fivi)

plot_fivi(data_poi_X_fivi, xi_big = xi_big)

# Initialize for fixed xi
set.seed(33524)
nsim <- 1

xi <- 1
xi_big <- 10000
v_values <- seq(-10, 10, by = 1)

data_table <- def_poi_X
data_pert_table <- def_poi_X_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- poisson
type <- "pearson"

# Simulate
sim_data_poi_X_fixi <- simulate_fixi(nsim = nsim, nobs = 300,
                                     xi = xi, xi_big = xi_big,
                                     v_values = v_values, 
                                     data_table = data_table,
                                     data_pert_table = data_pert_table, 
                                     formula = formula, A_formula = A_formula,
                                     family = family, type = type,
                                     quant_value = 0.9)

data_poi_X_states_fixi <- sim_data_poi_X_fixi$states
data_poi_X_fixi <- sim_data_poi_X_fixi$sim_data

head(data_poi_X_fixi)
summary(data_poi_X_fixi)

# red = GLM, green = glare, blue = xi_big
plot_fixi(data_poi_X_fixi)

# Anchor on X, H and Y poisson ------------------------------------------------

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
# HOW TO INTERVERNE?
def_poi_XHY_pert <- defData(varname = "A", dist = "normal", 
                            variance = 0.25,
                            formula = 0)
def_poi_XHY_pert <- defData(def_poi_XHY_pert, varname = "H", dist = "normal",
                            variance = 0.25,
                            formula = "0.6 + 0.2") # set perturbation
def_poi_XHY_pert <- defData(def_poi_XHY_pert, varname = "X", dist = "normal", 
                            variance = 0.25,
                            formula = "H - 1") # set perturbation 
def_poi_XHY_pert <- defData(def_poi_XHY_pert, varname = "Y",
                            dist = "poisson", link = "log", 
                            formula = "0.8 * X - H - 2") # set perturbation 

dd <- genData(300, def_poi_XHY)
dd_pert <- genData(300, def_poi_XHY_pert)
hist(dd$Y)

# Initialize
set.seed(32416)
nsim <- 100

xi_values <- seq(0, 10, by = 0.1)
xi_values <- xi_values[-1]
xi_big <- 10000

data_table <- def_poi_XHY
data_pert_table <- def_poi_XHY_pert
formula <- Y ~ X - 1
A_formula <- ~ A - 1
family <- poisson
type <- "pearson"

# Simulate
sim_data_poi_XHY_fivi <- simulate_fivi(nsim = nsim, nobs = 300,
                                       xi_values = xi_values, xi_big = xi_big,
                                       data_table = data_table,
                                       data_pert_table = data_pert_table, 
                                       formula = formula, A_formula = A_formula,
                                       family = family, type = type,
                                       quant_value = 0.9)

data_poi_XHY_states_fivi <- sim_data_poi_XHY_fivi$states
data_poi_XHY_fivi <- sim_data_poi_XHY_fivi$sim_data

head(data_poi_XHY_fivi)
summary(data_poi_XHY_fivi)

plot_fivi(data_poi_XHY_fivi, xi_big = xi_big)










# -----------------------------------------------------------------------------

#path_name <- "C:/Users/maicr/Desktop/Github/glare/data sets/simulation 1/data1"
#write.table(dd1, file=paste(path_name,Sys.Date(),sep = "_"))