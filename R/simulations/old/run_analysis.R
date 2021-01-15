# Run simulation analysis and plots
#
# Date: 14.12.20     Author: Maic Rakitta
###############################################################################

# sim1: Rothenhaeusler Comparison ---------------------------------------------

data_rot <- read.table(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/sim1/sim1_comp_data_2020-12-14.csv")

head(data_rot)
summary(data_rot)

plot_rot(data_rot)

# sim2: Anchor on X normal Rothenhaeusler (IV setting) ------------------------

# Fixed vi
data_rot_X_fivi <- read.table(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/sim2/sim2_fivi_data_2020-12-14.csv")

head(data_rot_X_fivi)
summary(data_rot_X_fivi)

plot_rot_fivi(data_rot_X_fivi)

# Fixed xi
data_rot_X_fixi <- read.table(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/sim2/sim2_fixi_data_2020-12-14.csv")

head(data_rot_X_fixi)
summary(data_rot_X_fixi)

plot_fixi(data_rot_X_fixi)

# sim3: Anchor on X poisson (IV setting) --------------------------------------

# Fixed vi
data_poi_X_fivi <- read.table(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/sim3/sim3_fivi_data_2020-12-14.csv")

head(data_poi_X_fivi)
summary(data_poi_X_fivi)

mean(data_poi_X_fivi$b[data_poi_X_fivi$xi==10000])

plot_fivi_X(data_poi_X_fivi, xi_big = 10000)


# Several quantiles
load("C:/Users/maicr/Desktop/Github/MScAnchor/data sets/sim3/2020-12-24/sim_data_poi_X_fivi_list.Rdata")

data_poi_X_fivi_list <- sim_data_poi_X_fivi_list$sim_data

str(data_poi_X_fivi_list)

# data_poi_X_fivi_list[[1]] # logLiks and coefs storage
# data_poi_X_fivi_list[[1]][[2]] # logLiks of iteration 1
# 
# data_poi_X_fivi_list[[1]][[2]][, 1] # first column is used xi vector

q_values <- seq(0, 1, by = 0.01)
plot_fivi_quantiles(data_poi_X_fivi_list, q_values, xi_big = 10000)


# Fixed xi
data_poi_X_fixi <- read.table(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/sim3/sim3_fixi_data_2020-12-15.csv")

head(data_poi_X_fixi)
summary(data_poi_X_fixi)

plot_fixi(data_poi_X_fixi)

# sim4: Anchor on X, H and Y poisson ------------------------------------------

# Fixed vi
data_poi_XHY_fivi <- read.table(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/sim4/sim4_fivi_data_2020-12-14.csv")

data_poi_XHY_fivi <- read.table(
  "C:/Users/maicr/Desktop/Github/MScAnchor/data sets/sim4/sim4_fivi_data_2020-12-21.csv")

head(data_poi_XHY_fivi)
summary(data_poi_XHY_fivi)

mean(data_poi_XHY_fivi$b[data_poi_XHY_fivi$xi==10000])

plot_fivi_XHY(data_poi_XHY_fivi, xi_big = xi_big)


# Several quantiles
load("C:/Users/maicr/Desktop/Github/MScAnchor/data sets/sim4/2021-01-04/sim_data_poi_XHY_fivi_list.Rdata")

data_poi_XHY_fivi_list <- sim_data_poi_XHY_fivi_list$sim_data

str(data_poi_XHY_fivi_list)

q_values <- seq(0, 1, by = 0.01)
plot_fivi_quantiles(data_poi_XHY_fivi_list, q_values, xi_big = 10000)






















# IV investigation ------------------------------------------------------------

# Using deviance residuals
load("C:/Users/maicr/Desktop/Github/MScAnchor/data sets/sim_IV/2021-01-13/IV_b_deviance.Rdata")
IV_b_data <- IV_b$sim_data

plot_IV(IV_b_data)

# Using pearson residuals
load("C:/Users/maicr/Desktop/Github/MScAnchor/data sets/sim_IV/2021-01-13/IV_b_pearson.Rdata")

IV_b_data <- IV_b$sim_data
plot_IV(IV_b_data)





IV_b_data <- IV_b$sim_data
plot_IV(IV_b_data, causal_parameter = 1)





