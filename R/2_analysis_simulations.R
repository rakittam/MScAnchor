
# Investigate data ------------------------------------------------------------

# Intervention on X
head(sim_data_bin_X)
summary(sim_data_bin_X)
plot_fixi(sim_data_bin_X)

# Intervention on H
plot_fixi(sim_data_bin_H)

# Intervention on Y
plot_fixi(sim_data_bin_Y)

# Investigate data ------------------------------------------------------------

# Intervention on X, H & Y
head(sim_data_bin_XHY)
head(sim_data_bin_XHY_add)

plot_fivi(sim_data_bin_XHY, sim_data_bin_XHY_add)
