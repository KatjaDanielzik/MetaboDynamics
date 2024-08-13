# data("intra")
# # we are taking t=0h for every radiation dose
# # we will duplicate the measurements for 0h for 2Gy and 10Gy
# temp <- intra[intra$time == "0h", ]
# temp$dose <- "2Gy"
# temp2 <- temp
# temp2$dose <- "10Gy"
# intra <- rbind(intra, temp, temp2)
# rm(temp, temp2)
# intra <- intra %>% mutate(log_cpc = log10(cpc))
# # get time as numeric variable
# intra$time <- gsub("h", "", intra$time)
# intra$time <- as.numeric(intra$time)
# intra <- intra %>% arrange(as.numeric(time))
# intra <- intra %>%
#   group_by(dose, metabolite) %>%
#   mutate(log_cpc_stand = ((log_cpc - mean(log_cpc)) / sd(log_cpc)))

# # fit model
# fits_dynamics <- fit_dynamics_model(
#   data = intra, cpc = "log_cpc_stand",
#   condition = "dose", max_treedepth = 14,
#   adapt_delta = 0.999, iter = 4000, cores = 7)

# # extract diagnostics
# diagnostics_dynamics <- extract_diagnostics_dynamics(
#   data = intra, iter = 4000,
#   fits = fits_dynamics
# )
# diagnostics_dynamics <- diagnostics_dynamics["model_diagnostics"]
# # extract estimates
# estimates_dynamics <- extract_estimates_dynamics(
#   data = intra, fits = fits_dynamics,
#   iter = 4000
# )
