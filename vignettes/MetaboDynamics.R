## ----style, echo = FALSE, results = 'asis'------------------------------------
BiocStyle::markdown()

## ----eval=FALSE---------------------------------------------------------------
#  if (!require("BiocManager", quietly = TRUE)) {
#    install.packages("BiocManager")
#  }
#  
#  # The following initializes usage of Bioc devel
#  BiocManager::install(version = "devel")
#  
#  BiocManager::install("MetaboDynamics")

## ----setup--------------------------------------------------------------------
library(MetaboDynamics)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)

## ----fig.wide=TRUE------------------------------------------------------------
data("longitudinalMetabolomics")
# convert to dataframe
longitudinalMetabolomics <- as.data.frame(SummarizedExperiment::colData(longitudinalMetabolomics))
ggplot(longitudinalMetabolomics, aes(x = measurement)) +
  geom_density() +
  theme_bw() +
  facet_grid(cols = vars(time), rows = vars(condition)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ggtitle("raw data", "raw measurements")

## ----fig.wide=TRUE------------------------------------------------------------
# we standardize to a mean of zero and a standard deviation of one of log-transformed data
ggplot(longitudinalMetabolomics, aes(x = log_m)) +
  geom_density() +
  theme_bw() +
  facet_grid(cols = vars(time), rows = vars(condition)) +
  ggtitle("data", "log-transformed values")

## ----fig.wide=TRUE------------------------------------------------------------
ggplot(longitudinalMetabolomics) +
  geom_line(aes(
    x = time, y = log_m, col = metabolite,
    group = interaction(metabolite, replicate)
  )) +
  theme_bw() +
  xlab("timepoint") +
  theme(legend.position = "none") +
  facet_grid(rows = vars(condition)) +
  ggtitle("raw metabolite dynamics", "color=metabolite")

## ----fig.wide=TRUE------------------------------------------------------------
ggplot(longitudinalMetabolomics) +
  geom_line(aes(
    x = time,
    y = m_scaled, col = metabolite,
    group = interaction(metabolite, replicate)
  )) +
  theme_bw() +
  xlab("timepoint") +
  theme(legend.position = "none") +
  facet_grid(rows = vars(condition)) +
  ggtitle("standardized dynamics", "color=metabolite")

## ----fig.wide=TRUE------------------------------------------------------------
# we can hand a SummarizedExperiment object to the function
data("longitudinalMetabolomics")
# we only use a subsection of the simulated data (1 condition and subsample of
# the whole dataset) for demonstration purposes
samples <- c("UMP","Taurine","Succinate","Phosphocreatine","PEP","Malic acid","L-Cystine","CMP","Citramalic acid","2-Aminomuconate")
longitudinalMetabolomics <- as.data.frame(SummarizedExperiment::colData(longitudinalMetabolomics))
longitudinalMetabolomics <- longitudinalMetabolomics[longitudinalMetabolomics$condition == "A", ]
longitudinalMetabolomics <- longitudinalMetabolomics[longitudinalMetabolomics$metabolite %in% samples, ]

# fit model
fits_dynamics <- fit_dynamics_model(
  data = longitudinalMetabolomics, scaled_measurement = "m_scaled", time = "time",
  condition = "condition", max_treedepth = 10,
  adapt_delta = 0.95, # default 0.95
  iter = 5000,
  cores = 1,
  chains = 2 # only set to 2 for vignette, default = 4
)

## -----------------------------------------------------------------------------
# extract diagnostics
diagnostics <- diagnostics_dynamics(
  data = longitudinalMetabolomics,
  iter = 5000, # number of iterations used for model fitting
  # the dynamic model
  scaled_measurement = "m_scaled",
  fits = fits_dynamics,
  chains = 2 # number of chains used for model fitting
)

plot_diagnostics(
  diagnostics = diagnostics[["model_diagnostics"]],
  data = longitudinalMetabolomics
)

# PPCs can be accessed with
plot_PPC(
  posterior = diagnostics[["posterior_A"]], data = longitudinalMetabolomics,
  scaled_measurement = "m_scaled"
)

## ----fig.wide=TRUE------------------------------------------------------------
# #extract estimates
estimates <- estimates_dynamics(
  condition = "condition",
  data = longitudinalMetabolomics, fits = fits_dynamics, samples = 1,
  iter = 5000, # number of iterations used for model fitting
  chains = 2 # number of chains used for model fitting
)

## ----fig.wide=TRUE------------------------------------------------------------
# 1) the differences between two timepoints
plot_estimates(
  estimates = estimates, data = longitudinalMetabolomics,
  dynamics = FALSE
)

## ----fig.wide=TRUE------------------------------------------------------------
# 2) dynamic profiles
plot_estimates(
  estimates = estimates, data = longitudinalMetabolomics,
  delta_t = FALSE
)

## ----fig.wide=TRUE------------------------------------------------------------
# get distances between vectors
dd_A <- dist(
  estimates[["A"]][, c(
    "mu1_mean", "mu2_mean",
    "mu3_mean", "mu4_mean"
  )],
  method = "euclidean"
)
# hierarchical clustering
clust <- hclust(dd_A, method = "ward.D2")
clust_cut <- cutree(clust, k = 8)
# adding cluster ID to estimates
clust_A <- estimates[["A"]][, c(
  "metabolite", "condition", "mu1_mean", "mu2_mean",
  "mu3_mean", "mu4_mean"
)]
clust_A$cluster <- clust_cut
clust_A

## ----fig.wide=TRUE------------------------------------------------------------
data("cluster")
temp <- cluster
temp <- temp %>% pivot_longer(
  cols = c(mu1_mean, mu2_mean, mu3_mean, mu4_mean),
  names_to = "timepoint", values_to = "mu_mean"
)
ggplot(temp, aes(
  x = as.factor(as.numeric(as.factor(timepoint))),
  y = mu_mean, group = metabolite
)) +
  geom_line() +
  xlab("timepoint") +
  ylab("estimated mean concentration") +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(rows = vars(condition), cols = vars(cluster)) +
  ggtitle("clustered dynamics", "panels=cluster ID")

## -----------------------------------------------------------------------------
data("metabolite_modules")
help("metabolite_modules")
head(metabolite_modules)
data("modules_compounds")
help("modules_compounds")
head(modules_compounds)

## ----fig.wide=TRUE------------------------------------------------------------
data("cluster")
ORA <- ORA_hypergeometric(
  background = modules_compounds,
  annotations = metabolite_modules,
  clusters = cluster, tested_column = "middle_hierarchy"
)
plot_ORA(ORA)

ORA_lower <- ORA_hypergeometric(
  background = modules_compounds,
  annotations = metabolite_modules,
  clusters = cluster[cluster$condition == "A", ],
  tested_column = "lower_hierarchy"
)
plot_ORA(ORA_lower)

