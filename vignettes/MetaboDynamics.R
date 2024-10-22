## ----style, echo = FALSE, results = 'asis'------------------------------------
BiocStyle::markdown()

## ----eval=FALSE---------------------------------------------------------------
#  if (!require("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  
#  # The following initializes usage of Bioc devel
#  BiocManager::install(version='devel')
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
  geom_line(aes(x = time, y = log_m, col = metabolite, 
                group = interaction(metabolite, replicate))) +
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
data(longitudinalMetabolomics)
# we only use a subsection of the simulated data (1 condition and subsample of
# the whole dataset) for demonstration purposes
samples <- sample(unique(longitudinalMetabolomics$metabolite),10)
longitudinalMetabolomics <- as.data.frame(SummarizedExperiment::colData(longitudinalMetabolomics))
longitudinalMetabolomics <- longitudinalMetabolomics[longitudinalMetabolomics$condition=="A",]
longitudinalMetabolomics <- longitudinalMetabolomics[longitudinalMetabolomics$metabolite%in%samples,]

# fit model
fits_dynamics <- fit_dynamics_model(
  data = longitudinalMetabolomics, scaled_measurement = "m_scaled", time = "time",
  condition = "condition", max_treedepth = 10,
  adapt_delta = 0.95, # default 0.95
  iter = 5000, 
  cores = 1, 
  chains = 2 # only set to 2 for vignette, default = 4
)

