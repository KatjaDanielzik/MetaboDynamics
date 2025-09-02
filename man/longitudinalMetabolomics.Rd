\name{longitudinalMetabolomics}
\alias{longitudinalMetabolomics}
\docType{data}
\title{
A simulated data set of longitudinal concentration tables of metabolites.
}
\description{
A simulated data set of 98 metabolites. 3 replicate measurements of 4 time points
and at 3 experimental conditions. Metabolites are in 8 dynamics groups per
experimental condition. 4 groups have varying dynamics between conditions.
Is represented as a SummarizedExperiment object where concentration tables of 
each experimental condition are stored in assays (raw concentrations in "concentrations",
log-transformed transformations in "log_con" and scaled log-transformed concentrations"
in "scaled_log") and metabolite names, KEGG IDs,
experimental conditions and clustering solutions per experimental condition 
are stored in colData and timepoint specifications in rowData.
(\link[SummarizedExperiment]{SummarizedExperiment}).

}
\usage{data("longitudinalMetabolomics")}
\format{
  A SummarizedExperiment object with concentration tables in assays. RowData
  contains the time point specification. ColData as specified below. 
  \describe{
    \item{\code{condition}}{experimental condition}
    \item{\code{metabolite}}{metabolite name}
    \item{\code{KEGG}}{KEGG ID of metabolites}
    \item{\code{replicate}}{column that specifies the measurement replicate}
    \item{\code{cluster}}{cluster ID that is condition specific for every metabolite}
  }
}
\source{
 \bold{Script used to create simulated data}
# load KEGG database for assignment of metabolite names:
data("metabolite_modules")

# metabolite_db <- metabolite_modules
# Group <- middle_hierarchy

library(dplyr)
library(SummarizedExperiment)
library(tidyverse)
# Parameters (as before)

n_features <- 98

n_groups <- 8 # Number of groups (randomly choose between 6-8)

n_time_points <- 4 # Number of time points

n_replicates <- 3 # Number of replicates for all features and time points

n_conditions <- 3 # Number of experimental conditions

x_varying_groups <- 4 # Number of groups with varying dynamics across conditions

condition_names <- c("A","B","C")

# Probability matrix for assigning metabolites from different database groups to dynamic groups
# For simplicity, we assume equal probability; customize as needed

group_probabilities <- matrix(c(0.8,rep(0.01,7), #amino acid metabolism
                                rep(0.01,7),0.8, #nucleotide metabolism
                                0.1,0.8,0.8,rep(0.1,5), # energy and carbohydrate metabolism
                                runif(5 * length(unique(metabolite_modules$middle_hierarchy)))),
                              nrow = n_groups,
                              ncol = length(unique(metabolite_modules$middle_hierarchy)))

# Generate group dynamics (base trends over time) for each condition

group_dynamics <- list()

# Define the base group dynamics for condition 1

group_dynamics[[1]] <- lapply(1:n_groups, function(g) {
  trend <- rnorm(n_time_points, mean = g * 2, sd = 0.5)
  return(trend)
})

# Define varying dynamics for selected groups across other conditions

varying_groups <- sample(1:n_groups, x_varying_groups, replace = FALSE)

for (cond in 2:n_conditions) {
  group_dynamics[[cond]] <- group_dynamics[[1]]
  for (g in varying_groups) {
    group_dynamics[[cond]][[g]] <- rnorm(n_time_points, mean = g * 2, sd = 1)
  }
}

# Assign each feature to a group

feature_to_group <- sample(1:n_groups, n_features, replace = TRUE)

# Initialize a list to store the simulated data

simulated_data <- list()

# Assign metabolite names to features

available_metabolites <- metabolite_modules # Copy of metabolite database to keep track of unused names

# Simulate data for each feature across all conditions

for (feature in 1:n_features) {
  
  # Get the group for this feature
  
  group <- feature_to_group[feature]
  
  # Determine probability of each metabolite database group for this dynamic group
  
  group_probs <- group_probabilities[group, ]
  
  # Subset the metabolite database for selection based on group probabilities
  
  metabolite_candidates <- available_metabolites %>%
    group_by(middle_hierarchy) %>%
    mutate(Probability = group_probs[match(middle_hierarchy, unique(metabolite_modules$middle_hierarchy))]) %>%
    ungroup() %>%
    filter(metabolite %in% available_metabolites$metabolite) # Ensure the metabolite is still available
  
  # Randomly sample a metabolite based on these probabilities
  
  metabolite_name <- sample(metabolite_candidates$metabolite, 1, prob = metabolite_candidates$Probability)
  
  # Remove this metabolite from available pool
  
  available_metabolites <- available_metabolites[available_metabolites$metabolite != metabolite_name, ]
  
  # Generate a random base mean for this feature between 0.001 and 1000
  
  base_mean <- runif(1, min = 0.1, max = 1e5)
  
  # Generate feature-specific variances for each time point
  
  feature_variances <- runif(n_time_points, min = 0.1, max = 2)
  
  # Store data for each condition
  
  for (cond in 1:n_conditions) {
    trend <- group_dynamics[[cond]][[group]]
    feature_means <- base_mean * trend / max(abs(trend))
    
    feature_data <- data.frame(
      metabolite = metabolite_name, # Assign metabolite name here
      condition = paste0(condition_names[[cond]]),
      time = rep(1:n_time_points, each = n_replicates),
      replicate = rep(1:n_replicates, times = n_time_points)
    )
    
    # Generate the actual data points with strictly positive concentrations
    
    feature_data$measurement <- unlist(lapply(1:n_time_points, function(t) {
      rlnorm(n_replicates, meanlog = log(feature_means[t]), sdlog = feature_variances[t])
    }))
    
    simulated_data[[length(simulated_data) + 1]] <- feature_data
  }
}

rm(base_mean,cond,feature,feature_means,feature_to_group,feature_variances,
   g,group,group_probs,metabolite_name,n_conditions,n_features,n_groups,n_replicates,
   n_time_points,trend,varying_groups,x_varying_groups,available_metabolites,feature_data,group_dynamics,
   group_probabilities,metabolite_candidates)

# Combine all features and conditions into one data frame

simulated_data_df <- do.call(rbind, simulated_data)


simulated_data_df <- simulated_data_df %>%
  group_by(metabolite, condition) %>%
  mutate(
    log_m = log10(measurement),
    m_scaled = (log_m - mean(log_m)) / sd(log_m)
  )

# add KEGG IDs
name_map_HMDB_CAS <- readr::read_csv("name_map_HMDB_CAS.csv")
longitudinalMetabolomics <- dplyr::left_join(simulated_data_df,name_map_HMDB_CAS[,c("Query","KEGG")],join_by("metabolite"=="Query"))

## concentrations
temp <- longitudinalMetabolomics
temp <- temp%>%
  select(condition,metabolite,KEGG,time,measurement,replicate)%>%
  pivot_wider(id_cols=c(condition,metabolite,KEGG,replicate),
              names_from = time,
              values_from = measurement)
concentrations <- temp%>%ungroup%>%select(-c(condition,metabolite,KEGG,replicate))
## transform matrix so that conditions are in columns to facilitate access
## with colData -> se[,se$condition="A"]
concentrations <- t(as.matrix(concentrations))
row.names(concentrations) <- NULL
# prepare log_transformed data
temp <- longitudinalMetabolomics
temp <- temp%>%
  select(condition,metabolite,KEGG,time,log_m,replicate)%>%
  pivot_wider(id_cols=c(condition,metabolite,KEGG,replicate),
              names_from = time,
              values_from = log_m)
log_con <- temp%>%ungroup%>%select(-c(condition,metabolite,KEGG,replicate))
log_con <- t(as.matrix(log_con))
row.names(log_con) <- NULL
# prepare scaled log_transformed data
temp <- longitudinalMetabolomics
temp <- temp%>%
  select(condition,metabolite,KEGG,time,m_scaled,replicate)%>%
  pivot_wider(id_cols=c(condition,metabolite,KEGG,replicate),
              names_from = time,
              values_from = m_scaled)
scaled_log <- temp%>%ungroup%>%select(-c(condition,metabolite,KEGG,replicate))
scaled_log <- t(as.matrix(scaled_log))
row.names(scaled_log) <- NULL

# simulate cell counts
## I am assuming that all cell counts are from a poisson distribution with the same
## location parameter (not the case for real life applications). One parameter is just used to show application of the dynamics model including cell counts

counts <- rbind(rpois(ncol(concentrations),1e5),
                rpois(ncol(concentrations),1e5),
                rpois(ncol(concentrations),1e5),
                rpois(ncol(concentrations),1e5))  # 882 columns, three entries for same metabolite, four time points

# prepare row and colData ####
row_data <- DataFrame(time_points=c("time_point_1","time_point_2","time_point_3",
                                    "time_point_4"))
col_data <- DataFrame(condition=temp$condition,metabolite=temp$metabolite,
                      KEGG=temp$KEGG,replicate=temp$replicate)

se <- SummarizedExperiment(assays=SimpleList(concentrations=concentrations,
                                             log_con=log_con,
                                             scaled_log=scaled_log,
                                             cell_counts = counts),
                           rowData = row_data,
                           colData = col_data)
# set row and colnames ####
rownames(se) <- c("time_point_1","time_point_2","time_point_3",
                  "time_point_4")
colnames(se) <- temp$metabolite

# add metadata ####
metadata(se)[["data origin"]] <- "Simulated data of 98 metabolites with three
concentration observations at four time points and at three different biological
conditions. Script to construct dataset can be seen with ?longitudinalMetabolomics"

data <- se
data <- fit_dynamics_model(
  data = data,
  scaled_measurement = "m_scaled", assay = "scaled_log",
  max_treedepth = 14, adapt_delta = 0.95, iter = 2000, cores = 4, chains = 4
)

data <- estimates_dynamics(
  data = data, iter = 2000,
  chains = 4, condition = "condition"
)
data <- cluster_dynamics(data)

# save shortened analysis in metadata of se
metadata(se)[["estimates_dynamics"]] <- metadata(data)[["estimates_dynamics"]]
metadata(se)[["cluster"]] <- metadata(data)[["cluster"]]

longitudinalMetabolomics <- se
usethis::use_data(longitudinalMetabolomics,overwrite = TRUE)
}

\keyword{datasets}
