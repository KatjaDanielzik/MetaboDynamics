\name{data_sim}
\alias{data_sim}
\docType{data}
\title{
A simulated data set
}
\description{
A simulated data set of 98 metabolites. 3 replicate measurements of 4 time points
and at 3 experimental conditions. Metabolites are in 8 dynamics groups per
experimental condition. 4 groups have varying dynamics between conditions.
}
\usage{data("data_sim")}
\format{
  A data frame with 3528 observations on the following 7 variables.
  \describe{
    \item{\code{metabolite}}{column that specifies the metabolite name}
    \item{\code{condition}}{column that specifies the experimental condition}
    \item{\code{time}}{column that specifies the observed time point}
    \item{\code{replicate}}{column that specifies the measurement replicate}
    \item{\code{measurement}}{column that specifies the metabolite concentration, p.e. from peak area/peak intensity}
    \item{\code{log_m}}{log transformation of "Measurement"}
    \item{\code{m_scaled}}{metabolite and Condition specific scaled (mean=0,sd=1) "log_m"}
    \item{\code{KEGG}}{KEGG ID of metabolites}
  }
}
\source{
 \bold{Script used to create simulated data}

# load KEGG database for assignment of metabolite names:
data("metabolite_modules")

# metabolite_db <- metabolite_modules
# Group <- middle_hierarchy

library(dplyr)
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

  base_mean <- runif(1, min = 0.001, max = 1000)

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
data_sim <- dplyr::left_join(data_sim,name_map_HMDB_CAS[,c("Query","KEGG")],join_by("metabolite"=="Query"))


}

\examples{
data(data_sim)
head(data_sim)
}
\keyword{datasets}
