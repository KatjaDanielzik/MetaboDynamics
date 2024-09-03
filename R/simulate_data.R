# Set seed for reproducibility
set.seed(123)

# Parameters
n_features <- 100 # Number of features
n_groups <- sample(6:8, 1) # Number of groups (randomly choose between 6-8)
n_time_points <- 4 # Number of time points
n_replicates <- 3 # Number of replicates for all features and time points
n_conditions <- 3 # Number of experimental conditions
x_varying_groups <- sample(2:4, 1) # Number of groups with varying dynamics across conditions

# Generate group dynamics (base trends over time) for each condition
group_dynamics <- list()

# Define the base group dynamics for condition 1
group_dynamics[[1]] <- lapply(1:n_groups, function(g) {
  # Simulate a general trend for each group over time points
  trend <- rnorm(n_time_points, mean = g * 2, sd = 0.5)
  return(trend)
})

# Define varying dynamics for selected groups across other conditions
varying_groups <- sample(1:n_groups, x_varying_groups, replace = FALSE) # Randomly select groups to vary

for (cond in 2:n_conditions) {
  group_dynamics[[cond]] <- group_dynamics[[1]] # Start with the base dynamics
  # Introduce different dynamics for the selected groups
  for (g in varying_groups) {
    group_dynamics[[cond]][[g]] <- rnorm(n_time_points, mean = g * 2, sd = 1) # Different trend for these groups
  }
}

# Assign each feature to a group
feature_to_group <- sample(1:n_groups, n_features, replace = TRUE)

# Initialize a list to store the simulated data
simulated_data <- list()

# Simulate data for each feature across all conditions
for (feature in 1:n_features) {
  # Get the group for this feature
  group <- feature_to_group[feature]

  # Generate a random base mean for this feature between 0.001 and 1000
  base_mean <- runif(1, min = 0.001, max = 1000)

  # Generate feature-specific variances for each time point
  feature_variances <- runif(n_time_points, min = 0.1, max = 2) # Random variances for each time point

  # Store data for each condition
  for (cond in 1:n_conditions) {
    # Get the trend for this group in this condition
    trend <- group_dynamics[[cond]][[group]]

    # Adjust group trend by base mean to set feature means
    feature_means <- base_mean * trend / max(abs(trend))

    # Simulate measurements for each time point with replicates
    feature_data <- data.frame(
      Feature = paste0("Feature_", feature),
      Condition = paste0("Condition_", cond),
      TimePoint = rep(1:n_time_points, each = n_replicates),
      Replicate = rep(1:n_replicates, times = n_time_points)
    )

    # Generate the actual data points with strictly positive concentrations
    feature_data$Measurement <- unlist(lapply(1:n_time_points, function(t) {
      rlnorm(n_replicates, meanlog = log(feature_means[t]), sdlog = feature_variances[t])
    }))

    # Store the data for this feature and condition
    simulated_data[[length(simulated_data) + 1]] <- feature_data
  }
}

# Combine all features and conditions into one data frame
simulated_data_df <- do.call(rbind, simulated_data)


simulated_data_df <- simulated_data_df %>%
  group_by(Feature, Condition) %>%
  mutate(
    log_m = log10(Measurement),
    m_scaled = (log_m - mean(log_m)) / sd(log_m)
  )

ggplot(simulated_data_df, aes(x = as.factor(TimePoint), y = log_m, group = Feature)) +
  geom_line(aes(col = Feature)) +
  guides(col = "none") +
  facet_wrap(~Condition)

ggplot(simulated_data_df, aes(x = as.factor(TimePoint), y = m_scaled, group = Feature)) +
  geom_line(aes(col = Feature)) +
  guides(col = "none") +
  facet_wrap(~Condition)
