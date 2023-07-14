# Load required packages
library(readr)
library(psy)
library(psych)
library(lavaan)
library(tidyverse)
library(GPArotation)
library(rstudioapi)
library(boot)
library(EFA.dimensions)
library(EFAtools)

# Set seed
set.seed(123)

# Set number of bootstrap iterations
boot_iterations <- 1000 # Set to 1000 or more for analysis

# Minimum number of subjects per item for EFA
min_subjects_item <- 20

# Set cutoff for factor loadings
threshold <- 0.3

# Set the file path
path <- "/Users/davidevans/Library/CloudStorage/OneDrive-Personal/My Projects/UBham/Disability attribution/Analyses/Factor analysis/ADLInterferenceAndAt-UniversalDisabilityI_DATA_2023-07-14_0928.csv"

# Import data from csv
data <- read_csv(path)

# Get questionnaire variables (retain all variables except 'record_id')
questionnaire_vars <- names(data)[!(names(data) %in% "record_id")]

# Count the number of questionnaire variables
num_vars <- length(questionnaire_vars)
print(num_vars)

# Remove rows (cases) with more than 2 missing values in questionnaire items
data <- data[rowSums(is.na(data[, questionnaire_vars])) <= 2, ]

# Replace missing values with column mean for rows (cases) with 2 or fewer missing values
data <- mutate(data, across(all_of(questionnaire_vars), ~ifelse(is.na(.), mean(., na.rm = TRUE), .), .names = "{.col}"))

# Count the number of subjects (rows)
num_subjects <- nrow(data)
print(num_subjects)

# Calculate the minimum number of subjects required per dataset (EFA and CFA)
min_subjects <- num_vars * min_subjects_item

# Calculate the split ratio if/when splitting for EFA:CFA (if not splitting, set to 1.0)
split_ratio <- min_subjects / num_subjects
print(paste("The minimum split ratio is", split_ratio))

# Check if split ratio is less than 0.5 (assuming at least 50% of the data for EFA)
if (split_ratio < 0.5) {
  print(paste("The split ratio has been increased from", split_ratio, "to 0.5"))
  split_ratio <- 0.5
} else if (split_ratio > 0.8) {  # too many items for the number of subjects
  print(paste("Too few available subjects for the number of items to spilt the data."))
  split_ratio <- 1.0
}
print(paste("The split ratio has been set to", split_ratio))

# Split the data into two subsets: one for EFA and one for CFA
data_ids <- data$record_id
split_size <- ceiling(length(data_ids) * split_ratio) # use 'ceiling' to round up to an integer
split_ids <- sample(data_ids, size = split_size)
efa_data <- data[data$record_id %in% split_ids, ]
cfa_data <- data[!(data$record_id %in% split_ids), ]

# Function to get number of factors based on Parallel Analysis (PA)
get_nfactors_pa <- function(data, indices) {
  bootstrap_data <- data[indices, ]  # Create the bootstrap sample
  pa <- psych::fa.parallel(bootstrap_data, fm = "minres", plot = FALSE)
  nfactors_pa <- pa$nfact
  return(nfactors_pa)
}

# Function to get number of factors based on Comparative Data (CD)
get_nfactors_cd <- function(data, indices) {
  bootstrap_data <- data[indices, ]  # Create the bootstrap sample
  cd <- EFAtools::CD(bootstrap_data, n_factors_max = NA, N_pop = 10000, N_samples = 500, 
                     alpha = 0.3, use = "pairwise.complete.obs", 
                     cor_method = "pearson", max_iter = 50)
  nfactors_cd <- cd$n_factors
  return(nfactors_cd)
}

# Function to get number of factors based on Minimum Average Partial (MAP) criterion
get_nfactors_map <- function(data, indices) {
  bootstrap_data <- data[indices, ]  # Create the bootstrap sample
  vss_results <- psych::VSS(bootstrap_data, rotate = "none", fm = "minres", plot = FALSE)
  nfactors_map <- which.min(vss_results$map)
  return(nfactors_map)
}

# Function to get number of factors based on Eigenvalues (and store eigenvalues)
get_nfactors_eigen <- function(data, indices) {
  bootstrap_data <- data[indices, ]  # Create the bootstrap sample
  correlation_matrix <- cor(bootstrap_data)
  eigenvalues <- eigen(correlation_matrix)$values
  nfactors_eigen <- sum(eigenvalues > 1)  # Keep eigenvalues > 1
  return(c(nfactors_eigen, eigenvalues))  # Return a vector
}

# Perform bootstrap with each of the four methods
boot_results_pa <- boot(data = efa_data[, questionnaire_vars], statistic = get_nfactors_pa, R = boot_iterations)
boot_results_cd <- boot(data = efa_data[, questionnaire_vars], statistic = get_nfactors_cd, R = boot_iterations)
boot_results_map <- boot(data = efa_data[, questionnaire_vars], statistic = get_nfactors_map, R = boot_iterations)
boot_results_eigen <- boot(data = efa_data[, questionnaire_vars], statistic = get_nfactors_eigen, R = boot_iterations)

# Extract the number of factors and eigenvalues from the boot_results_eigen
boot_nfactors_eigen <- boot_results_eigen$t[, 1]  # The first element is the number of factors
boot_eigenvalues <- boot_results_eigen$t[, -1]  # The remaining elements are eigenvalues

# Calculate the mean and 95% confidence interval for each eigenvalue
mean_eigenvalues <- colMeans(boot_eigenvalues, na.rm = TRUE)
ci_eigenvalues <- apply(boot_eigenvalues, 2, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))

# Create a data frame for plotting
plot_data <- data.frame(
  eigenvalue_index = 1:length(mean_eigenvalues),
  mean_eigenvalue = mean_eigenvalues,
  ci_lower = ci_eigenvalues[1, ],
  ci_upper = ci_eigenvalues[2, ]
)

# Create a scree plot with 95% confidence intervals
df <- data.frame(nfactors = boot_nfactors_eigen)
ggplot(plot_data, aes(x = eigenvalue_index, y = mean_eigenvalue)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  labs(
    title = "Scree Plot with Bootstrap Eigenvalues",
    x = "Eigenvalue Index",
    y = "Mean Eigenvalue"
  ) +
  theme_minimal()

# Function to create histogram
plot_histogram <- function(nfactors, title) {
  df <- data.frame(nfactors = nfactors)
  
  ggplot(df, aes(x=nfactors)) +
    geom_histogram(binwidth=1, fill="skyblue", color="black") +
    labs(title=paste("Histogram of Bootstrapped Results (", title, ")"),
         x="Number of Factors",
         y="Frequency") +
    scale_x_continuous(breaks = seq(floor(min(df$nfactors)), ceiling(max(df$nfactors)), by = 1))  # Ensure x-axis only displays integers
}

# Plot histograms of bootstrapped results
plot_histogram(boot_results_pa$t, "Parallel Analysis")
plot_histogram(boot_results_cd$t, "Comparative Data")
plot_histogram(boot_results_map$t, "Minimum Average Partial")
plot_histogram(boot_nfactors_eigen, "Eigenvalues > 1")

# Get mode (or median) of bootstrap estimates for number of factors
nfactors_final_pa <- as.integer(names(which.max(table(boot_results_pa$t))))  # or median(boot_results_pa$t)
nfactors_final_cd <- as.integer(names(which.max(table(boot_results_cd$t))))  # or median(boot_results_cd$t)
nfactors_final_map <- as.integer(names(which.max(table(boot_results_map$t))))  # or median(boot_results_map$t)
nfactors_final_eigen <- as.integer(names(which.max(table(boot_results_eigen$t))))  # or median(boot_results_eigen$t)

# Print out the results with a descriptive message
print(paste("Based on the bootstrapped PA criterion, the modal number of factors is", nfactors_final_pa))
print(paste("Based on the bootstrapped CD criterion, the modal number of factors is", nfactors_final_cd))
print(paste("Based on the bootstrapped MAP criterion, the modal number of factors that minimise the average squared partial correlation is", nfactors_final_map))
print(paste("Based on the bootstrapped Eigenvalue > 1 criterion, the modal number of factors for Eigenvalues > 1 is", nfactors_final_eigen))

# Make final decision on number of factors by asking for user input
nfactors_input <- readline(prompt="How many factors do you want to retain for the EFA? ")















## Part 2: Run EFA

# Set number of factors for EFA
nfactors <- as.integer(nfactors_input)
print(paste("The number of factors to be retained for the EFA is", nfactors))

# Helper function to perform EFA
run_factor_analysis <- function(data, nfactors) {
  if(nfactors == 1) {
    rotate_method = 'none'
  } else {
    rotate_method = 'oblimin'
  }
  
  fa_result <- tryCatch({
    psych::fa(r = data, nfactors = nfactors, rotate = rotate_method, 
              scores = "regression", fm = "minres")
  }, error = function(e) {
    print(paste("EFA failed:", e$message))
    return(NULL)
  })
  #print(fa_result)
  return(fa_result)
}

# Function to perform Procrustes rotation of loadings from the bootstrap EFA to match loadings of the initial EFA
procrustes_rotate_loadings <- function(bootstrap_efa, initial_efa) {
  procrustes_result <- tryCatch({
    EFA.dimensions::PROCRUSTES(loadings = bootstrap_efa$loadings, target = initial_efa$loadings, type = 'oblique', verbose = FALSE)
  }, error = function(e) {
    print(paste("Procrustes rotation failed:", e$message))
    return(NULL)
  })
  if(!is.null(procrustes_result)){
    return(procrustes_result$loadingsPROC)
  } else {
    return(NULL)
  }
}

# Function to perform EFA on a bootstrap sample, rotate loadings and return them in vector form
efa_func <- function(data, indices, nfactors, initial_efa) {
  bootstrap_data <- data[indices, ]
  
  bootstrap_efa <- run_factor_analysis(bootstrap_data, nfactors)
  
  if (is.null(bootstrap_efa)) return(list(rotated = NULL, unrotated = NULL))
  
  # Remove NA rows from loadings matrix before converting to data frame
  bootstrap_efa$loadings <- bootstrap_efa$loadings[rowSums(is.na(bootstrap_efa$loadings)) != ncol(bootstrap_efa$loadings), ]
  
  unrotated_loadings_df <- as.data.frame(bootstrap_efa$loadings)
  unrotated_loadings_df$item <- colnames(data)
  #print("Iteration loadings after EFA")
  #print(unrotated_loadings_df)
  
  if(nfactors == 1) {
    loadings_df <- as.data.frame(bootstrap_efa$loadings)
    #print("Iteration final loadings")
  } else {
    rotated_loadings <- procrustes_rotate_loadings(bootstrap_efa, initial_efa)
    if (is.null(rotated_loadings)) return(list(rotated = NULL, unrotated = unrotated_loadings_df))
    
    loadings_df <- as.data.frame(rotated_loadings)
    #print("Iteration final loadings after Procrustes rotation")
  }
  loadings_df$item <- colnames(data)
  #print(loadings_df)
  
  return(list(rotated = loadings_df, unrotated = unrotated_loadings_df))
}

# Function to perform bootstrap EFA
run_bootstrap_efa <- function(efa_data, questionnaire_vars, boot_iterations, nfactors, initial_efa) {
  results_list_rotated <- vector("list", boot_iterations)
  results_list_unrotated <- vector("list", boot_iterations)
  
  for (i in 1:boot_iterations) {
    #print(paste("Bootstrap EFA iteration", i))
    
    result <- efa_func(efa_data[, questionnaire_vars], sample(1:nrow(efa_data), replace = TRUE), nfactors, initial_efa)
    
    if (!is.null(result$rotated)) {
      results_list_rotated[[i]] <- result$rotated
    }
    if (!is.null(result$unrotated)) {
      results_list_unrotated[[i]] <- result$unrotated
    }
  }
  return(list(rotated = results_list_rotated, unrotated = results_list_unrotated))
}

# Function returning summary of factor loadings from bootstrap EFA
loadings_summary <- function(bootstrap_results, nfactors) {
  # Create an empty list to store factor loadings for each item and factor
  factor_loadings <- vector("list", nfactors)
  names(factor_loadings) <- paste0("MR", 1:nfactors)
  
  # Loop over the list of bootstrapped factor analyses
  for(i in 1:length(bootstrap_results)) {
    for(j in 1:nfactors) {
      # If the first iteration, initialize a data.frame
      if(i == 1) {
        factor_loadings[[j]] <- data.frame(item = bootstrap_results[[i]]$item, 
                                           loading = bootstrap_results[[i]][,j],
                                           factor = names(factor_loadings)[j])
      } else {
        # If not the first iteration, bind new loadings to the existing data.frame
        factor_loadings[[j]] <- rbind(factor_loadings[[j]], 
                                      data.frame(item = bootstrap_results[[i]]$item, 
                                                 loading = bootstrap_results[[i]][,j],
                                                 factor = names(factor_loadings)[j]))
      }
    }
  }
  
  # Combine all factor loadings into one data frame
  all_loadings <- do.call(rbind, factor_loadings)
  
  # Calculate the mean and 95% confidence intervals
  results <- all_loadings %>%
    group_by(factor, item) %>%
    summarise(mean = mean(loading),
              ci_lower = quantile(loading, 0.025),
              ci_upper = quantile(loading, 0.975)) %>%
    arrange(factor, desc(mean)) %>%
    select(item, factor, mean, ci_lower, ci_upper)  # reorder columns
  
  return(results)
}

# Function returning summary of factor communalities from bootstrap EFA
communalities_summary <- function(loadings_list) {
  # Use do.call() and rbind to combine the list of data frames into one data frame
  communalities_df <- do.call(rbind, lapply(loadings_list, function(loadings) {
    item_names <- loadings$item # Store item names separately
    loadings <- loadings %>% select(-item) # Remove the item column from loadings data
    communalities <- rowSums(loadings^2) # Calculate communalities
    # Create a data frame containing item names and communalities
    data.frame(item = item_names, communality = communalities, uniqueness = 1 - communalities)
  }))
  
  # Summarise communalities and uniquenesses by calculating mean and CI
  summarized_metrics <- communalities_df %>%
    group_by(item) %>%
    summarise(
      communality = paste0(round(mean(communality), 2), " (", round(quantile(communality, 0.025), 2), ", ", round(quantile(communality, 0.975), 2), ")"),
      uniqueness = paste0(round(mean(uniqueness), 2), " (", round(quantile(uniqueness, 0.025), 2), ", ", round(quantile(uniqueness, 0.975), 2), ")"),
      .groups = "drop"
    )
  return(summarized_metrics)
}

format_loadings <- function(loadings_summary) {
  
  # Transform numeric values into character with specified format
  loadings_summary <- loadings_summary %>%
    mutate(value = paste0(round(mean, 2), " (", round(ci_lower, 2), ", ", round(ci_upper, 2), ")")) %>%
    select(item, factor, value)
  
  # Use pivot_wider to spread factor values into separate columns
  formatted_loadings <- loadings_summary %>%
    pivot_wider(names_from = factor, values_from = value) %>%
    arrange(item)
  
  return(formatted_loadings)
}

# Run initial EFA to create factor structure
initial_efa <- run_factor_analysis(efa_data[, questionnaire_vars], nfactors)
if(is.null(initial_efa)){
  print("Initial EFA failed. Check your data.")
} else {
  print("Initial EFA:")
  print(initial_efa)
  
  # Run boostrap EFA
  bootstrap_results <- run_bootstrap_efa(efa_data, questionnaire_vars, boot_iterations, nfactors, initial_efa)
  
  # Print summaries
  loadings_summary_rotated <- format_loadings(loadings_summary(bootstrap_results$rotated, nfactors))
  print("Final loadings summary after Procrustes rotations:")
  print(loadings_summary_rotated)
  
  loadings_summary_unrotated <- format_loadings(loadings_summary(bootstrap_results$unrotated, nfactors))
  print("Loadings summary without Procrustes rotations:")
  print(loadings_summary_unrotated)
  
  communalities_summary_rotated <- communalities_summary(bootstrap_results$rotated)
  print("Communalities summary after Procrustes rotations:")
  print(communalities_summary_rotated)
  
  communalities_summary_unrotated <- communalities_summary(bootstrap_results$unrotated)
  print("Communalities summary without Procrustes rotations:")
  print(communalities_summary_unrotated)
}


