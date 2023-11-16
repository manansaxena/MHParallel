library(MHParallel)
library(Rcpp)
library(RcppEigen)

# Define the initial states, number of chains, number of steps, etc.
num_chains <- 4  # Number of chains
num_steps <- 15000  # Number of steps per chain
seed <- 1  # Seed for random number generation
n_cores <- 4  # Number of cores for parallel processing

# Creating initial states for each chain
initial_states <- lapply(1:num_chains, function(x) {
  as.numeric(0)  # Assuming a 1-dimensional state
})


# initial_states <- lapply(1:num_chains, function(x) {
#   numeric(3)  # 3-dimensional state initialized to zeros
# })

# Call the metro_interface function
results <- metro_interface(num_chains = num_chains,
                           initial_states = initial_states,
                           num_steps = num_steps,
                           seed = seed,
                           n_cores = n_cores,
                           scale_factor_of_proposal = 0.05)

# Extracting and analyzing the results
final_states <- results$particles
timings <- results$timer

# Example of analyzing the results: print the first few steps of the first chain
print(final_states[[1]][, 1:10])

samples <- results[[1]][[3]]  # Assuming results is a list of matrices, one for each chain

# Plot the histogram of the samples for the first chain
hist(samples, probability = TRUE, breaks = 50, main = "Metropolis Samples", xlab = "Value")

# Overlay the true Beta distribution
x=seq(0,1,length=1000)
alpha <- 4  # Set the alpha parameter for the Beta distribution
beta <- 2   # Set the beta parameter for the Beta distribution
curve(dbeta(x, alpha, beta), add = TRUE, col = "blue", lwd = 2)
# Add a legend
# legend("topright", legend = "True Beta Distribution", col = "blue", lwd = 2)
# 
