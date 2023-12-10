num_chains <- 10  # Number of chains
num_steps <- 15000  # Number of steps per chain
seed <- 1  # Seed for random number generation
scale_value_of_proposal = 0.1
# Creating initial states for each chain
initial_states <- matrix(,nrow=1,ncol=num_chains)
for(i in 1:num_chains){
  initial_states[1,i] = 0
}

target <- function(x) {
  alpha <- 4
  beta <- 2
  result <- 0
  for (xi in x) {
    if (is.na(xi) || xi < 0 || xi > 1) {
      return(-Inf)
    }
    result <- result + log(dbeta(xi, alpha, beta))
  }
  return(result)
}

metropolis_hastings <- function(num_chains, initial_states, num_steps, seed, scale_factor_of_proposal, covar_matrix) {
  set.seed(seed)
  final_states <- array(dim = c(nrow(initial_states), num_steps + 1, num_chains))
  accepted_count <- numeric(num_chains)
  
  for (chain_id in 1:num_chains) {
    interim_states <- matrix(nrow = nrow(initial_states), ncol = num_steps + 1)
    interim_states[, 1] <- initial_states[, chain_id]
    proposal_accepted_cnt <- 0
    
    for (i in 1:num_steps) {
      current_state <- interim_states[, i]
      proposed_state <- current_state + scale_factor_of_proposal * covar_matrix %*% rnorm(length(current_state))
      acceptance_ratio <- exp(target(proposed_state) - target(current_state))
      acceptance_ratio <- min(1, acceptance_ratio)
      if (is.nan(acceptance_ratio) || is.infinite(acceptance_ratio)) {
        next
      }else {
        if (runif(1) < acceptance_ratio) {
          interim_states[, i + 1] <- proposed_state
          proposal_accepted_cnt <- proposal_accepted_cnt + 1
        } else {
          interim_states[, i + 1] <- current_state
        }
      }
    }
    final_states[, , chain_id] <- interim_states
    accepted_count[chain_id] <- proposal_accepted_cnt
  }
  list(final_states = final_states, accepted_counts = accepted_count)
}

metro_interface <- function(num_chains, initial_states, num_steps, seed, scale_factor_of_proposal, covar_object) {
  start_time <- Sys.time()
  results <- metropolis_hastings(num_chains, initial_states, num_steps, seed, scale_factor_of_proposal, covar_object)
  end_time <- Sys.time()
  duration <- end_time - start_time
  return(list(results, duration = duration))
}

results <- metro_interface(num_chains = num_chains,
                           initial_states = initial_states,
                           num_steps = num_steps,
                           seed = seed,
                           scale_factor_of_proposal = scale_value_of_proposal,
                           covar_object = matrix(c(0.1))
)