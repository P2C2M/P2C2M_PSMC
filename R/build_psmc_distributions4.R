##### Function for preparing simulations #####

simulate_master <- function(boots_ms, total_bp, s, n_sims, num_cores){
  
  scaff = length(total_bp) # get number of scaffolds
  
  boots_times <- list()
  boots_nes <- list()
  
  for (i in seq(1, (boots_ms[[1]]@n_free_lambdas - 1))){
    boots_times[[length(boots_times) + 1]] <- lapply(boots_ms, function(x) x@ms_times[i]) # get each bootstrap time for interval i
    
    boots_nes[[length(boots_nes) + 1]] <- lapply(boots_ms, function(x) x@ms_nes[i]) # get each bootstrap ne for interval i
  }
  
  boots_thetas <- lapply(boots_ms, function(x) x@macs_theta) # get all thetas (correcting for sequence length)
  boots_rhos <- lapply(boots_ms, function(x) x@macs_rho) # get all rhos (correcting for sequence_length)
  
  
  times_dists <- lapply(boots_times, function(x) rnorm(1000, mean(unlist(x)), sd(unlist(x)))) # get normal distributions for ms times
  nes_dists <- lapply(boots_nes, function(x) rnorm(1000, mean(unlist(x)), sd(unlist(x)))) # get normal distributions for ms nes
  
  theta_dist <- rnorm(1000, mean(unlist(boots_thetas)), sd(unlist(boots_thetas))) # get normal distribution for theta
  rho_dist <- rnorm(1000, mean(unlist(boots_rhos)), sd(unlist(boots_rhos))) # get normal distribution for rho
  
  
  file_name <- lapply(seq(1, n_sims, 1), function(x) paste0("ms_sim", x, ".psmcfa")) # create list of file names for simulations
  
  if (num_cores == 0){ # if number of cores not specified
    num_cores <- detectCores() - 1 # use total cores - 1
  }
  
  mclapply(file_name, function(x) simulate_msfa(times_dists, nes_dists, theta_dist, rho_dist, scaff, total_bp, s, boots_ms[[1]]@n_free_lambdas, x), mc.cores = num_cores) # run all simulations/file conversions in parallel
  
  system("cat *.cmd > ms_cmds.csv") # combine ms command files
  system(" rm *.cmd") # get rid of command files
}


