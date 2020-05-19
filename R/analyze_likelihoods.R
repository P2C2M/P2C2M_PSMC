##### Function to set up reading in results #####

#boots_file
#n_iterations

analyze_likelihoods <- function(boots_file, n_iterations){

  emp_likes <- get_likelihoods(boots_file) 
  emp_likes <- emp_likes[1:n_iterations] # get likelihoods of whole empirical analysis. IMPORTANT: THIS ASSUMES THE RUN THAT INCLUDES ALL DATA IS THE FIRST ENTRY IN THE BOOTSTRAP RESULTS FILE!!!!!
  
  sim_file_list <- list.files(pattern = "ms_sim\\d+.psmc$") # get names of simulation results
  sim_likes <- lapply(sim_file_list, function(x) get_likelihoods(x)) # get likelihhod values for all simulated files
  
  
  emp_like_max <- min(emp_likes) # calculate empirical max
  sim_like_maxes <- lapply(sim_likes, function(x) min(x)) # calculate simulated maxes
  max_p <- compare_onetomany(emp_like_max, sim_like_maxes) # calculate p value
  hist_writer("Max_Likelihood", emp_like_max, sim_like_maxes) # write histogram
  
  emp_like_sd <- sd(emp_likes) # calculate empirical standard deviation
  sim_like_sds <- lapply(sim_likes, function(x) sd(x)) # calculate simulated sds
  sd_p <- compare_onetomany(emp_like_sd, sim_like_sds) # calculate p value
  hist_writer("SD_Likelihood", emp_like_sd, sim_like_sds) # write histogram
  
  emp_like_range <- max(emp_likes) - min (emp_likes) # calculate empirical range
  sim_like_ranges <- lapply(sim_likes, function(x) max(x) - min(x)) # calculate simulated range
  range_p <- compare_onetomany(emp_like_range, sim_like_ranges) # calculate p value
  hist_writer("Range_Likelihood", emp_like_range, sim_like_ranges) # write histogram
  

  print(paste0("Empirical Max: ", emp_like_max))
  print(paste0("Sim Max (mean, min, max): ", mean(unlist(sim_like_maxes)), ", ", min(unlist(sim_like_maxes)), ", ", max(unlist(sim_like_maxes))))
  print(paste0("Empirical sd: ", emp_like_sd))
  print(paste0("Sim SD (mean, min, max): ", mean(unlist(sim_like_sds)), ", ", min(unlist(sim_like_sds)), ", ", max(unlist(sim_like_sds))))
  print(paste0("Empirical Range: ", emp_like_range))
  print(paste0("Sim Range (mean, min, max): ", mean(unlist(sim_like_ranges)), ", ", min(unlist(sim_like_ranges)), ", ", max(unlist(sim_like_ranges))))
  
  return(list(max_p, sd_p, range_p))
}



