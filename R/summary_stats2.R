##### Function for calculating and comparing summary statistics #####

summary_stats <- function(boots_file, n_iterations, num_cores){
  
  if (num_cores == 0){ # if number of cores not specified
    num_cores <- detectCores() - 1 # use total cores - 1
  }
  
  stats <- list()
  
  ### Likelihood ###
  like_stats <- analyze_likelihoods(boots_file, n_iterations) # get likelihood stats
  stats$max_like <- like_stats[1]
  stats$sd_like <- like_stats[2]
  stats$range_like <- like_stats[3]
  
  
  ### SNP windows ###
  win_stats <- analyze_snp_windows(num_cores) # get snp window stats
  stats$dist_mean <- win_stats[1]
  stats$dist_sd <- win_stats[2]
  stats$dist_range <- win_stats[3]
  stats$ratio_mean <- win_stats[4]
  stats$ratio_sd <- win_stats[5]
  stats$ratio_range <- win_stats[6]
  
  
  
  return(stats)
  
}
