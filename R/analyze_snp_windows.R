##### Function to get snp window summary stats #####


analyze_snp_windows <- function(num_cores){
  
  psmcfa_files <- list.files(path = ".", pattern = "*.psmcfa") # get list of psmcfa files
  
  snp_stats <- mclapply(psmcfa_files, get_snp_windows, mc.cores = num_cores) # get snp stats
  
  emp_stats <- snp_stats[[1]] # get empirical values
  snp_stats <- snp_stats[2:length(snp_stats)] # remove empirical values
  
  dist_means <- lapply(snp_stats, function(x) x[[1]]) # get means
  dist_sds <- lapply(snp_stats, function(x) x[[2]]) # get sds
  dist_ranges <- lapply(snp_stats, function(x) x[[3]]) # get ranges
  
  ratio_means <- lapply(snp_stats, function(x) x[[4]]) # get snp ratio means
  ratio_sds <- lapply(snp_stats, function(x) x[[5]]) # get snp ratio sds
  ratio_ranges <- lapply(snp_stats, function(x) x[[6]]) # get snp ratio ranges
  
  
  ### Calculate p-values
  dist_mean_p <- compare_onetomany(unlist(emp_stats)[1], unlist(dist_means))
  hist_writer("Mean Distance", unlist(emp_stats)[1], unlist(dist_means)) # write histogram
  print(paste0("Empirical Mean Distance: ", unlist(emp_stats)[1]))
  print(paste0("Sim Mean Distance (mean, min, max): ", mean(unlist(dist_means)), ", ", min(unlist(dist_means)), ", ", max(unlist(dist_means))))
  
  dist_sd_p <- compare_onetomany(unlist(emp_stats)[2], unlist(dist_sds))
  hist_writer("SD Distance", unlist(emp_stats)[2], unlist(dist_sds)) # write histogram
  print(paste0("Empirical SD Distance: ", unlist(emp_stats)[2]))
  print(paste0("Sim SD Distance (mean, min, max): ", mean(unlist(dist_sds)), ", ", min(unlist(dist_sds)), ", ", max(unlist(dist_sds))))
  
  dist_range_p <- compare_onetomany(unlist(emp_stats)[3], unlist(dist_ranges))
  hist_writer("Range Distance", unlist(emp_stats)[3], unlist(dist_ranges)) # write histogram
  print(paste0("Empirical Range Distance: ", unlist(emp_stats)[3]))
  print(paste0("Sim Range Distance (mean, min, max): ", mean(unlist(dist_ranges)), ", ", min(unlist(dist_ranges)), ", ", max(unlist(dist_ranges))))
  
  ratio_mean_p <- compare_onetomany(unlist(emp_stats)[4], unlist(ratio_means))
  hist_writer("Mean Percent SNPs", unlist(emp_stats)[4], unlist(ratio_means)) # write histogram
  print(paste0("Empirical Mean Percent SNPs: ", unlist(emp_stats)[4]))
  print(paste0("Sim Mean Percent SNPs (mean, min, max): ", mean(unlist(ratio_means)), ", ", min(unlist(ratio_means)), ", ", max(unlist(ratio_means))))
  
  ratio_sd_p <- compare_onetomany(unlist(emp_stats)[5], unlist(ratio_sds))
  hist_writer("SD Percent SNPs", unlist(emp_stats)[5], unlist(ratio_sds)) # write histogram
  print(paste0("Empirical SD Perecent SNPs: ", unlist(emp_stats)[5]))
  print(paste0("Sim SD Percent SNPs (mean, min, max): ", mean(unlist(ratio_sds)), ", ", min(unlist(ratio_sds)), ", ", max(unlist(ratio_sds))))
  
  ratio_range_p <- compare_onetomany(unlist(emp_stats)[6], unlist(ratio_ranges))
  hist_writer("Range Percent SNPs", unlist(emp_stats)[6], unlist(ratio_ranges)) # write histogram
  print(paste0("Empirical Range Percent SNPs: ", unlist(emp_stats)[6]))
  print(paste0("Sim Range Percent SNPs (mean, min, max): ", mean(unlist(ratio_ranges)), ", ", min(unlist(ratio_ranges)), ", ", max(unlist(ratio_ranges))))
  
  
  snp_sigs <- list(dist_mean_p, dist_sd_p, dist_range_p, ratio_mean_p, ratio_sd_p, ratio_range_p)
  
  return(snp_sigs)
  
  
  
}