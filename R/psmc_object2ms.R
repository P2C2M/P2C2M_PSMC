##### Function to convert PSMC history to ms/MACS command line #####

psmc_object2ms <- function(psmc_object, n = 2, L = 30000000, s, R = 10, r = 1, print_cmd){ # psmc_object = psmc object, n = number of chromosomes (individuals) to simulate, L = length of each chromosome, s = skip (window size?) used in PSMC run, R = recombination rate in hotspots are R times larger, r = replicates, print_cmd = whether to print the commands to the console 
  
  N <- psmc_object@theta / s / (4 * psmc_object@u)
  
  N0 <- psmc_object@base_nes[1] * N
  
  psmc_object@ms_theta <- 4 * N0 * psmc_object@u * L
  psmc_object@ms_rho <- psmc_object@ms_theta * (psmc_object@rho / psmc_object@theta)
  psmc_object@macs_theta <- psmc_object@ms_theta / L
  psmc_object@macs_rho <- psmc_object@ms_rho / L
  
  psmc_object@ms_times <- (psmc_object@base_times / psmc_object@base_nes[1] / 2)[2:psmc_object@n_free_lambdas]
  
  psmc_object@ms_nes <- (psmc_object@base_nes / psmc_object@base_nes[1])[2:psmc_object@n_free_lambdas]
  
  eN_cmd <- paste(sprintf("%s 1 %s", psmc_object@ms_times, psmc_object@ms_nes), collapse = " ")
  
  psmc_object@ms_cmd <- paste("msHot-lite", n, r, "-t", psmc_object@ms_theta, "-r", psmc_object@ms_rho, L, "-l", eN_cmd, sep = " ")
  
  psmc_object@macs_cmd <- paste("macs", n, L, "-t", psmc_object@macs_theta, "-r", psmc_object@macs_rho, eN_cmd, sep = " ")
  
  if (print_cmd == TRUE){
    print(psmc_object@ms_cmd)
    print(psmc_object@macs_cmd)
  }
    
  return(psmc_object)
}
