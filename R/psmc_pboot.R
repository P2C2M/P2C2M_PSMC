##### Function to read in PSMC bootstrap file, get scaffold length distribution, extract parameters, simulate new datasets, and convert to psmcfa format #####


psmc_pboot <- function(boots_file, psmcfa_file, n_iterations, gen, u, s, n_sims, num_cores){ # boots_file = psmc bootstrap output file, psmcfa_file = psmc input file, n_iterations = number of iterations used for psmc analysis, gen = generation time in years, u = mutation rate in subs/site/gen, s = window size, n_sims = number of simulations, num_cores = number of cores to use for simulations
  
  library(parallel)
  
  ### create PSMC class - may not be necessary?
  psmc <- setClass("psmc", slots = list(name="character", truncated="character", theta="numeric", rho="numeric", gen="numeric",
                                        u="numeric", pattern="character", n_free_lambdas="numeric", n_iterations="numeric", 
                                        N0="numeric", lk0="numeric", base_times="numeric", base_nes="numeric", times="numeric", 
                                        nes="numeric", ms_theta="numeric", ms_rho="numeric", ms_times = "numeric", ms_nes = "numeric", 
                                        ms_cmd="character", macs_theta="numeric", macs_rho="numeric", macs_cmd="character"))
  
  
  ### read in bootstraps and extract parameter information
  boots_ms <- bootstrap2ms(boots_file, n_iterations, gen, u, s = s)
  
  
  ### get scaffold lengths
  total_bp <- scaff_length(s = s, psmcfa_file)
  
  
  ### simulate datasets and convert to psmcfa format
  simulate_master(boots_ms, total_bp, s, n_sims, num_cores)
  
}
 
  

