##### Function for simulating datasets with ms and converting them to psmcfa format #####

simulate_msfa <- function(times_dists, nes_dists, theta_dist, rho_dist, scaff, total_bp, s, n_free_lambdas, file_name){
  
  sim_nes <- lapply(nes_dists, function(x) sample(x, 1, replace = TRUE)) # get nes to simulate
  
  sim_times <- list()
  sim_times <- append(sim_times, 0, after = length(sim_times)) # add 0 for if statement
  for (t in seq(1, n_free_lambdas - 1)){

    success <- FALSE
    while (!success){ # while selected time is less than previous time
      time <- sample(unlist(times_dists[t]), 1, replace = TRUE) # sample a time for interval t
      success <- (time > sim_times[t]) # is time greater than previous time?
    }

    sim_times <- append(sim_times, time, after = length(sim_times)) # append to list

  }

  sim_times <- sim_times[2:length(sim_times)] # remove initial 0
  
  
  eN_cmd <- paste(sprintf("-eN %s %s", sim_times, sim_nes), collapse = " ") # get ms population size change commands
  theta_n <- sample(theta_dist, 1, replace = TRUE) # get theta value
  rho_n <- sample(rho_dist, 1, replace = TRUE) # get rho value
  
  
  sim_file <- file(file_name, "w+") # create file for simulation data
  
  cmd_filename <- paste(strsplit(file_name, "\\.")[[1]][1], "cmd", sep = ".") # create file name for ms command file
  cmd_file <- file(cmd_filename, "w+") # create file for logging ms information
  ms_vals <- paste(2, 1, theta_n, rho_n, paste(sprintf("%s,%s", sim_times, sim_nes), collapse = ","), sep = ",") # format ms values to log
  writeLines(ms_vals, con = cmd_file) # log ms values
  
  
  for (z in seq(1:scaff)){ # for each scaffold
    
    bp <- total_bp[z] # get sequence length for scaffold
    
    theta_s <- theta_n * bp # adjust theta for scaffold size
    rho_s <- rho_n * bp # adjust rho for scaffold size
    
    ms_cmd <- paste("./ms", 2, 1, "-t", theta_s, "-r", rho_s, format(bp, scientific = FALSE), eN_cmd, sep = " ") # create ms command line
    
    positions_list <- system(ms_cmd, intern = TRUE)[6] # run ms command and get positions of segsites
    positions_list <- as.numeric(strsplit(strsplit(positions_list, "positions:")[[1]][2], " ")[[1]][-1]) # reformat and convert to numeric
    
    segsites <- length(positions_list)
    
    nbins <- ceiling(bp / s) # get number of bins
    scaff_fa <- rep("T", nbins) # create bin sequence where T = homozygous
    
    if (1 %in% positions_list){ # if the last position is a snp
      pos <- which(positions_list == 1)
      positions_list[pos] <- .999999999 # this prevents incorrect indices below
    }
    
    
    positions_list2 <- as.integer(format(positions_list * bp / s, 100, 10)) # get bins in sequence
    bin_list <- unique(positions_list2) # get unique bins
    bin_list <- bin_list + 1 # add one because as.integer rounds down
    scaff_fa[bin_list] <- "K" # convert those bins to snp bins
    
    
    scaff_line <- paste0(">", z) # create scaffold label line
    scaff_fa <- paste0(scaff_fa, collapse = "")
    scaff_fa <- substring(scaff_fa, seq(1, nchar(scaff_fa), 60), unique(c(seq(60, nchar(scaff_fa), 60), nchar(scaff_fa)))) # break into 60 character lines
    writeLines(scaff_line, con = sim_file) # write scaff line
    writeLines(scaff_fa, con = sim_file) # write bin lines
  }
  
  
  close(sim_file)
  close(cmd_file)
}


