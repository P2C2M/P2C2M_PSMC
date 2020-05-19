##### Function for converting psmc output file to history file #####


### Testing ###
# file_name = "French_upper.psmc" # name of psmc output file
# n_iterations = 20
# u = 1.1e-08

psmc2history <- function(file_name, n_iterations, gen, u, from_file, save_output){ # file_name = name of PSMC output file or list of lines from output, n_iterations = which iteration of the PSMC to pull values from, gen = generation time (years), u = mutation rate (subs/site/gen), from_file = whether PSMC information should be read from a file, save_output = whether to save history to a file. Note: don't use save_output with from_file = FALSE
  
  if (from_file == TRUE){  
    psmc_file <- file(file_name, "r+")
    psmc_lines <- suppressWarnings(readLines(psmc_file)) # read in psmc ouput file
    close(psmc_file)
    psmc_obj <- psmc(name = file_name, gen = gen, u = u, n_iterations = n_iterations) # create psmc class object
    psmc_obj@truncated <- unlist(strsplit(psmc_obj@name, ".psmc")[[1]])
    
  } else{
    psmc_lines <- file_name
    psmc_obj <- psmc(gen = gen, u = u, n_iterations = n_iterations) # create psmc class object
  }
  
  psmc_meta <- psmc_lines[which(lapply(psmc_lines, function(x) strsplit(x, "\t")[[1]][1]) == "MM")][1:5] # get metadata lines
  psmc_meta <- lapply(psmc_meta, function(x) strsplit(x, "\t")[[1]][2]) # remove "MM"
  
  psmc_obj@pattern <- unlist(strsplit(unlist(strsplit(psmc_meta[[2]][1], ","))[1], ":")[[1]][2]) # get pattern of atomic intervals
  psmc_obj@n_free_lambdas <- as.numeric(unlist(strsplit(unlist(strsplit(psmc_meta[[2]][1], ","))[3], ":")[[1]][2])) # get number of free lambda variables
  #n_iterations <- as.numeric(unlist(strsplit(unlist(strsplit(psmc_meta[[3]][1], ","))[1], ":")[[1]][2])) # get number of iterations
  
  iteration <- psmc_lines[(which(psmc_lines == "//")[psmc_obj@n_iterations] + 1): (which(psmc_lines == "//")[psmc_obj@n_iterations + 1] - 1)] # get lines from last iteration
  
  psmc_obj@theta <- as.numeric(strsplit(iteration[6], "\t")[[1]][2]) # get theta 
  psmc_obj@rho <- as.numeric(strsplit(iteration[6], "\t")[[1]][3]) # get rho
  
  curve_lines <- iteration[9:(length(iteration) -1)] # get lines with PSMC values
  curve_values <- lapply(curve_lines, function(x) as.list(unlist(strsplit(x, "\t")))) # split lines by \t
  curve_df <- as.data.frame(do.call("rbind", curve_values)) # combine PSMC values into df
  
  lambda_rows <- lapply(unique(curve_df[,4]), function(x) which(curve_df[,4] == x)[1]) # get indices of start times for each atomic interval
  lambda_df <- curve_df[unlist(lambda_rows),] # get df with start times for each atomic interval
  
  psmc_obj@lk0 <- as.numeric(unlist(lambda_df[1,4])) # get lambda at time 0
  psmc_obj@N0 <- as.numeric((psmc_obj@theta / 4 / psmc_obj@u) * psmc_obj@lk0) # get Ne at time 0
  
  calc_time <- function(t_k, N0, lk0, u){ # t_k = time at interval k, N0 = Ne at time 0, lk0, lambda value at time 0, u = mutation rate
    if (u <= 0){
      time <- t_k
    } else{
      time <- (as.numeric(t_k) * (2 * N0)) / lk0 # calculate time of interval
    }
    
    return(time)
  }
  
  calc_ne <- function(lk, N0, lk0, u){ # lk = lambda at time interval k, N0 = Ne at time 0, lk0 = lambda value at time 0, u = mutation rate
    if (u <= 0){
      ne <- lk
    } else{
      ne <- (as.numeric(lk) * N0) / lk0 # calculate Ne of interval
    }
    
    return(ne)
  }
  
  psmc_obj@base_times <- as.numeric(sapply(seq(1, psmc_obj@n_free_lambdas), function(x) calc_time(lambda_df[x,3], psmc_obj@N0, psmc_obj@lk0, 0))) 
  psmc_obj@base_nes <- as.numeric(sapply(seq(1, psmc_obj@n_free_lambdas), function(x) calc_ne(lambda_df[x,4], psmc_obj@N0, psmc_obj@lk0, 0)))
  psmc_obj@times <- as.numeric(sapply(seq(1, psmc_obj@n_free_lambdas), function(x) calc_time(lambda_df[x,3], psmc_obj@N0, psmc_obj@lk0, psmc_obj@u))) # calculate times
  psmc_obj@nes <- as.numeric(sapply(seq(1, psmc_obj@n_free_lambdas), function(x) calc_ne(lambda_df[x,4], psmc_obj@N0, psmc_obj@lk0, psmc_obj@u))) # calculate Ne values
  
  if (save_output == TRUE){
    param_lines <- c(sprintf("T %s", psmc_obj@theta), sprintf("R %s", psmc_obj@rho), sprintf("N %s", psmc_obj@n_free_lambdas)) # lines for head of file
    est_lines <- mapply(function(x,y) sprintf("H %s\t%s", x, y), psmc_obj@times, psmc_obj@nes, SIMPLIFY = TRUE) # create interval lines
    hist_lines <- append(param_lines, est_lines, after = length(param_lines)) # combine lines
  
    new_file_name <- paste0(psmc_obj@truncated, sprintf("_n%s.history", n_iterations)) # create new file name
  
    history <- file(new_file_name, "w+") # create history file
    writeLines(hist_lines, history) # write history lines
    close(history)
  }
    
  return(psmc_obj)
}


  

  
  
  
  
  
  