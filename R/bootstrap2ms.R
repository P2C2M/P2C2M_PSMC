### Function to convert psmc output from bootstrap file to psmc objects with ms values ###


# boots_file <- "combined.psmc"
# n_iterations <- 25
# gen <- 1
# u <- 1.1e-08
# n = 2
# L = 30000000
# s = 100
# R = 10
# r = 1

bootstrap2ms <- function(boots_file, n_iterations, gen, u, n = 2, L = 30000000, s, R = 10, r = 1){

  psmc_file <- file(boots_file, "r+")
  psmc_lines <- readLines(psmc_file) # read in psmc results and all bootstraps
  close(psmc_file)
  
  
  cc_lines <- which(psmc_lines == "CC")[c(TRUE,FALSE)] # get first line of each PSMC rep
  end_lines <- (cc_lines - 1)[2:length(cc_lines)] # get end line for each PSMC rep
  end_lines <- append(end_lines, length(psmc_lines), after = length(end_lines)) # add in final line
  
  boots <- mapply(function(x,y) psmc_lines[x:y], cc_lines, end_lines) # separate bootstraps
  
  boots_objs <- apply(boots, MARGIN = 2, function(x) psmc2history(file_name = x, n_iterations = n_iterations, gen = gen, u = u, from_file = FALSE, save_output = FALSE)) # calculate history for each boot
  
  boots_objs <- lapply(boots_objs, function(x) psmc_object2ms(psmc_object = x, n = n, L = L, s = s, R = R, r = r, print_cmd = FALSE)) # calculate ms values
  
  #empirical_psmc <- boots_objs[1]
  #boots_objs <- boots_objs[2:length(boots_objs)] # remove empirical results
  
  return(boots_objs)
}






