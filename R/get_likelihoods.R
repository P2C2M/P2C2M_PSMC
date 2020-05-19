##### Function to get likelihood values from PSMC results #####

#file_name <- sim_file_list[1]

get_likelihoods <- function(file_name){ # file_name = psmc results file

  file_con <- file(file_name, "r") # open file connection
  likelihood_lines <- c()
  
  while (TRUE){
    line <- readLines(file_con, n = 1) # read in a line
    if (length(line) == 0) { # end of file
      break
    }else if (grepl("LK *", line)) { # if line contains likelihood indicator
      likelihood_lines <- c(likelihood_lines,line)
    }
  }
  
  close(file_con)
  
  likelihood_lines <- likelihood_lines[2:length(likelihood_lines)] # get rid of likelihood 0
  
  likelihoods <- as.numeric(lapply(likelihood_lines, function(x) as.numeric(strsplit(x, "LK\t")[[1]][2]))) # isolate likelihood values
  
  return(likelihoods)

}