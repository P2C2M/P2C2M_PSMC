##### Function to read in psmcfa file line by line and calculate the length, #homozyg regions, and #heterozyg regions of each scaffold #####

scaff_length <- function(s = 100, file_name){ # s = window size, file_name = name of empricial psmcfa file
  fa_file <- file(file_name, "r") # open file
  
  ns <- rep(0, 100000)
  homs <- rep(0, 100000)
  hets <- rep(0, 100000)
  scaff_names <- list()
  
  scaff <- 0
  
  while(TRUE){
    line <- readLines(fa_file, n = 1) # read a line
    line_split <- unlist(strsplit(line, "")) # split line by character
    
    if (length(line) == 0){ # if at end of file
      break
      
    }else if (">" %in% line_split == TRUE){
      scaff = scaff + 1
      scaff_names <- append(scaff_names, unlist(strsplit(line, ">"))[2], after = length(scaff_names)) # add scaffold name to list
      
    }else{
      ns[scaff] <- ns[scaff] + length(which(line_split == "N")) # add number of Ns to scaffold count
      homs[scaff] <- homs[scaff] + length(which(line_split == "T")) # add number of Ts to scaffold count
      hets[scaff] <- hets[scaff] + length(which(line_split == "K")) # add number of Ks to scaffold count
    }
  }
  
  
  close(fa_file)
  
  ns <- ns[1:scaff] # remove unused 0s from list
  homs <- homs[1:scaff] # remove unused 0s from list
  hets <- hets[1:scaff] # remove unused 0s from list
  
  total_wins <- ns + homs + hets # get total number of windows per scaffold
  total_bp <- total_wins * s # get total number of bps per scaffold *Note* since scaffold length is not always neatly divisible by s, these scaffold lengths can be off by as much as (s-1)bp
  
  return(total_bp)
}
