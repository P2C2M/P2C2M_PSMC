##### Function to get snp window summary stats #####


#file_name <- "DD14_trunc.psmcfa"

get_snp_windows <- function(file_name){ # file_name = scaffold file
  
  file_con <- file(file_name, "r") # open file connection
  scaff_wins <- list()
  scaff_snps <- list()
  dist <- list()
  scaff_bases <- ""

  readLines(file_con, n=1) # read the first line
  
  while (TRUE){
    line <- readLines(file_con, n = 1) # read in a line
    if (length(line) == 0) { # end of file
      split_bases <- strsplit(scaff_bases, split = "")[[1]] # split string by base
      snps <- which(split_bases == "K") # get positions of snps
      ### Calculate number of windows between snps
      scaff_dist <- list()
      for (i in seq(1,length(snps) - 1)){
        scaff_dist[i] <- snps[i+1] - snps[i] # calculate distances between snp windows
      }
      dist <- append(dist, scaff_dist, after = length(dist))
      
      scaff_wins <- append(scaff_wins, length(split_bases)) # get number of windows
      scaff_snps <- append(scaff_snps, length(snps), after = length(scaff_snps)) # get total number of snp windows
      break
    }else if (grepl(">", line)) { # if line contains scaffold indicator
      split_bases <- strsplit(scaff_bases, split = "")[[1]] # split string by base
      snps <- which(split_bases == "K") # get positions of snps
      ### Calculate number of windows between snps
      scaff_dist <- list()
      for (i in seq(1,length(snps) - 1)){
        scaff_dist[i] <- snps[i+1] - snps[i] # calculate distances between snp windows
      }
      dist <- append(dist, scaff_dist, after = length(dist)) # add scaffold distances to list
      
      scaff_wins <- append(scaff_wins, length(split_bases)) # get number of windows
      scaff_snps <- append(scaff_snps, length(snps), after = length(scaff_snps)) # get total number of snp windows
      
      scaff_bases <- ""
    }else { # if line contains bases
      scaff_bases <- paste0(scaff_bases, line) # add bases to string
      
    }
  }
  
  
  close(file_con)
  
  
  dist <- dist[!is.na(dist)] # remove NA values
  dist_mean <- mean(unlist(dist))
  dist_sd <- sd(unlist(dist))
  dist_range <- max(unlist(dist)) - min(unlist(dist))
  
  scaff_ratio <- unlist(scaff_snps) / unlist(scaff_wins) # get number of snp windows divided by total number of windows
  ratio_mean <- mean(scaff_ratio)
  ratio_sd <- sd(scaff_ratio)
  ratio_range <- max(scaff_ratio) - min(scaff_ratio)
  
  window_stats <- list(dist_mean, dist_sd, dist_range, ratio_mean, ratio_sd, ratio_range) # package stats to return
  
  
  return(window_stats)
  
}