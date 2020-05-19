### Function for writing histograms of summary statistics ###

hist_writer <- function(summary_statistic, emp_value, sim_values){

  grDevices::pdf(file = paste0(summary_statistic, ".pdf")) # start pdf
  graphics::hist(unlist(sim_values), breaks = 20, main = summary_statistic, xlab = "") # plot histogram
  graphics::abline(v = emp_value, col = "red") # add line for posterior comparison
  grDevices::dev.off() # close pdf

}
