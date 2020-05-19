##### Function to compare one empirical value to many simulated values #####

compare_onetomany <- function(emp_value, sim_values){
  
  upper <- length(which(unlist(sim_values) > emp_value)) # calculate how many simulated values are >= empirical value
  lower <- length(which(unlist(sim_values) <= emp_value)) # calculate how many simulated values are <= empirical value
  p_violations <- min(c(upper,lower)) * 2 /length(sim_values) # treat as two-tailed test
  
  return(p_violations)
}