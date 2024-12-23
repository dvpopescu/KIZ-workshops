# The function below, 'occu_data' simulates ocupancy datasets given the following 
# parameters: 
# - occupancy probability - psi
# - detection probability - p
# - number of sites - n_sites
# - number of visits - n_vis

# The function returns simulated true and observed occupancy values, given the 
# psi and p probabilities. The function also returns several other statistics
# about the performance of the occupancy monitoring program 

# RUN the function below; DO NOT CHANGE ANYTHING!!!!!! 

occu_data <- function(psi, p, n_sites, n_vis) {
  
  tocc <- rep(0, n_sites)
  tocc <- as.data.frame(tocc)
  names(tocc) <- c("T1")
  tocc
  
  ### set up an empty observed occupancy matrix (OBSERVATION PROCESS)
  obsocc <- matrix(rep(0, n_sites*n_vis), ncol=(n_vis),
                   dimnames=list(paste("site", 1:n_sites, sep=""), 
                                 paste("v", 1:n_vis, sep="")))
  obsocc <- as.data.frame(obsocc)
  obsocc
  
  # fill the observed occupancy matrix
  for(i in 1:n_sites) {
    tocc[i,1] = rbinom(1, 1, psi)                       # Simulate site occupancy at t1
    obsocc[i,1:n_vis] = rbinom(n_vis, 1, tocc[i,1]*p)   # Simulate observed P/A at t1 using detection prob. p
  }
  
  return(list(obsocc, tocc,
              
              print(paste("true number of sites occupied  =", sum(tocc))),
              
              print(paste("when detection probability = ", p, "and number of visits =", n_vis,
                          ", the number of sites with at least 1 detection =", length(which(rowSums(obsocc)>0)))),
              
              print(paste("the number of sites missed due to low survey effort is = ", 
                          sum(tocc) - length(which(rowSums(obsocc)>0)))),
              
              print(paste("when detection probability = ", p, "and number of visits =", n_vis,
                          ", the total number of detections =", sum(obsocc))),
              
              print(paste("when detection probability = ", p, "and number of visits =", n_vis,
                          "the proportion of ocassions with detections is", round(sum(obsocc)/(sum(tocc)*n_vis), 3)))
              
  ))
  
}


# EXERCISE: using a set psi and p, change the number of visits (the last parameter 
# of the function) to evaluate how many occupied sites would be completely missed 
# (no detections) due to low survey effort. 
# Enter a different value for p  (second parameter of the function) and repeat the 
# procedure above

# example
occu_data(psi = 0.3, p = 0.3, n_sites = 100, n_vis = 7)



