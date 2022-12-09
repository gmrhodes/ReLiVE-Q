#This file computes the validation value estimate for the observed treatment regime and the no treatment regime.
start_time = Sys.time()
set.seed(5678)
library(dplyr)
library(foreach)
library(doParallel)
library(abind)
library(zoo)
source("D:/Research/sim_functions.R")

########################################## USER-DEFINED INPUT ########################################## 
#Number of available cores
num_cores = 65

#Number of iterations
num_iter = 500

#Path of directory containing validation data 
valid_path = "D:/Research/data/"

#Path of directory to write results
result_path = "D:/Research/results/"


########################################## MAIN ########################################## 
#Read in baseline validation data 
valid_covs = readRDS(paste0(valid_path, "valid_covs_base.RDS"))

#Create vector of decision point times
tVec = seq(0,9,by=3)

#Create matrices to store "observed" and no treatment decisions
noTrtMat = matrix(0, nrow=length(valid_covs$base_cov), ncol=length(tVec))
colnames(noTrtMat) = paste0("trt_",c(1:4))
obsTrtMat = matrix(NA, nrow=length(valid_covs$base_cov), ncol=length(tVec))
colnames(obsTrtMat) = paste0("trt_",c(1:4))
valid_covs_no = valid_covs_obs = valid_covs

#Construct biomarker measurements for "observed" (random) and no treatment regime
for(dp in 1:(length(tVec)-1)) {
  #Create "observed" (random) treatment decisions at decision point dp
  obsTrtMat[,dp] = rbinom(length(valid_covs$base_cov), 1, 0.5) 
  
  #Generate longitudinal biomarker measurements collected between decision point dp and dp+1
  ##Observed treatment regime
  valid_upd_obs = intermed_meas(dp, valid_covs_obs, obsTrtMat)
  valid_covs_obs = valid_upd_obs[[2]]
  ##No treatment regime
  valid_upd_no = intermed_meas(dp, valid_covs_no, noTrtMat)
  valid_covs_no = valid_upd_no[[2]]
}
#Create "observed" (random) treatment decisions at last decision point
obsTrtMat[,length(tVec)] = rbinom(length(valid_covs$base_cov), 1, 0.5) 

#Function to compute validation value estimate for observed & no treatment regimes using seed 1233+i
iterative_fn = function(i, tVec, valid_covs_obs, valid_covs_no, obsTrtMat, noTrtMat) {
  source("/home/gmrhodes/Project2/simulations/sim_functions.R")
  source("/home/gmrhodes/Project2/shared_code/functions.R")
  #source("D:/Research/Project2/simulations/sim_functions.R")
  #source("D:/Research/Project2/shared_code/functions.R")
  set.seed(1233+i)
  
  #Compute values for AFT & Cox populations
  ##Observed regime
  obs_vals = rep(NA, 2)
  obs_vals[1] = compute_valid_value(valid_covs_obs, obsTrtMat, "aft", tVec)
  obs_vals[2] = compute_valid_value(valid_covs_obs, obsTrtMat, "cox", tVec)
  ##No treatment regime
  no_trt_vals = rep(NA, 2)
  no_trt_vals[1] = compute_valid_value(valid_covs_no, noTrtMat, "aft", tVec)
  no_trt_vals[2] = compute_valid_value(valid_covs_no, noTrtMat, "cox", tVec)
  
  #Return results
  print(paste("Iteration", i, "complete"))
  return(list(obs_vals, no_trt_vals))
}

#Compute validation value estimate for observed regime & no treatment regime 'num_iter' times
myCluster = makeCluster(num_cores, type="PSOCK", outfile="") 
registerDoParallel(myCluster)
results = foreach(i=1:num_iter) %dopar% {iterative_fn(i, tVec, valid_covs_obs, valid_covs_no, obsTrtMat, noTrtMat)}
stopCluster(myCluster)

#Create matrix to store validation value estimates of observed regime 
obs_valMat = matrix(NA, nrow=num_iter, ncol=2)
colnames(obs_valMat) = c("aft","cox")

#Create matrix to store validation value estimates of no treatment regime 
noTrt_valMat = matrix(NA, nrow=num_iter, ncol=2)
colnames(noTrt_valMat) = c("aft","cox")

#Retrieve results
for(i in 1:num_iter) {
  obs_valMat[i,] = results[[i]][[1]]
  noTrt_valMat[i,] = results[[i]][[2]]
}

#Write results
write.table(obs_valMat, paste0(result_path,"obs_valid_value.csv"), sep=",", row.names=F)
write.table(noTrt_valMat, paste0(result_path,"noTrt_valid_value.csv"), sep=",", row.names=F)

end_time = Sys.time()
run_time = end_time - start_time
print(run_time)
