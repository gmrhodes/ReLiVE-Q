# This file conducts the testing procedure (and optional validation procedure) to estimate the value of 
# a user-specified treatment regime for patients in the simulated data set using a user-specified longitudinal summary vector. 
start_time = Sys.time()
library(anytime)
library(dplyr)
library(foreach)
library(doParallel)
source("D:/Research/functions.R")

########################################## USER-DEFINED INPUT ########################################## 
#Survival model ("aft" or "cox")
surv_mod = "aft"

#Method to construct longitudinal summary 
## baseline, average, lvcf, context_vectors
mainEff = "baseline"

#Regime for which to estimate the value
## 1=optimal, 2=observed, 3=no treatment
regime = 1

#Function for validation procedure
if(regime==1) {
  ## optRegime_validate, or NULL for no validation procedure
  valid_funct = optRegime_validate
} else {
  ## ALWAYS NULL - NEVER validate observed/no treatment regime in this program 
  valid_funct = NULL
}

#Number of cores to parallelize across
num_cores = 20

#Number of iterations
num_iter = 500

#Proportion of data in training data set
train_prop = 0.5

#Path of directory containing simulation data
data_path = "D:/Research/data/"

#If using context vector method, path of directory containing context vectors for simulation data
contVec_path = paste0("D:/Research/context_vectors/", surv_mod, "/")

#If valid_funct is not null, path of directory containing validation data
if(!is.null(valid_funct)) {
  valid_path = "D:/Research/data/"
} else {
  valid_path = NULL
}

#If valid_funct is not null & using the context vector method,
#path of PYTHON PROGRAM to create context vectors for validation data
if(!is.null(valid_funct)) {
  contVec_path_valid = "D:/Research/04_create_validation_contextVecs.py"
} else {
  contVec_path_valid = NULL
}

#Path of directory to write results
result_path = paste0("D:/Research/results/", surv_mod, "/", mainEff, "/") 


########################################## MAIN ########################################## 
#Read in and format data 
df = read.table(paste0(data_path, surv_mod, "_df.csv"), header=T, sep=",")
flatDF = df[which(!duplicated(df$id)),]

#Create vector of decision point times
tVec = c(0,3,6,9)

#Create vector of baseline covariate names
base_covs = c("base_cov")

#Create vector of longitudinal covariate names
long_covs = c("long_cov1","long_cov2","long_cov3")

#If valid_funct is not null
if(!is.null(valid_funct)) {
  #Read in validation data at baseline
  valid_df = read.table(paste0(valid_path, "valid_df_base.csv"), sep=",", header=T)
  valid_covs = readRDS(paste0(valid_path, "valid_covs_base.RDS"))
} else {
  valid_df = NULL
  valid_covs = NULL
}

#Conduct Q-learning on 'num_iter' unique training/testing data divisions  
myCluster = makeCluster(num_cores, type="PSOCK", outfile="") 
registerDoParallel(myCluster)
results = foreach(i=1:num_iter) %dopar% {iterative_fn(i, train_prop, mainEff, regime, contVec_path, df, flatDF, tVec, base_covs, long_covs, 
                                                      valid_funct, valid_df, valid_covs, contVec_path_valid, surv_mod)}
stopCluster(myCluster)

#Create lists to store treatment decisions and pseudo-outcomes
trt_iters_train = vector(mode="list", length=num_iter)
trt_iters_test = vector(mode="list", length=num_iter)
value_iters_train = vector(mode="list", length=num_iter)
value_iters_test = vector(mode="list", length=num_iter)

#Create matrix to store estimated values of the treatment regime 
regime_value_est = matrix(NA, nrow=num_iter, ncol=3)
colnames(regime_value_est) = c("Value_Train", "Value_Test", "Value_Valid")

#Retrieve results
for(i in 1:num_iter) {
  value_iters_train[[i]] = results[[i]][[1]]
  value_iters_test[[i]] = results[[i]][[2]]
  trt_iters_train[[i]] = results[[i]][[3]]
  trt_iters_test[[i]] = results[[i]][[4]]
  regime_value_est[i,] = results[[i]][[5]]
}

#Save results
regime_text = ifelse(regime==1, "opt", ifelse(regime==2, "obs", ifelse(regime==3, "no", "err")))
saveRDS(value_iters_train, paste0(result_path, regime_text, "_value_train.RDS"))
saveRDS(value_iters_test, paste0(result_path, regime_text, "_value_test.RDS"))
saveRDS(trt_iters_train, paste0(result_path, regime_text, "_trt_train.RDS"))
saveRDS(trt_iters_test, paste0(result_path, regime_text, "_trt_test.RDS"))
write.table(regime_value_est, paste0(result_path, regime_text, "_regime_value.csv"), sep=",", row.names=F)

end_time = Sys.time()
run_time = end_time - start_time
print(run_time)
