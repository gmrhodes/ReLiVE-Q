#This file generates the data set for the simulations. 
library(MASS)
library(zoo)
set.seed(5678)
source("D:/Research/sim_functions.R")


################################################### MAIN ################################################### 
start = Sys.time()
#Define output folder where data will be written to
output_path = "D:/Research/data/"

#Define total sample size (training + testing)
m = 10000

#Define restricted lifetime
L = 50

#Define survival time model parameters
aftBetas = c(0.25,-0.25,0.25,0.5,-1.5,-0.25,0.25,0.25,0.25)
coxBetas = c(0.25,-0.5,0.25,0.5,-1.5,-0.5,-0.5,0.25,-0.25)
lambda = 0.1

#Define longitudinal covariate parameters
aVec = c(-2,3,5)
cMat = matrix(NA, nrow=10, ncol=3)
cMat[,1] = c(4,-7,5,-2.5,3.5,-5,1.5,2,-2,1)
cMat[,2] = c(-2,3,1,-3,-1,3.5,-2,-0.5,0.5,0.75)
cMat[,3] = c(-0.5,-2,2,-2,2,3,-2,-1,2,-0.5)
dMat = matrix(NA, nrow=10, ncol=3)
dMat[,1] = rep(0.5, 10)
dMat[,2] = rep(-0.25, 10)
dMat[,3] = c(rep(0.25,5), rep(-0.25,5))
D1 = matrix(0, 11, 11)
diag(D1) = c(1, rep(0.025,5), rep(0.01,5))
D2 = matrix(0, 10, 10)
diag(D2) = c(1, rep(0.025,5), rep(0.01,4))

#Define sequence of measurement times
measTime = seq(0, 9, by=1)

#Generate covariates
covs = gen_covs(measTime, m, aVec, cMat, dMat, D1, D2)
covs_long = reshape_w2l(covs) 

#Create matrix of treatments at 10 measurement times
trt10 = matrix(NA, m, length(measTime))
trt10[,1] = trt10[,2] = trt10[,3] = covs$trt[,1]
trt10[,4] = trt10[,5] = trt10[,6] = covs$trt[,2]
trt10[,7] = trt10[,8] = trt10[,9] = covs$trt[,3]
trt10[,10] = covs$trt[,4]

#Generate and merge survival times from accelerated failure time model
aft_time = gen_aft(m, aVec, cMat, dMat, aftBetas, covs$b.randeffs, covs$f.randeffs, covs$base_cov, trt10)
aft_time = analysis_vars(aft_time, L)
aftDF = merge(covs_long, aft_time, by="id")
write.table(aftDF, paste0(output_path, "aft_df.csv"), sep=",", row.names=F)

#Generate and merge survival times from Cox model
cox_time = gen_cox(m, aVec, cMat, dMat, coxBetas, lambda, covs$b.randeffs, covs$f.randeffs, covs$base_cov, trt10)
cox_time = analysis_vars(cox_time, L)
coxDF = merge(covs_long, cox_time, by="id")
write.table(coxDF, paste0(output_path, "cox_df.csv"), sep=",", row.names=F)

end = Sys.time()
print(end-start)
