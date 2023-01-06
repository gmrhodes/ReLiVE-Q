#This file generates baseline data for the validation procedure.
library(MASS)
source("D:/Research/sim_functions.R")
set.seed(6789)


################################################### USER-DEFINED INPUT ################################################### 
start = Sys.time()
#Define output folder where data will be written to
output_path = "D:/Research/data/"

#Define sample size
m = 100000


################################################### MAIN ###################################################
#Longitudinal covariate parameters
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

#Generate baseline covariates, longitudinal measurements at baseline, and random effects
covs = gen_baseline(m, aVec, cMat, dMat, D1, D2)

#Reformat data
long_cov_mat = matrix(t(covs$long_cov[,,1]), nrow=m, ncol=3)
covs_long = as.data.frame(cbind(covs$base_cov, long_cov_mat))
names(covs_long) = c("base_cov","long_cov1","long_cov2","long_cov3")
covs_long$id = c(1:nrow(covs_long))
covs_long$time = 0

#Write dataframe of covariates & associated data for method validation available at Decision Point 1
write.table(covs_long, paste0(output_path,"valid_df_base.csv"), sep=",", row.names=F)
saveRDS(covs, paste0(output_path,"valid_covs_base.RDS"))

end = Sys.time()
print(end-start)
