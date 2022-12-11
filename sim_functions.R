#This file contains functions used to generate the data for the simulation studies. 

#Function to compute longitudinal covariate 
#Measurement time: t
#Fixed intercept of piecewise constant function: a
#Vector of 10 fixed slopes of piecewise constant function: cVec
#Matrix of 11 random effects for m patients: b.randeff
#Vector of 10 fixed, treatment-dependent slopes of piecewise constant function: dVec
#Matrix of 10 treatment-dependent random effects for m patients: f.randeff
#Matrix of 4 treatments for m patients: trt
Bt = function(t, a, cVec, b.randeff, dVec, f.randeff, trt) {
  a + b.randeff[,1] + (cVec[1]+b.randeff[,2] + trt[,1]*(dVec[1]+f.randeff[,1]))*(t) +
    as.integer(t-1>0)*(cVec[2]+b.randeff[,3] + trt[,1]*(dVec[2]+f.randeff[,2]))*(t-1) +
    as.integer(t-2>0)*(cVec[3]+b.randeff[,4] + trt[,1]*(dVec[3]+f.randeff[,3]))*(t-2) +
    as.integer(t-3>0)*(cVec[4]+b.randeff[,5] + trt[,2]*(dVec[4]+f.randeff[,4]))*(t-3) +
    as.integer(t-4>0)*(cVec[5]+b.randeff[,6] + trt[,2]*(dVec[5]+f.randeff[,5]))*(t-4) +
    as.integer(t-5>0)*(cVec[6]+b.randeff[,7] + trt[,2]*(dVec[6]+f.randeff[,6]))*(t-5) +
    as.integer(t-6>0)*(cVec[7]+b.randeff[,8] + trt[,3]*(dVec[7]+f.randeff[,7]))*(t-6) +
    as.integer(t-7>0)*(cVec[8]+b.randeff[,9] + trt[,3]*(dVec[8]+f.randeff[,8]))*(t-7) +
    as.integer(t-8>0)*(cVec[9]+b.randeff[,10] + trt[,3]*(dVec[9]+f.randeff[,9]))*(t-8)  +
    as.integer(t-9>0)*(cVec[10]+b.randeff[,11] + trt[,4]*(dVec[10]+f.randeff[,10]))*(t-9) 
}


#Function to generate baseline covariates & longitudinal covariates 
#Sequence of measurement times: measTime 
#Sample size: m
#Vector of fixed intercepts for 3 piecewise constant functions: aVec
#Matrix of 10 fixed slopes for 3 piecewise constant functions: cMat
#Matrix of 10 fixed, treatment-dependent slopes for 3 piecewise constant functions: dMat
#Variance-covariance matrix for 11 random effects: D1
#Variance-covariance matrix for 10 treatment-dependent random effects: D2
gen_covs = function(measTime, m, aVec, cMat, dMat, D1, D2) {
  #Generate baseline covariates 
  xVec = runif(m, 0, 1) 
  
  #Generate random effects for each covariate 
  b_i = vector(mode="list", length=3)
  mu1 = rep(0,11)
  for(l in 1:3) {
    b_i[[l]] = mvrnorm(m, mu1, D1)
  }
  f_i = vector(mode="list", length=3)
  mu2 = rep(0,10)
  for(l in 1:3) {
    f_i[[l]] = mvrnorm(m, mu2, D2)
  }
  
  #Matrix to store treatments 
  trtMat = matrix(NA, nrow=m, ncol=4)
  colnames(trtMat) = paste0("trt_",c(1:4))
  #Matrix to store longitudinal covariates WITH measurement errors
  obsMat = array(NA, dim=c(3,m,length(measTime)))
  dimnames(obsMat)[[3]] = paste0("obs_",c(1:length(measTime)))
  #Matrix to store longitudinal covariates WITHOUT measurement errors
  BtMat = array(NA, dim=c(3,m,length(measTime)))
  dimnames(BtMat)[[3]] = paste0("val_",c(1:length(measTime)))
  
  #At each decision point
  for(j in 1:4) {
    #Generate treatment variables
    trtMat[,j] = rbinom(m, 1, 0.5) 
  }
  
  #At each measurement time
  for(j in 1:length(measTime)) {
    #Generate 3 longitudinal covariates
    for(l in 1:3) {
      BtMat[l,,j] = Bt(measTime[j], aVec[l], cMat[,l], b_i[[l]], dMat[,l], f_i[[l]], trtMat)
      obsMat[l,,j] = BtMat[l,,j] + rnorm(m, 0, 0.5)
    }
  }
  
  #Return results
  res = list(xVec, trtMat, obsMat, b_i, f_i, BtMat)
  names(res) = c("base_cov", "trt", "long_cov", "b.randeffs", "f.randeffs", "err_free_long")
  return( res )
}


#Function to compute nu*(t;K)
#Scalar interval: K
#Scalar time: t
#Vector of fixed intercepts for 3 piecewise constant functions: aVec
#Matrix of 10 fixed slopes for 3 piecewise constant functions: cMat
#Matrix of 10 fixed, treatment-dependent slopes for 3 piecewise constant functions: dMat
#Survival time model parameters: betaVec, lambda (lambda=1 for AFT)
#List for 3 biomarkers of 11 random effects for m patients: b_i
#List for 3 biomarkers of 10 treatment-dependent random effects for m patients: f_i
#Vector of m baseline covariates: xVec
#Matrix of m treatments at 10 measurement times: trtMat
nu_K_t = function(K, t, aVec, cMat, dMat, betaVec, lambda=1, b_i, f_i, xVec, trtMat) {
  #Compute summations
  sum0 = 0
  for(l in 1:3) {
    sum0 = sum0 + (betaVec[l]+betaVec[l+5]*trtMat[,K])*(aVec[l]+b_i[[l]][,1])
  }
  
  sum1 = 0
  for(l in 1:3) {
    sum1a = 0
    for(j in 1:K) {
      sum1a = sum1a + (cMat[j,l]+b_i[[l]][,j+1] + trtMat[,K]*(dMat[j,l]+f_i[[l]][,j]))*(j-1)
    }
    sum1 = sum1 + (betaVec[l]+betaVec[l+5]*trtMat[,K])*sum1a
  }
  
  sum2 = 0
  for(l in 1:3) {
    sum2a = 0
    for(j in 1:K) {
      sum2a = sum2a + (cMat[j,l]+b_i[[l]][,j+1] + trtMat[,K]*(dMat[j,l]+f_i[[l]][,j]))
    }
    sum2 = sum2 + (betaVec[l]+betaVec[l+5]*trtMat[,K])*sum2a
  }
  
  #Compute multiplier
  mult = lambda*exp((betaVec[4]+betaVec[9]*trtMat[,K])*xVec + betaVec[5]*trtMat[,K] + sum0 - sum1) / (sum2)
  
  #Compute exp terms
  term1 = exp(t*sum2)
  term2 = exp((K-1)*sum2)
  
  #Compute nu*(t;K)
  return( mult*(term1-term2) )
}


#Function to compute nu(t)
#Scalar time: t
#Vector of fixed intercepts for 3 piecewise constant functions: aVec
#Matrix of 10 fixed slopes for 3 piecewise constant functions: cMat
#Matrix of 10 fixed, treatment-dependent slopes for 3 piecewise constant functions: dMat
#Survival time model parameters: betaVec, lambda (lambda=1 for AFT)
#List for 3 covariates of 11 random effects for m patients: b_i
#List for 3 covariates of 10 treatment-dependent random effects for m patients: f_i
#Vector of m baseline covariates: xVec
#Matrix of m treatments at 10 measurement times: trtMat
nu_t = function(t, aVec, cMat, dMat, betaVec, lambda=1, b_i, f_i, xVec, trtMat) {
  #Identify interval
  K = 1
  for(j in 1:8) {
    K = K + j*as.integer(t>j)*as.integer(t<=(j+1))
  }
  K = K + 9*as.integer(t>9)
  
  #Compute nu(t)
  sum = 0
  if(K>1) {
    for(j in 1:(K-1)) {
      sum = sum + nu_K_t(j, j, aVec, cMat, dMat, betaVec, lambda, b_i, f_i, xVec, trtMat)
    }
  }
  nu_K_t_curr = nu_K_t(K, t, aVec, cMat, dMat, betaVec, lambda, b_i, f_i, xVec, trtMat)
  return( sum + nu_K_t_curr )
}


#Function to compute inverse of nu(t)
#Scalar nu for patient i: nu
#Vector of fixed intercepts for 3 piecewise constant functions: aVec
#Matrix of 10 fixed slopes for 3 piecewise constant functions: cMat
#Matrix of 10 fixed, treatment-dependent slopes for 3 piecewise constant functions: dMat
#Survival time model parameters: betaVec, lambda (lambda=1 for AFT)
#List for 3 biomarkers of 11 random effects for patient i: b_ind
#List for 3 biomarkers of 10 treatment-dependent random effects for patient i: f_ind
#Scalar baseline covariate for patient i: x_ind
#Matrix of treatments for patient i at 10 measurement times: trt_ind
#Vector of 9 knots for patient i: knot_ind
nu_t_inv = function(nu, aVec, cMat, dMat, betaVec, lambda=1, b_ind, f_ind, x_ind, trt_ind, knot_ind) {
  #Identify interval
  R = 1
  for(j in 1:8) {
    R = R + j*as.integer(nu>knot_ind[j])*as.integer(nu<=knot_ind[j+1])
  }
  R = R + 9*as.integer(nu>knot_ind[9])
  
  #Compute inverse of nu_i(t)
  ##Compute summations
  sum0 = 0
  for(l in 1:3) {
    sum0 = sum0 + (betaVec[l]+betaVec[l+5]*trt_ind[,R])*(aVec[l]+b_ind[[l]][,1])
  }
  
  sum1 = 0
  for(l in 1:3) {
    sum1a = 0
    for(j in 1:R) {
      sum1a = sum1a + (cMat[j,l]+b_ind[[l]][,j+1] + trt_ind[,R]*(dMat[j,l]+f_ind[[l]][,j]))*(j-1)
    }
    sum1 = sum1 + (betaVec[l]+betaVec[l+5]*trt_ind[,R])*sum1a
  }
  
  sum2 = 0
  for(l in 1:3) {
    sum2a = 0
    for(j in 1:R) {
      sum2a = sum2a + (cMat[j,l]+b_ind[[l]][,j+1] + trt_ind[,R]*(dMat[j,l]+f_ind[[l]][,j]))
    }
    sum2 = sum2 + (betaVec[l]+betaVec[l+5]*trt_ind[,R])*sum2a
  }
  
  sum3 = 0
  if(R>1) {
    for(j in 1:(R-1)) {
      sum3 = sum3 + nu_K_t(j, j, aVec, cMat, dMat, betaVec, lambda, b_ind, f_ind, x_ind, trt_ind)
    }
  }
  
  ##Compute terms
  term1 = ((nu-sum3)*sum2) / (lambda*exp((betaVec[4]+betaVec[9]*trt_ind[,R])*x_ind + betaVec[5]*trt_ind[,R] + sum0 - sum1))
  term2 = exp((R-1)*sum2)
  
  ##Compute inverse
  nu_inv = log(term1+term2) / sum2
  
  return( nu_inv )
}


#Function to generate survival times from a Cox proportional hazards model with an exponential baseline hazard
#Sample size: m
#Vector of fixed intercepts for 3 piecewise constant functions: aVec
#Matrix of 10 fixed slopes for 3 piecewise constant functions: cMat
#Matrix of 10 fixed, treatment-dependent slopes for 3 piecewise constant functions: dMat
#Cox model parameters: betaVec, lambda 
#List for 3 biomarkers of 11 random effects for m patients: b_i
#List for 3 biomarkers of 10 treatment-dependent random effects for m patients: f_i
#Vector of m baseline covariates: xVec
#Matrix of m treatments at 10 measurement times: trtMat
gen_cox = function(m, aVec, cMat, dMat, betaVec, lambda, b_i, f_i, xVec, trtMat) {
  #Compute knots of linear piecewise function for each patient
  knotMat = matrix(NA, nrow=m, ncol=9)
  for(j in 1:9) {
    knotMat[,j] = nu_t(j, aVec, cMat, dMat, betaVec, lambda, b_i, f_i, xVec, trtMat)
  }
  
  #Generate survival times
  ##Generate random uniform variables
  thetaVec = runif(m, 0, 1) 
  nuVec = -log(1-thetaVec)
  
  ##Compute survival times
  surv_time = rep(NA, m)
  for(i in 1:m) {
    b_i_list = list(b_i[[1]][i,,drop=F], b_i[[2]][i,,drop=F], b_i[[3]][i,,drop=F])
    f_i_list = list(f_i[[1]][i,,drop=F], f_i[[2]][i,,drop=F], f_i[[3]][i,,drop=F])
    surv_time[i] = nu_t_inv(nuVec[i], aVec, cMat, dMat, betaVec, lambda, b_i_list, f_i_list, xVec[i], trtMat[i,,drop=F], knotMat[i,])
  }
  surv_time[which(is.na(surv_time))] = 9999
  
  #Generate censoring times 
  cens_time = runif(m, 0, 55)
  evt_times = cbind(surv_time, cens_time)
  
  #Return results
  return( evt_times )
}


#Function to generate survival times from an accelerated failure time model
#Sample size: m
#Vector of fixed intercepts for 3 piecewise constant functions: aVec
#Matrix of 10 fixed slopes for 3 piecewise constant functions: cMat
#Matrix of 10 fixed, treatment-dependent slopes for 3 piecewise constant functions: dMat
#AFT model parameters: betaVec 
#List for 3 biomarkers of 11 random effects for m patients: b_i
#List for 3 biomarkers of 10 treatment-dependent random effects for m patients: f_i
#Vector of m baseline covariates: xVec
#Matrix of m treatments at 10 measurement times: trtMat
gen_aft = function(m, aVec, cMat, dMat, betaVec, b_i, f_i, xVec, trtMat) {
  #Compute knots of linear piecewise function for each patient
  knotMat = matrix(NA, nrow=m, ncol=9)
  for(j in 1:9) {
    knotMat[,j] = nu_t(j, aVec, cMat, dMat, betaVec, 1, b_i, f_i, xVec, trtMat)
  }
  
  #Generate survival times
  ##Generate random normal variables
  thetaVec = rnorm(m, 3, 1)
  nuVec = exp(thetaVec)
  
  ##Compute survival times
  surv_time = rep(NA, m)
  for(i in 1:m) {
    b_i_list = list(b_i[[1]][i,,drop=F], b_i[[2]][i,,drop=F], b_i[[3]][i,,drop=F])
    f_i_list = list(f_i[[1]][i,,drop=F], f_i[[2]][i,,drop=F], f_i[[3]][i,,drop=F])
    surv_time[i] = nu_t_inv(nuVec[i], aVec, cMat, dMat, betaVec, 1, b_i_list, f_i_list, xVec[i], trtMat[i,,drop=F], knotMat[i,])
  }
  surv_time[which(is.na(surv_time))] = 9999
  
  #Generate censoring times 
  cens_time = runif(m, 0, 55)
  evt_times = cbind(surv_time, cens_time)
  
  #Return results
  return( evt_times )
}


#Function to re-shape dataframe from wide-format to long-format
#List returned from gen_covs: covs
reshape_w2l = function(covs) {
  #Merge treatments and baseline covariates; convert to long format
  wideDF = cbind(covs[[1]], as.data.frame(covs[[2]]))
  wideDF = wideDF[,which(colSums(!is.na(wideDF))>0)]
  names(wideDF)[1] = "base_cov"
  vary_names = names(wideDF)[grepl("trt_",names(wideDF))]
  df = reshape(wideDF, varying=vary_names, v.names="trt", timevar="time", 
               times=seq(0,9,by=3)[1:length(vary_names)], direction="long")
  
  #Merge longitudinal covariates in long format 
  for(i in 1:3) {
    longDF = as.data.frame(covs[[3]][i,,])
    longDF = reshape(longDF, varying=names(longDF), v.names=paste0("long_cov",i), timevar="time",
                     times=seq(0,9,by=1)[1:length(names(longDF))], direction="long")
    df = merge(df, longDF, by=c("id","time"), all=T)
  }
  df = df[order(df$id, df$time),]
  
  #Carry treatments & baseline covariates forward
  df$trt = na.locf(df$trt)
  df$base_cov = na.locf(df$base_cov)
  
  return(df)
}


#Function to create analysis variables
#Matrix of survival/censoring times from gen_aft or gen_cox: timeMat
#Restricted lifetime: L
analysis_vars = function(timeMat, L) {
  #Add id 
  timeDF = as.data.frame(timeMat)
  timeDF$id = c(1:nrow(timeDF))
  
  #Create analysis variables
  timeDF$U = pmin(timeDF$surv_time, timeDF$cens_time) #U=min(T,C)
  timeDF$delta = as.integer(timeDF$surv_time<timeDF$cens_time) #delta=I(T<C)
  timeDF$U_L = pmin(timeDF$U, L) #U_L=min(T,C,L)
  timeDF$delta_L = timeDF$delta + as.integer(L<timeDF$U)*(1-timeDF$delta) #delta_L=I(min(T,L)<C)
  
  return(timeDF)
}


#Function to generate baseline covariates, longitudinal biomarker measurements at baseline, and random effects
#Sample size: m
#Vector of fixed intercepts for 3 piecewise constant functions: aVec
#Matrix of 10 fixed slopes for 3 piecewise constant functions: cMat
#Matrix of 10 fixed, treatment-dependent slopes for 3 piecewise constant functions: dMat
#Variance-covariance matrix for 11 random effects: D1
#Variance-covariance matrix for 10 treatment-dependent random effects: D2
gen_baseline = function(m, aVec, cMat, dMat, D1, D2) {
  #Generate baseline covariates 
  xVec = runif(m, 0, 1) 
  
  #Generate random effects for each biomarker 
  b_i = vector(mode="list", length=3)
  mu1 = rep(0,11)
  for(l in 1:3) {
    b_i[[l]] = mvrnorm(m, mu1, D1)
  }
  f_i = vector(mode="list", length=3)
  mu2 = rep(0,10)
  for(l in 1:3) {
    f_i[[l]] = mvrnorm(m, mu2, D2)
  }
  
  #Matrix to store longitudinal biomarkers WITH measurement errors
  obsMat = array(NA, dim=c(3,m,1))
  dimnames(obsMat)[[3]] = "obs_0"
  #Matrix to store longitudinal biomarkers WITHOUT measurement errors
  BtMat = array(NA, dim=c(3,m,1))
  dimnames(BtMat)[[3]] = "val_0"
  
  #Generate 3 longitudinal biomarker measurements at baseline
  for(l in 1:3) {
    BtMat[l,,1] = aVec[l] + b_i[[l]][,1] 
    obsMat[l,,1] = BtMat[l,,1] + rnorm(m, 0, 0.5)
  }
  
  #Return results
  res = list(xVec, obsMat, b_i, f_i, BtMat)
  names(res) = c("base_cov", "long_cov", "b.randeffs", "f.randeffs", "err_free_long")
  return( res )
}


#Function to generate longitudinal biomarker measurements between decision points dp and dp+1
#Last completed decision point: dp
#Sequence of measurement times between dp and dp+1: measTimeDP
#Sample size: m
#Vector of fixed intercepts for 3 piecewise constant functions: aVec
#Matrix of 10 fixed slopes for 3 piecewise constant functions: cMat
#Matrix of 10 fixed, treatment-dependent slopes for 3 piecewise constant functions: dMat
#List of validation covariate attributes: valid_covs
#Matrix of treatment decisions: trtMat
gen_intermediate = function(dp, measTimeDP, m, aVec, cMat, dMat, valid_covs, trtMat) {
  #Matrix to store longitudinal biomarkers WITH measurement errors
  obsMat = array(NA, dim=c(3,m,length(measTimeDP)))
  dimnames(obsMat)[[3]] = paste0("obs_",measTimeDP)
  #Matrix to store longitudinal biomarkers WITHOUT measurement errors
  BtMat = array(NA, dim=c(3,m,length(measTimeDP)))
  dimnames(BtMat)[[3]] = paste0("val_",measTimeDP)
  
  #Update trtMat for computational purposes
  trtMat[,c((dp+1):4)] = 0
  
  #At each measurement time between dp and dp+1
  for(j in 1:length(measTimeDP)) {
    #Generate 3 longitudinal biomarkers
    for(l in 1:3) {
      BtMat[l,,j] = Bt(measTimeDP[j], aVec[l], cMat[,l], valid_covs$b.randeffs[[l]], dMat[,l], valid_covs$f.randeffs[[l]], trtMat)
      obsMat[l,,j] = BtMat[l,,j] + rnorm(m, 0, 0.5)
    }
  }
  
  #Combine with previous measurements & revert trtMat
  obsMat = abind(valid_covs$long_cov, obsMat, along=3)
  BtMat = abind(valid_covs$err_free_long, BtMat, along=3)
  trtMat[,c((dp+1):4)] = NA
  
  #Return results
  res = list(valid_covs$base_cov, trtMat, obsMat, valid_covs$b.randeffs, valid_covs$f.randeffs, BtMat)
  names(res) = c("base_cov", "trt", "long_cov", "b.randeffs", "f.randeffs", "err_free_long")
  return( res )
}


#Function to add longitudinal biomarker measurements generated between decision points dp and dp+1 to validation data
#Last completed decision point: dp
#List of validation covariate attributes: valid_covs
#Matrix of treatment decisions: trtMat
intermed_meas = function(dp, valid_covs, trtMat) {
  #Longitudinal biomarker parameters
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
  
  #Define measurement times between dp and dp+1
  if(dp==1) {
    measTimeDP = c(1,2,3)
  } else if(dp==2) {
    measTimeDP = c(4,5,6)
  } else {
    measTimeDP = c(7,8,9)
  } 
  
  #Generate longitudinal biomarker measurements between dp and dp+1
  new_covs = gen_intermediate(dp, measTimeDP, length(valid_covs[[1]]), aVec, cMat, dMat, valid_covs, trtMat)
  new_covs_long = reshape_w2l(new_covs) 
  
  #Return results
  return(list(new_covs_long, new_covs))
}

