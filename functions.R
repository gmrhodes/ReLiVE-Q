#This file contains functions to implement the RL Q-learning method.

#Function to construct dataframe of main effects for baseline method
#Resulting dataframe includes id, base_covs, long_covs
baseline_mainEff = function(df, tVec, dp, base_covs, long_covs) {
  df_t = df[which((df$time<=tVec[dp]) & (df$U_L>tVec[dp])),]
  baseDF = df_t[which(df_t$time==0), c("id",base_covs,long_covs)]
  return(baseDF)
}


#Function to construct dataframe of main effects for average method
#Resulting dataframe includes id, base_covs, long_covs
average_mainEff = function(df, tVec, dp, base_covs, long_covs) {
  df_t = df[which((df$time<=tVec[dp]) & (df$U_L>tVec[dp])),]
  avgDF = df_t %>% group_by(id) %>% summarise_at(vars(all_of(long_covs)), mean, na.rm=T)
  avgDF = merge(avgDF, df_t[which(df_t$time==0),c("id",base_covs)], by="id")
  return(avgDF)
}


#Function to construct dataframe of main effects for last-value carried forward method
#Resulting dataframe includes id, base_covs, long_covs 
lvcf_mainEff = function(df, tVec, dp, base_covs, long_covs) {
  df_t = df[which((df$time<=tVec[dp]) & (df$U_L>tVec[dp])),]
  df_t = df_t[order(df_t$id,df_t$time),]
  df_t[,c(base_covs,long_covs)] = apply(df_t[,c(base_covs,long_covs)], 2, na.locf)
  lvcfDF = subset(df_t %>% group_by(id) %>% top_n(n=1, wt=time), select=c("id",base_covs,long_covs))
  return(lvcfDF)
}


#Function to construct dataframe of main effects for context vector method
#Resulting dataframe includes id, base_covs, context vector elements for long_covs
contextVec_mainEff = function(contVec_path, tVec, dp, base_covs, long_covs) {
  contVecDF = read.table(paste0(contVec_path, "contextVec_t", tVec[dp], ".csv"), sep=",", header=T) 
  pattern = paste(long_covs, collapse="|")
  long_covs_names = names(contVecDF)[grepl(pattern, names(contVecDF))]  
  contVecDF = contVecDF[,c("id",base_covs,long_covs_names)]
  return(contVecDF)
}


#Function to compute Kaplan-Meier estimate of P(C > U_L)
#Resulting dataframe includes id, prob_cGTUl
prob_C_gt_Ul = function(df) {
  df$censored_ind = 1 - df$delta
  kmf = survfit(Surv(df$U, df$censored_ind) ~ 1, type="kaplan-meier")
  kmfSumm = summary(kmf, times=df$U_L)
  kmfDF = data.frame(days=kmfSumm$time, prob=kmfSumm$surv)
  prob_cGTUl = kmfDF[match(df$U_L, kmfDF$days),"prob"]
  probDF = data.frame("id"=df$id, "prob_cGTUl"=prob_cGTUl)
  return( probDF )
}


#Function to compute Kaplan-Meier estimate of P(C > tau_k)
#Returns scalar value
prob_C_gt_tauk = function(df, tau_k) {
  df$censored_ind = 1 - df$delta
  kmf = survfit(Surv(df$U, df$censored_ind) ~ 1, type="kaplan-meier")
  tauEst = summary(kmf, times=tau_k)$surv
  return( tauEst )
}


#Function to compute the estimated optimal treatment decision for each patient in predDF 
#Resulting dataframe includes id, trt
optimal_trt_dec = function(predDF) {
  predDF$trt = as.integer(predDF$pred1 > predDF$pred0)
  return(predDF[,c("id","trt")])
}


#Function to save observed treatment decision for each patient in designDF
#Resulting dataframe includes id, trt
observed_trt_dec = function(designDF) {
  resultDF = data.frame(id=designDF$id, trt=designDF$trt)
  return(resultDF)
}


#Function to generate treatment=0 (no treatment) for each patient in designDF
#Resulting dataframe includes id, trt
no_trt_dec = function(designDF) {
  resultDF = data.frame(id=designDF$id, trt=rep(0,nrow(designDF)))
  return(resultDF)
}


#Function to construct dataframe of main effects for validation data using baseline method
#Resulting dataframe includes id, base_covs, long_covs
baseline_mainEff_valid = function(df, tVec, dp, base_covs, long_covs) {
  baseDF = df[which(df$time==0), c("id",base_covs,long_covs)]
  return(baseDF)
}


#Function to construct dataframe of main effects for validation data using average method
#Resulting dataframe includes id, base_covs, long_covs
average_mainEff_valid = function(df, tVec, dp, base_covs, long_covs) {
  df_t = df[which(df$time<=tVec[dp]),]
  avgDF = merge(df_t %>% group_by(id) %>% summarise_at(vars(all_of(long_covs)), mean),
                df_t[which(df_t$time==0),c("id","base_cov")], by="id")
  return(avgDF)
}


#Function to construct dataframe of main effects for validation data using last-value carried forward method
#Resulting dataframe includes id, base_covs, long_covs
lvcf_mainEff_valid = function(df, tVec, dp, base_covs, long_covs) {
  df_t = df[which(df$time<=tVec[dp]),]
  lvcfDF = subset(df_t %>% group_by(id) %>% top_n(n=1, wt=time), select=c("id",base_covs,long_covs))
  return(lvcfDF)
}


#Function to compute validation value estimate for patients in valid_covs using treatments in trtMat 
#Returns scalar validation value estimate
compute_valid_value = function(valid_covs, trtMat, surv_mod, tVec) {
  #Define simulation parameters
  ##Restricted lifetime
  L = 50
  ##Survival time model parameters
  aftBetas = c(0.25,-0.25,0.25,0.5,-1.5,-0.25,0.25,0.25,0.25)
  coxBetas = c(0.25,-0.5,0.25,0.5,-1.5,-0.5,-0.5,0.25,-0.25)
  lambda = 0.1
  ##Longitudinal biomarker parameters
  aVec = c(-2,3,5)
  cMat = matrix(NA, nrow=10, ncol=3)
  cMat[,1] = c(4,-7,5,-2.5,3.5,-5,1.5,2,-2,1)
  cMat[,2] = c(-2,3,1,-3,-1,3.5,-2,-0.5,0.5,0.75)
  cMat[,3] = c(-0.5,-2,2,-2,2,3,-2,-1,2,-0.5)
  dMat = matrix(NA, nrow=10, ncol=3)
  dMat[,1] = rep(0.5, 10)
  dMat[,2] = rep(-0.25, 10)
  dMat[,3] = c(rep(0.25,5), rep(-0.25,5))
  
  #Create matrix of treatments at 10 measurement times
  trt10 = matrix(NA, length(valid_covs$base_cov), 10)
  trt10[,1] = trt10[,2] = trt10[,3] = trtMat[,1]
  trt10[,4] = trt10[,5] = trt10[,6] = trtMat[,2]
  trt10[,7] = trt10[,8] = trt10[,9] = trtMat[,3]
  trt10[,10] = trtMat[,4]
  
  #Generate survival times corresponding to treatments in trtMat
  if(surv_mod=="cox") {
    timeDF = gen_cox(length(valid_covs$base_cov), aVec, cMat, dMat, coxBetas, lambda, valid_covs$b.randeffs, valid_covs$f.randeffs, valid_covs$base_cov, trt10)
  } else {
    timeDF = gen_aft(length(valid_covs$base_cov), aVec, cMat, dMat, aftBetas, valid_covs$b.randeffs, valid_covs$f.randeffs, valid_covs$base_cov, trt10)
  }
  timeDF = as.data.frame(timeDF)
  timeDF$min_Td_L = pmin(timeDF$surv_time, L) 
  
  #Compute each patient's value 
  val = rep(0, nrow(timeDF))
  for(dp in 1:length(tVec)) {
    val = val + as.integer(timeDF$min_Td_L>tVec[dp])*(timeDF$min_Td_L-tVec[dp])
  }
  
  #Compute & return validation value for regime
  return( mean(val) )
}


#Function to compute validation value estimate for estimated optimal treatment regime 
#Returns scalar validation value estimate
optRegime_validate = function(qModel_trt0, qModel_trt1, valid_df, valid_covs, contVec_path_valid, 
                                     mainEff_funct_valid, tVec, base_covs, long_covs, surv_mod, contVec_path, df) {
  source("D:/Research/sim_functions.R")
  
  #If using context vector method
  if(identical(mainEff_funct_valid, "contVec_mainEff_valid")) {
    #Create variables for scaling
    minVec = maxVec = rep(NA, length(long_covs))
    for(i in 1:length(long_covs)) {
      minVec[i] = min(df[,paste0("long_cov",i)])
      maxVec[i] = max(df[,paste0("long_cov",i)])
    }
    
    #Source python code to create context vectors
    use_python("/usr/bin/python3")
    source_python(contVec_path_valid)
  }
  
  #Create a matrix to store estimated optimal treatment decisions
  trtMat = matrix(NA, nrow=length(unique(valid_df$id)), ncol=length(tVec))
  colnames(trtMat) = paste0("trt_",c(1:length(tVec)))
  
  #Compute estimated optimal treatment for each decision point
  for(dp in 1:length(tVec)) {
    #Create design matrix
    if(dp==1) {
      #Main effects equal the observed values at decision point 1 (i.e. baseline)
      designDF = baseline_mainEff_valid(valid_df, tVec, dp, base_covs, long_covs)
    } else {
      #Main effects are derived according to user-specified function
      if(identical(mainEff_funct_valid, "contVec_mainEff_valid")) {
        designDF = contVecDF
      } else {
        designDF = mainEff_funct_valid(valid_df, tVec, dp, base_covs, long_covs)
      }
    }
    designDF = designDF[order(designDF$id),]
    
    #Predict values for Trt 0 & Trt 1
    pred0 = predict(qModel_trt0[[dp]], newdata=designDF)
    pred1 = predict(qModel_trt1[[dp]], newdata=designDF)
    
    #Estimate optimal treatment decision
    trtMat[,dp] = as.integer(pred1 > pred0)
    
    #If not the last decision point
    if(dp<length(tVec)) {
      #Generate longitudinal biomarker measurements collected between decision point dp and dp+1
      valid_upd = intermed_meas(dp, valid_covs, trtMat)
      valid_df = valid_upd[[1]]
      valid_covs = valid_upd[[2]]
      
      #If using context vector method
      if(identical(mainEff_funct_valid, "contVec_mainEff_valid")) {
        #Construct context vectors for decision point dp+1 using fitted LSTM autoencoder
        contVecDF = create_valid_contVecs(valid_df, tVec[dp+1], contVec_path, minVec, maxVec)
        
        #Format context vector dataframe
        pattern = paste(long_covs, collapse="|")
        long_covs_names = names(contVecDF)[grepl(pattern, names(contVecDF))]
        contVecDF = contVecDF[,c("id",base_covs,long_covs_names)]
      }
    }
  }
  
  #Compute validation value for estimated optimal treatment regime
  value_valid = compute_valid_value(valid_covs, trtMat, surv_mod, tVec) 
  return(value_valid)
}


#Function to compute random forest dependent variable for dataframe df
#Returns dataframe df with new column 'outcome' containing the dependent variable
randomForest_outcome = function(df, prob_cGTUl_df, tVec, dp, resultDF=NULL) {
  #Compute censoring weight for first summand, scriptK(U_L) 
  df = merge(df, prob_cGTUl_df, by="id")
  prob_cGTtauk = prob_C_gt_tauk(df, tVec[dp])
  df$scriptK_Ul = df$prob_cGTUl / prob_cGTtauk
  
  #Compute first summand
  df$summand1 = (df$delta_L * df$RRL) / df$scriptK_Ul
  
  if(dp!=length(tVec)) { #If not at last decision point
    #Compute censoring weight for second summand, scriptK(tau_k+1) 
    prob_cGTtauk1 = prob_C_gt_tauk(df, tVec[dp+1])
    df$scriptK_tauk1 = prob_cGTtauk1 / prob_cGTtauk
    
    #Merge pseudo-outcomes, vTilde_kPlus1
    names(resultDF) = c("id", "trt_k1", "v_k1")
    df = merge(df, resultDF, by="id", all=T)
    
    #Create indicator I(U_L > tau_k+1)
    df$ind_Ul_gt_tauk1 = as.integer(df$U_L > tVec[dp+1])
    
    #Compute second summand
    df$summand2 = (df$ind_Ul_gt_tauk1 * df$v_k1) / df$scriptK_tauk1
    df[which(df$ind_Ul_gt_tauk1==0),"summand2"] = 0
    
  } else {
    #Set summand2=0 for last decision point
    df$summand2 = 0
  }
  
  #Compute random forest dependent variable
  df$outcome = df$summand1 + df$summand2
  return(df)
}


#Function to estimate value of regime using random forest Q-models
#Resulting train & test dataframes include id, trt, value; also returns training Q-models
randomForest_qModels = function(designDF, trainDF, testDF, base_covs, long_covs_names, regime, dp, tVec) {
  #Merge train/test design matrices with outcomes, treatments, and risk indicators
  designTrain = merge(designDF, trainDF[,c("id","ind_Ul_gt_tauk","trt","outcome")], by="id", all.x=F, all.y=T)
  designTest = merge(designDF, testDF[,c("id","ind_Ul_gt_tauk","trt","outcome")], by="id", all.x=F, all.y=T)
  
  #Remove patients not at risk & sort by id
  riskTrain = designTrain[which(designTrain$ind_Ul_gt_tauk==1),]
  riskTrain = riskTrain[order(riskTrain$id),]
  riskTest = designTest[which(designTest$ind_Ul_gt_tauk==1),]
  riskTest = riskTest[order(riskTest$id),]
  
  #Estimate Q-function for Trt 0 via random forest using training data 
  riskTrain0 = riskTrain[which(riskTrain$trt==0),]
  form = as.formula(paste("outcome ~ ", paste(c(base_covs, long_covs_names), collapse=" + ")))
  randfor0 = randomForest(form, data=riskTrain0, ntree=250, mtry=ceiling((length(long_covs_names)+1)/3), 
                          sampsize=nrow(riskTrain0), replace=T, nodesize=50)
  
  #Estimate Q-function for Trt 1 via random forest using training data 
  riskTrain1 = riskTrain[which(riskTrain$trt==1),]
  randfor1 = randomForest(form, data=riskTrain1, ntree=250, mtry=ceiling((length(long_covs_names)+1)/3), 
                          sampsize=nrow(riskTrain1), replace=T, nodesize=50)
  
  #Compute predicted outcome corresponding to Trt 0 for training & testing risk sets using training Q-model
  predTrain0.0 = randfor0$predicted #OOB predictions
  predTrain0.1 = predict(randfor0, newdata=riskTrain1)
  predTest0 = predict(randfor0, newdata=riskTest)
  
  #Compute predicted outcome corresponding to Trt 1 for training & testing risk sets using training Q-model
  predTrain1.0 = predict(randfor1, newdata=riskTrain0)
  predTrain1.1 = randfor1$predicted #OOB predictions
  predTest1 = predict(randfor1, newdata=riskTest)
  
  #Create dataframes of predicted pseudo-outcomes
  trainPreds = data.frame(id = c(riskTrain0$id, riskTrain1$id), 
                          pred0 = c(predTrain0.0, predTrain0.1),
                          pred1 = c(predTrain1.0, predTrain1.1))
  trainPreds = trainPreds[order(trainPreds$id),]
  testPreds = data.frame(id=riskTest$id, pred0=predTest0, pred1=predTest1)
  
  #Compute treatment decision for each patient in training & testing risk sets according to user-specified regime
  if(regime==1) {
    #Compute estimated optimal treatment
    resultTrain = optimal_trt_dec(trainPreds)
    resultTest = optimal_trt_dec(testPreds)
    
  } else if(regime==2) {
    #Save observed treatment
    resultTrain = observed_trt_dec(riskTrain)
    resultTest = observed_trt_dec(riskTest)
    
  } else if(regime==3) {
    #Generate treatment=0 (no treatment)
    resultTrain = no_trt_dec(riskTrain)
    resultTest = no_trt_dec(riskTest)
  }
  
  #Estimate value of the treatment decision for each patient in training risk set 
  resultTrain = merge(resultTrain, trainPreds, by="id")
  resultTrain$value = ifelse(resultTrain$trt==1, resultTrain$pred1, resultTrain$pred0)
  
  #Estimate Q-function for Trt 0 via random forest using testing data 
  riskTest0 = riskTest[which(riskTest$trt==0),]
  randfor0.test = randomForest(form, data=riskTest0, ntree=250, mtry=ceiling((length(long_covs_names)+1)/3), 
                               sampsize=nrow(riskTest0), replace=T, nodesize=50)
  
  #Estimate Q-function for Trt 1 via random forest using testing data 
  riskTest1 = riskTest[which(riskTest$trt==1),]
  randfor1.test = randomForest(form, data=riskTest1, ntree=250, mtry=ceiling((length(long_covs_names)+1)/3), 
                               sampsize=nrow(riskTest1), replace=T, nodesize=50)
  
  #Compute predicted outcome corresponding to Trt 0 for testing risk set using testing Q-model
  predTest0.0.test = randfor0.test$predicted #OOB predictions
  predTest0.1.test = predict(randfor0.test, newdata=riskTest1)
  
  #Compute predicted outcome corresponding to Trt 1 for testing risk set using testing Q-model
  predTest1.0.test = predict(randfor1.test, newdata=riskTest0)
  predTest1.1.test = randfor1.test$predicted #OOB predictions
  
  #Create dataframe of predicted pseudo-outcomes from Q-model
  testPreds.test = data.frame(id = c(riskTest0$id, riskTest1$id), 
                              pred0 = c(predTest0.0.test, predTest0.1.test),
                              pred1 = c(predTest1.0.test, predTest1.1.test))
  testPreds.test = testPreds.test[order(testPreds.test$id),]
  
  #Estimate value of the treatment decision for each patient in testing risk set
  resultTest = merge(resultTest, testPreds.test, by="id")
  resultTest$value = ifelse(resultTest$trt==1, resultTest$pred1, resultTest$pred0)
  
  #Return result matrices of id, trt, value, as well as the fitted training Q-models
  return(list(resultTrain[,c("id","trt","value")], resultTest[,c("id","trt","value")], randfor0, randfor1))
}


#Function to conduct Q-learning on a single train/test split using seed 1233+i
#Returns dataframes of training & testing pseudo-outcomes and treatment decisions
#Returns scalar training, testing, and validation value estimates
iterative_fn = function(i, train_prop, mainEff, regime, contVec_path, df, flatDF, tVec, base_covs, long_covs, 
                        valid_funct=NULL, valid_df=NULL, valid_covs=NULL, contVec_path_valid=NULL, valid_param=NULL) {
  library(survival)
  library(randomForest)
  library(splitstackshape)
  library(data.table)
  library(anytime)
  library(dplyr)
  library(abind)
  library(reticulate)
  library(zoo)
  source("D:/Research/functions.R")
  set.seed(1233+i)
  
  #Define main effect function
  mainEff_funct = ifelse(mainEff=="baseline", baseline_mainEff, 
                         ifelse(mainEff=="average", average_mainEff, 
                                ifelse(mainEff=="lvcf", lvcf_mainEff, 
                                       ifelse(mainEff=="context_vectors", contextVec_mainEff, "error"))))
  
  #Retrieve ids of patients in training & testing sets
  trainIDs = sample(flatDF$id, size=round(nrow(flatDF)*train_prop), replace=F)
  testIDs = flatDF[which(!(flatDF$id %in% trainIDs)),"id"]
  
  #Create dataframes to store treatment decisions 
  trtDF_train = data.frame(id=trainIDs)
  trtDF_test = data.frame(id=testIDs)
  
  #Create dataframes to store value estimates
  valueDF_train = data.frame(id=trainIDs)
  valueDF_test = data.frame(id=testIDs)
  
  #Create list to store random forest Q-models for Trt 0 & Trt 1
  qModel_trt0 = vector(mode="list", length=length(tVec))
  qModel_trt1 = vector(mode="list", length=length(tVec))
  
  #Compute Kaplan-Meier estimate of P(C > U_L) for training & testing data 
  prob_cGTUl_train = prob_C_gt_Ul(flatDF[which(flatDF$id %in% trainIDs),])
  prob_cGTUl_test = prob_C_gt_Ul(flatDF[which(flatDF$id %in% testIDs),])
  
  #For k=K,K-1,...,1
  for(dp in c(length(tVec):1)) {
    #Create dataframe of id, U, delta, U_L, delta_L, and treatment administered at decision point k
    dpDF = df[which(df$time==tVec[dp]),c("id","U","delta","U_L","delta_L","trt")]
    
    #Create indicator I(U_L > tau_k)
    dpDF$ind_Ul_gt_tauk = as.integer(dpDF$U_L > tVec[dp])
    
    #Compute restricted residual life
    dpDF$RRL = ifelse(dpDF$ind_Ul_gt_tauk==1, (dpDF$U_L - tVec[dp]), NA) 
    
    #Divide data into a training & testing sets 
    trainDF = dpDF[which(dpDF$id %in% trainIDs),]
    testDF = dpDF[which(dpDF$id %in% testIDs),]
    
    #Compute random forest dependent variable
    if(dp!=length(tVec)) {
      trainDF = randomForest_outcome(trainDF, prob_cGTUl_train, tVec, dp, resultTrain)
      testDF = randomForest_outcome(testDF, prob_cGTUl_test, tVec, dp, resultTest)
    } else {
      trainDF = randomForest_outcome(trainDF, prob_cGTUl_train, tVec, dp)
      testDF = randomForest_outcome(testDF, prob_cGTUl_test, tVec, dp)
    }
    
    #Create dataframe of main effects
    if(dp==1) {
      #Main effects equal the observed values at decision point 1 (i.e. baseline)
      designDF = baseline_mainEff(df, tVec, dp, base_covs, long_covs)
    } else {
      #Main effects are derived according to user-specified function
      if(identical(mainEff_funct, contextVec_mainEff)) {
        designDF = mainEff_funct(contVec_path, tVec, dp, base_covs, long_covs)
      } else {
        designDF = mainEff_funct(df, tVec, dp, base_covs, long_covs)
      }
    }
    
    #Define variable names of longitudinal covariates
    pattern = paste(long_covs, collapse="|")
    long_covs_names = sort(names(designDF)[grepl(pattern, names(designDF))])
    
    #Compute pseudo-outcomes for regime using random forest Q-models
    results = randomForest_qModels(designDF, trainDF, testDF, base_covs, long_covs_names, regime, dp, tVec) 
    resultTrain = results[[1]]
    resultTest = results[[2]]
    qModel_trt0[[dp]] = results[[3]]
    qModel_trt1[[dp]] = results[[4]]
    rm(designDF, trainDF, testDF)
    
    #Save treatment decisions & estimated pseudo-outcomes
    ##Training
    names(resultTrain)[2:3] = paste0(names(resultTrain)[2:3], "_dp", dp)
    trtDF_train = merge(trtDF_train, resultTrain[,c(1,2)], by="id", all=T)
    valueDF_train = merge(valueDF_train, resultTrain[,c(1,3)], by="id", all=T)
    ##Testing
    names(resultTest)[2:3] = paste0(names(resultTest)[2:3], "_dp", dp)
    trtDF_test = merge(trtDF_test, resultTest[,c(1,2)], by="id", all=T)
    valueDF_test = merge(valueDF_test, resultTest[,c(1,3)], by="id", all=T)
    rm(results)
  }
  
  #Compute train/test value estimate for regime
  value_train = mean(valueDF_train$value_dp1)
  value_test = mean(valueDF_test$value_dp1)
  
  #Compute validation value estimate for regime
  if(!is.null(valid_funct)) {
    #Define main effect function
    mainEff_funct_valid = ifelse(mainEff=="baseline", baseline_mainEff_valid, 
                                 ifelse(mainEff=="average", average_mainEff_valid, 
                                        ifelse(mainEff=="lvcf", lvcf_mainEff_valid, 
                                               ifelse(mainEff=="context_vectors", "contVec_mainEff_valid", "error"))))
    #Run validation procedure
    value_valid = valid_funct(qModel_trt0, qModel_trt1, valid_df, valid_covs, contVec_path_valid, 
                              mainEff_funct_valid, tVec, base_covs, long_covs, valid_param, contVec_path, df)
  } else {
    value_valid = NA
  }
  
  #Return results
  print(paste("Iteration", i, "complete"))
  return( list(valueDF_train, valueDF_test, trtDF_train, trtDF_test, c(value_train,value_test,value_valid)) )
}
