# This file contains functions to create the window-specific context vectors for the validation data.
#################################################### SET-UP ####################################################
import os
os.environ['PYTHONHASHSEED']='0'
os.environ['TF_CPP_MIN_LOG_LEVEL']='2'
import random
random.seed(5678)
import numpy as np
np.random.seed(5678)
import tensorflow as tf
tf.random.set_seed(5678)
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
import pandas as pd

#Padding value for missing time steps
padVal = -9999 

#LSTM autoencoder hyperparameters
epochs = 500 #training epochs
cv_dimension = 5 #context vector dimension


#################################################### FUNCTION TO SCALE COVARIATES ####################################################
def minMaxScale(varVec, minVal, maxVal):
    #Min-max scale varVec given minVal and maxVal
    scaleVec = (varVec-minVal) / (maxVal-minVal)
    return scaleVec


#################################################### FUNCTION TO CREATE 3D COVARIATE ARRAY ####################################################
def covArray_construct(covDF):    
    #Create padded 3D numpy arrays of covariates for Keras
    gb = covDF.groupby(['id'])
    maxBlocks = gb['id'].size().max() 
    covArr = np.array([np.pad(frame['time'].values, pad_width=(0, maxBlocks-len(frame)), mode='constant', constant_values=padVal) for _,frame in gb]).reshape(-1, maxBlocks, 1)
    for col in covDF.columns[2:]:
        newArr = np.array([np.pad(frame[col].values, pad_width=(0, maxBlocks-len(frame)), mode='constant', constant_values=padVal) for _,frame in gb]).reshape(-1, maxBlocks, 1)    
        covArr = np.dstack((covArr,newArr))
    return covArr
    

#################################################### FUNCTION TO CONSTRUCT WINDOW-SPECIFIC CONTEXT VECTORS #################################################### 
#Function to construct window-specific context vectors for longitudinal covariate 'cov' at decision point time 't' from trained LSTM autoencoder
def context_vector_construct_load(longCov, longName, cov, t, model_path):   
    #Construct context vector using fitted LSTM autoencoder 
    longCov_lstm = longCov[:,:,np.newaxis]
    contextVec_model = tf.keras.models.load_model('{}/model_t{}_cov{}'.format(model_path, int(t), cov))
    cv = contextVec_model.predict(longCov_lstm)[:,0,:]  
    del contextVec_model
    
    #Save variable names 
    contVecNames = list()
    for j in np.arange(1, cv_dimension+1, 1):
        contVecNames.append('{}_{}'.format(longName,j)) 
    
    #Return window-specific context vector & associated variable names
    return (cv, contVecNames)      
    
    
#################################################### FUNCTION TO RETURN CONTEXT VECTORS FOR GIVEN DATAFRAME ####################################################
def create_valid_contVecs(df, t, model_path, minVec, maxVec):    
    #Min-max scale longitudinal covariates
    for i in np.arange(0,len(minVec)):
        df[['long_cov{}'.format(i+1)]] = minMaxScale(df[['long_cov{}'.format(i+1)]].values, minVec[i], maxVec[i])
    
    #Keep only measurements taken prior to/at t 
    simDF_t = df[(df.time<=t)].copy()
    
    #Create id & baseline covariate dataframe 
    respDF = simDF_t[["id", "base_cov"]].copy()
    respDF = respDF.drop_duplicates(subset='id', keep='first')
    respDF.index = range(0,len(respDF))
    
    #Create 3D covariate array (necessary for Keras)
    covDF = simDF_t[["id", "time", "long_cov1", "long_cov2", "long_cov3"]].copy()
    covArr = covArray_construct(covDF)
    longArr = np.delete(covArr,[0],2)
    longNames = list(covDF.columns)[2:]
    
    #Create lists to save context vectors & corresponding variable names
    context_vector_list = list()    
    contVec_name_list = list()   
    
    #Create context vector for each longitudinal covariate
    for i in np.arange(0,longArr.shape[2]):
        result = context_vector_construct_load(longArr[:,:,i], longNames[i], i, t, model_path)
        context_vector_list.append(result[0])
        contVec_name_list.append(result[1])
            
    #Combine baseline covariates, context vectors, and response variables
    regressMatrix = pd.concat([respDF, pd.DataFrame(np.hstack(context_vector_list), columns=np.hstack(contVec_name_list))], axis=1)
        
    #Return dataframe
    return(regressMatrix)
