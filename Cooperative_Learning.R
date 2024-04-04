# Cooperative learning is used for the supervised learning task to build prediction models for the outcome with multimodal data.

library(multiview)
library(glmnet)
library (janitor)
library(dplyr)
library(tidyr)
library("data.table")
library('caret')

nested_cv_dual <- function(x,z,y){   # Assuming we have two data matrices, x and z, on the same set of n observations; y is the outcome
  nested_k_folds <- 5                # The final performance scores are averaged after a 5-fold nested CV
  data_tr <- data.table(x, z, y, fold = 0)
  set.seed()
  for (y_i in levels(y) ) {                                          # Create stratified folds
    nrow_i       <- nrow(data_tr[y == y_i,])
    n_per_fold_i <- ceiling(nrow_i / nested_k_folds)
    data_tr[y == y_i, fold := sample(rep(1:nested_k_folds, n_per_fold_i), nrow_i, replace = FALSE)]
  }
  data_tr$y = as.factor(data_tr$y)
  
  fold_results <- list()
  
  for (i in 1:nested_k_folds){
    print(paste('fold:',i))
    fold_result <- list()
    training <- data_tr[data_tr$fold != i,]
    testing <- data_tr[data_tr$fold == i,]
    y_train <- as.factor(training[,y])
    X_train <- apply(as.matrix(training[,1:ncol(x)]),2,as.numeric)
    Z_train <- apply(as.matrix(training[,(ncol(x)+1):(ncol(x)+ncol(z))]),2,as.numeric)
    y_test <- as.factor(testing[,y])
    X_test <- apply(as.matrix(testing[,1:ncol(x)]),2,as.numeric)
    Z_test <- apply(as.matrix(testing[,(ncol(x)+1):(ncol(x)+ncol(z))]),2,as.numeric)
    
    alphalist <- seq(0,1,by=0.2)      # alpha is for the elastic net mixing parameter
    rholist <- c(0,0.1,0.25,0.5,1)    # rho is the aggrement penalty to encourage alignment between predictions from different data modalities
    for (m in alphalist){
      for (n in rholist){
        hyper_param <- c(m,n)
        set.seed()
        cvfit <- cv.multiview(list(X_train,Z_train), y_train,family = binomial(),alpha=m,rho=n,nfolds=10,type.measure = 'deviance',keep = TRUE)
        cvm_min_index <- cvfit$index["min",]
        rocs_cv <- roc.glmnet(cvfit$fit.preval, newy = y_train)
        assess_cv <- assess.glmnet(cvfit$fit.preval, newy = y_train, family = "binomial")
        cvm_cv_min <- min(cvfit$cvm)
        auc_cv_max <- max(assess_cv$auc)
        auc_cv_max_index <- which.max(assess_cv$auc)
        
        pred_deviance <- predict(cvfit,newx = list(X_test,Z_test),s='lambda.1se')
        pred_auc <- predict(cvfit,newx = list(X_test,Z_test),s=cvfit$lambda[auc_cv_max_index])
        
        assess_pred_deviance <- assess.glmnet(pred_deviance, newy = y_test, family = 'binomial')
        assess_pred_auc <- assess.glmnet(pred_auc, newy = y_test, family = 'binomial')
        pred_deviance_auc <- assess_pred_deviance$auc
        pred_auc_auc <- assess_pred_auc$auc
        
        conf_pred_deviance <- confusion.glmnet(pred_deviance, newy = y_test, family = "binomial")
        conf_pred_auc <- confusion.glmnet(pred_auc, newy = y_test, family = "binomial")
        
        roc_pred_deviance <- roc.glmnet(pred_deviance, newy = y_test)
        roc_pred_auc <- roc.glmnet(pred_auc, newy = y_test)
        
        fold_result[[length(fold_result) + 1]] <- list(hyper_param,cvfit,rocs_cv,assess_cv,cvm_cv_min,auc_cv_max,pred_deviance,pred_auc,
                                                       assess_pred_deviance,assess_pred_auc,pred_deviance_auc,pred_auc_auc,
                                                       conf_pred_deviance,conf_pred_auc,roc_pred_deviance,roc_pred_auc)
      }
    }
    fold_results[[length(fold_results) + 1]] <- fold_result
  }
  return(fold_results) # pred_deviance_auc and conf_pred_deviance are the measured metrics in the paper. One can evaluate other metrics if needed


nested_cv_tetra <- function(x,z,a,b,y){   # Assuming we have four data matrices, x, z, a and b on the same set of n observations; y is the outcome
  nested_k_folds <- 5                     # The final performance scores are averaged after a 5-fold nested CV
  data_tr <- data.table(x, z, a,b,y, fold = 0)
  set.seed()
  for (y_i in levels(y) ) {                                          # Create stratified folds
    nrow_i       <- nrow(data_tr[y == y_i,])
    n_per_fold_i <- ceiling(nrow_i / nested_k_folds)
    data_tr[y == y_i, fold := sample(rep(1:nested_k_folds, n_per_fold_i), nrow_i, replace = FALSE)]
  }
  data_tr$y = as.factor(data_tr$y)
  
  fold_results <- list()
  
  for (i in 1:nested_k_folds){
    print(paste('fold:',i))
    fold_result <- list()
    training <- data_tr[data_tr$fold != i,]
    testing <- data_tr[data_tr$fold == i,]
    y_train <- as.factor(training[,y])
    X_train <- apply(as.matrix(training[,1:ncol(x)]),2,as.numeric)
    Z_train <- apply(as.matrix(training[,(ncol(x)+1):(ncol(x)+ncol(z))]),2,as.numeric)
    A_train <- apply(as.matrix(training[,(ncol(x)+ncol(z)+1):(ncol(x)+ncol(z)+ncol(a))]),2,as.numeric)
    B_train <- apply(as.matrix(training[,(ncol(x)+ncol(z)+ncol(a)+1):(ncol(x)+ncol(z)+ncol(a)+ncol(b))]),2,as.numeric)
    y_test <- as.factor(testing[,y])
    X_test <- apply(as.matrix(testing[,1:ncol(x)]),2,as.numeric)
    Z_test <- apply(as.matrix(testing[,(ncol(x)+1):(ncol(x)+ncol(z))]),2,as.numeric)
    A_test <- apply(as.matrix(testing[,(ncol(x)+ncol(z)+1):(ncol(x)+ncol(z)+ncol(a))]),2,as.numeric)
    B_test <- apply(as.matrix(testing[,(ncol(x)+ncol(z)+ncol(a)+1):(ncol(x)+ncol(z)+ncol(a)+ncol(b))]),2,as.numeric)
    
    alphalist <- seq(0,1,by=0.2)             # alpha is for the elastic net mixing parameter
    rholist <- c(0,0.1,0.25,0.5,1)           # rho is the aggrement penalty to encourage alignment between predictions from different data modalities
    for (m in alphalist){
      for (n in rholist){
        hyper_param <- c(m,n)
        set.seed()
        cvfit <- cv.multiview(list(X_train,Z_train,A_train,B_train), y_train,family = binomial(),alpha=m,rho=n,nfolds=10,type.measure = 'deviance',keep = TRUE)
        cvm_min_index <- cvfit$index["min",]
        rocs_cv <- roc.glmnet(cvfit$fit.preval, newy = y_train)
        assess_cv <- assess.glmnet(cvfit$fit.preval, newy = y_train, family = "binomial")
        cvm_cv_min <- min(cvfit$cvm)
        auc_cv_max <- max(assess_cv$auc)
        auc_cv_max_index <- which.max(assess_cv$auc)
        
        pred_deviance <- predict(cvfit,newx = list(X_test,Z_test,A_test,B_test),s='lambda.1se')
        pred_auc <- predict(cvfit,newx = list(X_test,Z_test,A_test,B_test),s=cvfit$lambda[auc_cv_max_index])
        
        assess_pred_deviance <- assess.glmnet(pred_deviance, newy = y_test, family = 'binomial')
        assess_pred_auc <- assess.glmnet(pred_auc, newy = y_test, family = 'binomial')
        pred_deviance_auc <- assess_pred_deviance$auc
        pred_auc_auc <- assess_pred_auc$auc
        
        conf_pred_deviance <- confusion.glmnet(pred_deviance, newy = y_test, family = "binomial")
        conf_pred_auc <- confusion.glmnet(pred_auc, newy = y_test, family = "binomial")
        
        roc_pred_deviance <- roc.glmnet(pred_deviance, newy = y_test)
        roc_pred_auc <- roc.glmnet(pred_auc, newy = y_test)
        
        fold_result[[length(fold_result) + 1]] <- list(hyper_param,cvfit,rocs_cv,assess_cv,cvm_cv_min,auc_cv_max,pred_deviance,pred_auc,
                                                       assess_pred_deviance,assess_pred_auc,pred_deviance_auc,pred_auc_auc,
                                                       conf_pred_deviance,conf_pred_auc,roc_pred_deviance,roc_pred_auc)
      }
    }
    fold_results[[length(fold_results) + 1]] <- fold_result 
  }
  return(fold_results) # pred_deviance_auc and conf_pred_deviance are the measured metrics in the paper. One can evaluate other metrics if needed
}
