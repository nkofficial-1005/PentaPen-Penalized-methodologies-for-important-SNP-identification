##############
# Clear environment
rm(list=ls())
gc()
options(warn=-1)

##############
# Import libraries
library(dplyr)
library(trio)
library(ROCR)
library(logr)
library(matrixStats)
library(Matrix)
library(glmnet)
library(SGL)
library(pROC)
library(dendextend)
library(InformationValue)
library(MLmetrics)
library(gglasso)
library(caTools)
library(caret)
library(tidyverse)
library(foreach)
library(doParallel)
library(grpreg)
##############
# Import pre-processed data
load("~/AtPolyDB Demo(Original) (PLINK)/PreprocessedData.RData")

##############
#Create 5 equally size folds
folds <- cut(seq(1,nrow(df_final)), breaks=5, labels=FALSE)

#############
#Define performance metrics
#Use appropriately for phenotype data type
#comment out in the loop accordingly
#1. Precision
precision = function(cm){
  diag(cm) / apply(cm, 2, sum)
}

#2. Recall
recall = function(cm){
  diag(cm) / apply(cm, 1, sum)
}

#3. RMSE
RMSE = function(y_actual,y_predict){
  sqrt(mean((y_actual-y_predict)^2))
}

#4. R SQUARED
RSQUARE = function(y_actual,y_predict){
  cor(y_actual,y_predict)^2
}

#############
#Perform 5 fold cross validation and parallel computation
for(i in 1:folds){
  
  print(paste("Fold:",i))
  #Segement your data by fold using the which() function 
  train_index <- sample(nrow(df_final), 0.8*nrow(df_final))
  train_data <- df_final[train_index,]
  val_data <- df_final[-train_index,]
  x_train <- train_data[,-ncol(train_data)]
  y_train  <- train_data[,ncol(train_data)]
  x_test <- val_data[,-ncol(val_data)]
  y_test <- val_data[,ncol(val_data)]
  nfolds <- 5
  
  source("BeforeSNPPool.R")
  #Parallel computing for Ridge, LASSO, Elastic net
  cl <- makeCluster(detectCores() - 1) # create a cluster with all but one core
  registerDoParallel(cl) # register the cluster with foreach
  ridge_lasso_elnet #Function to run Ridge, LASSO, and Elastic Net in parallel
  stopCluster(cl) # stop the cluster when done
  
  #Results from Parallel Computing of Ridge, LASSO, Elastic net
  #Ridge
  coef_ridge <- ridge_lasso_elnet[[1]]$coef[-1]#[length]
  coef_ridge <- abs(coef_ridge)
  coef_ridge <- sort(coef_ridge, decreasing = TRUE)
  selected_snps_ridge <- row.names(ridge_lasso_elnet[[1]]$coef)[-1]
  selected_nsnps_ridge <- length(coef_ridge[coef_ridge>=mean(coef_ridge)])
  selected_snps_ridge <- selected_snps_ridge[1:selected_nsnps_ridge]
  print(paste("Ridge snps",selected_nsnps_ridge))
  write(selected_snps_ridge, "ridgesnps.txt")
  
  #Train set
  #Predictions
  predictions_train_ridge <- ridge_lasso_elnet[[1]]$predictions_train_ridge
  #Precision
  # cm_ridge = as.matrix(table(Actual = y_train, Predicted = predictions_train_ridge))
  # precision_ridge_train <- precision(cm_ridge)
  # print(paste(mean(precision_ridge_train), "precision_ridge_train"))
  # 
  # #Recall
  # recall_ridge_train <- recall(cm_ridge)
  # print(paste(mean(recall_ridge_train), "recall_ridge_train"))
  # 
  # #F1-score
  # print(paste(F1_Score(predictions_train_ridge, y_train), "F1 train ridge"))
  # 
  # #AUC
  # auc_ridge_train <- auc(y_train, predictions_train_ridge)
  # print(paste(auc_ridge_train, "auc_ridge_train"))
  
  # #RMSE
  # rmse_ridge_train <- RMSE(y_train, predictions_train_ridge)
  # print(paste(rmse_ridge_train, "rmse_ridge_train"))
  # 
  # #R-squared
  # rsq_ridge_train <- RSQUARE(y_train, predictions_train_ridge)
  # print(paste(rsq_ridge_train, "rsq_ridge_train"))
  
  #Accuracy
  #print(paste(sum(y_train == as.integer(predictions_train_ridge)) / length(y_train), "Accuracy train ridge"))
  
  #Test set#
  #Predictions
  predictions_test_ridge <- ridge_lasso_elnet[[1]]$predictions_test_ridge
  
  #Precision
  # cm_ridge = as.matrix(table(Actual = y_test, Predicted = predictions_test_ridge))
  # precision_ridge_test <- precision(cm_ridge)
  # print(paste(mean(precision_ridge_test), "precision_ridge_test"))
  # 
  # #Recall
  # recall_ridge_test <- recall(cm_ridge)
  # print(paste(mean(recall_ridge_test), "recall_ridge_test"))
  # 
  # #F1-score
  # print(paste(F1_Score(predictions_test_ridge, y_test), "F1 test ridge"))
  # 
  # #AUC
  # auc_ridge_test <- auc(y_test, predictions_test_ridge)
  # print(paste(auc_ridge_test, "auc_ridge_test"))
  
  # #RMSE
  # rmse_ridge_test <- RMSE(y_test, predictions_test_ridge)
  # print(paste(rmse_ridge_test, "rmse_ridge_test"))
  # 
  # #R-squared
  # rsq_ridge_test <- RSQUARE(y_test, predictions_test_ridge)
  # print(paste(rsq_ridge_test, "rsq_ridge_test"))
  
  #Accuracy
  #print(paste(sum(y_test == as.integer(predictions_test_ridge)) / length(y_test), "Accuracy test ridge"))
  
  #Lasso
  coef_lasso <- ridge_lasso_elnet[[2]]$coef[-1]
  coef_lasso <- abs(coef_lasso)
  coef_lasso <- sort(coef_lasso, decreasing = TRUE)
  selected_snps_lasso <- row.names(ridge_lasso_elnet[[2]]$coef)[-1]
  selected_nsnps_lasso <- length(coef_lasso[coef_lasso>=mean(coef_lasso)])
  selected_snps_lasso <- selected_snps_lasso[1:selected_nsnps_lasso]
  print(paste("LASSO snps",selected_nsnps_lasso))
  write(selected_snps_lasso, "lassosnps.txt")
  
  #Train set#
  #Predictions
  predictions_train_lasso <- ridge_lasso_elnet[[2]]$predictions_train_lasso
  
  #Precision
  # cm_lasso = as.matrix(table(Actual = y_train, Predicted = predictions_train_lasso))
  # precision_lasso_train <- precision(cm_lasso)
  # print(paste(mean(precision_lasso_train), "precision_lasso_train"))
  # 
  # #Recall
  # recall_lasso_train <- recall(cm_lasso)
  # print(paste(mean(recall_lasso_train), "recall_lasso_train"))
  # 
  # #F1-score
  # print(paste(F1_Score(predictions_train_lasso, y_train), "F1 train lasso"))
  # 
  # #AUC
  # auc_lasso_train <- auc(y_train, predictions_train_lasso)
  # print(paste(auc_lasso_train, "auc_lasso_train"))
  
  # #RMSE
  # rmse_lasso_train <- RMSE(y_train, predictions_train_lasso)
  # print(paste(rmse_lasso_train, "rmse_lasso_train"))
  # 
  # #R-squared
  # rsq_lasso_train <- RSQUARE(y_train, predictions_train_lasso)
  # print(paste(rsq_lasso_train, "rsq_lasso_train"))
  
  #Accuracy
  #print(paste(sum(y_train == as.integer(predictions_train_lasso)) / length(y_train), "Accuracy train lasso"))
  
  #Test set#
  #Predictions
  predictions_test_lasso <- ridge_lasso_elnet[[2]]$predictions_test_lasso
  
  #Precision
  # cm_lasso = as.matrix(table(Actual = y_test, Predicted = predictions_test_lasso))
  # precision_lasso_test <- precision(cm_lasso)
  # print(paste(mean(precision_lasso_test), "precision_lasso_test"))
  # 
  # #Recall
  # recall_lasso_test <- recall(cm_lasso)
  # print(paste(mean(recall_lasso_test), "recall_lasso_test"))
  # 
  # #F1-score
  # print(paste(F1_Score(predictions_test_lasso, y_test), "F1 test lasso"))
  # 
  # #AUC
  # auc_lasso_test <- auc(y_test, predictions_test_lasso)
  # print(paste(auc_lasso_test, "auc_lasso_test"))
  
  # #RMSE
  # rmse_lasso_test <- RMSE(y_test, predictions_test_lasso)
  # print(paste(rmse_lasso_test, "rmse_lasso_test"))
  # 
  # #R-squared
  # rsq_lasso_test <- RSQUARE(y_test, predictions_test_lasso)
  # print(paste(rsq_lasso_test, "rsq_lasso_test"))
  
  #Accuracy
  #print(paste(sum(y_test == as.integer(predictions_test_lasso)) / length(y_test), "Accuracy test lasso"))
  ########
  
  #Elastic net
  coef_enet <- ridge_lasso_elnet[[3]]$coef[-1]
  coef_enet <- abs(coef_enet)
  coef_enet <- sort(coef_enet, decreasing = TRUE)
  selected_snps_enet <- row.names(ridge_lasso_elnet[[3]]$coef)[-1]
  selected_nsnps_enet <- length(coef_enet[coef_enet>=mean(coef_enet)])
  selected_snps_enet <- selected_snps_enet[1:selected_nsnps_enet]
  print(paste("Elastic net snps",selected_nsnps_enet))
  write(selected_snps_enet, "enetsnps.txt")
  
  #Train set#
  #Predictions
  predictions_train_elnet <- ridge_lasso_elnet[[3]]$predictions_train_enet
  
  #Precision
  # cm_elnet = as.matrix(table(Actual = y_train, Predicted = predictions_train_elnet))
  # precision_elnet_train <- precision(cm_elnet)
  # print(paste(mean(precision_elnet_train), "precision_elnet_train"))
  # 
  # #Recall
  # recall_elnet_train <- recall(cm_elnet)
  # print(paste(mean(recall_elnet_train), "recall_elnet_train"))
  # 
  # #F1-score
  # print(paste(F1_Score(predictions_train_elnet, y_train), "F1 train elnet"))
  # 
  # #AUC
  # auc_elnet_train <- auc(y_train, predictions_train_elnet)
  # print(paste(auc_elnet_train, "auc_elnet_train"))
  
  # #RMSE
  # rmse_elnet_train <- RMSE(y_train, predictions_train_elnet)
  # print(paste(rmse_elnet_train, "rmse_elnet_train"))
  # 
  # #R-squared
  # rsq_elnet_train <- RSQUARE(y_train, predictions_train_elnet)
  # print(paste(rsq_elnet_train, "rsq_elnet_train"))
  
  #Accuracy
  #print(paste(sum(y_train == as.integer(predictions_train_elnet)) / length(y_train), "Accuracy train elnet"))
  
  #Test set#
  #Predictions
  predictions_test_elnet <- ridge_lasso_elnet[[3]]$predictions_test_enet
  
  #Precision
  # cm_elnet = as.matrix(table(Actual = y_test, Predicted = predictions_test_elnet))
  # precision_elnet_test <- precision(cm_elnet)
  # print(paste(mean(precision_elnet_test), "precision_elnet_test"))
  # 
  # #Recall
  # recall_elnet_test <- recall(cm_elnet)
  # print(paste(mean(recall_elnet_test), "recall_elnet_test"))
  # 
  # #F1-score
  # print(paste(F1_Score(predictions_test_elnet, y_test), "F1 test elnet"))
  # 
  # #AUC
  # auc_elnet_test <- auc(y_test, predictions_test_elnet)
  # print(paste(auc_elnet_test, "auc_elnet_test"))
  
  # #RMSE
  # rmse_elnet_test <- RMSE(y_test, predictions_test_elnet)
  # print(paste(rmse_elnet_test, "rmse_elnet_test"))
  # 
  # #R-squared
  # rsq_elnet_test <- RSQUARE(y_test, predictions_test_elnet)
  # print(paste(rsq_elnet_test, "rsq_elnet_test"))
  
  #Accuracy
  #print(paste(sum(y_test == as.integer(predictions_test_elnet)) / length(y_train), "Accuracy test elnet"))
  ########
  
  #SNP Pooling
  selected_snps_union <- c(selected_snps_ridge, selected_snps_lasso, selected_snps_enet)
  selected_snps_union <- unique(selected_snps_union)
  #selected_snps_union <- sample(1:length(selected_snps_union), 100) #Comment this out if you want to see the reproducibility of the code in <4GB RAM
  print(paste("Snp Pool",length(selected_snps_union)))
  ##############
  
  #Data Preparation of next computation of group lasso and SGL
  data_subset <- train_data[,selected_snps_union]
  Z.train <- train_data[,selected_snps_union]
  y.train <- train_data[,ncol(train_data)]
  Z.test <- val_data[,selected_snps_union]
  y.test <- val_data[,ncol(val_data)]
  
  snps_df_sc <- as.data.frame(scale(t(Z.train)))
  dist_mat <- dist(snps_df_sc, method = 'euclidean')
  hclust_avg <- hclust(dist_mat, method = 'average')
  cut_avg <- cutree(hclust_avg, k = 10)
  
  Z.train.sd <- colSds(Z.train)
  Z.test.sd <- colSds(Z.test)
  Z.train    <- Z.train[, which(Z.train.sd > 1e-7)]
  Z.test     <- Z.test[, which(Z.test.sd > 1e-7)]
  grp <- cut_avg[which(Z.train.sd > 1e-7)]
  Z.train.sd <- Z.train.sd[which(Z.train.sd > 1e-7)]
  Z.test.sd <- Z.test.sd[which(Z.test.sd  > 1e-7)]
  ########
  
  #Data for SGL
  Z.train_sgl <- scale(Z.train, center = TRUE, scale = Z.train.sd)
  Z.test_sgl <- scale(Z.test, center = TRUE, scale = Z.test.sd)
  data.SGL <- list(x = Z.train_sgl, y = y.train)
  
  #Parallel computing of Group LASSO and SGL
  cl <- makeCluster(detectCores() - 1) # create a cluster with all but one core
  registerDoParallel(cl) # register the cluster with foreach
  source("AfterSNPPool.R")
  grp_sgl #Function to run Group LASSO and SGL in parallel
  stopCluster(cl) # stop the cluster when done
  
  #Results from parallel computation of Group LASSO and SGL
  #Group lasso
  grplasso.coef <- grp_sgl[[1]]$coef
  names(grplasso.coef) <- colnames(Z.train)
  grplasso.coef <- abs(grplasso.coef)
  grplasso.coef <- sort(grplasso.coef, decreasing = TRUE)
  grplasso.snps <- names(grplasso.coef)
  grplasso.snps <- grplasso.snps[grplasso.snps != ""]
  grplasso.nsnps <- length(grplasso.coef[grplasso.coef >= mean(grplasso.coef)])
  if(grplasso.nsnps > ncol(train_data[,-ncol(train_data)])){
    grplasso.nsnps <- ncol(train_data[,-ncol(train_data)])
  }
  grplasso.snps <- grplasso.snps[1:grplasso.nsnps]
  print(paste("Group Lasso snps",grplasso.nsnps))
  write(grplasso.snps, "grplassosnps.txt")
  
  #Using training set#
  #Predictions
  predictions_train_grplasso <- grp_sgl[[1]]$predictions_train_grplasso
  
  #Precision
  # cm_grplasso = as.matrix(table(Actual = y.train, Predicted = predictions_train_grplasso))
  # precision_grplasso_train <- precision(cm_grplasso)
  # print(paste(mean(precision_grplasso_train), "precision_grplasso_train"))
  # 
  # #Recall
  # recall_grplasso_train <- recall(cm_grplasso)
  # print(paste(mean(recall_grplasso_train), "recall_grplasso_train"))
  # 
  # #F1-score
  # print(paste(F1_Score(predictions_train_grplasso, y.train), "F1 train grplasso"))
  # 
  # #AUC
  # auc_grplasso_train <- auc(y.train, predictions_train_grplasso)
  # print(paste(auc_grplasso_train, "auc_grplasso_train"))
  
  # #RMSE
  # rmse_grplasso_train <- RMSE(y.train, predictions_train_grplasso)
  # print(paste(rmse_grplasso_train, "rmse_grplasso_train"))
  # 
  # #R-squared
  # rsq_grplasso_train <- RSQUARE(y.train, predictions_train_grplasso)
  # print(paste(rsq_grplasso_train, "rsq_grplasso_train"))
  
  #Accuracy
  #print(paste(sum(y.train == as.integer(predictions_train_grplasso)) / length(y.train), "Accuracy train grplasso"))
  #Using testing set#
  #Predictions
  predictions_test_grplasso <- grp_sgl[[1]]$predictions_test_grplasso
  
  #Precision
  # cm_grplasso = as.matrix(table(Actual = y.test, Predicted = predictions_test_grplasso))
  # precision_grplasso_test <- precision(cm_grplasso)
  # print(paste(mean(precision_grplasso_test), "precision_grplasso_test"))
  # 
  # #Recall
  # recall_grplasso_test <- recall(cm_grplasso)
  # print(paste(mean(recall_grplasso_test), "recall_grplasso_test"))
  # 
  # #F1-score
  # print(paste(F1_Score(predictions_test_grplasso, y.test), "F1 test grplasso"))
  # 
  # #AUC
  # auc_grplasso_test <- auc(y.test, predictions_test_grplasso)
  # print(paste(auc_grplasso_test, "auc_grplasso_test"))
  
  # #RMSE
  # rmse_grplasso_test <- RMSE(y.test, predictions_test_grplasso)
  # print(paste(rmse_grplasso_test, "rmse_grplasso_test"))
  # 
  # #R-squared
  # rsq_grplasso_test <- RSQUARE(y.test, predictions_test_grplasso)
  # print(paste(rsq_grplasso_test, "rsq_grplasso_test"))
  
  #Accuracy
  #print(paste(sum(y.test == as.integer(predictions_test_grplasso)) / length(y.test), "Accuracy test grplasso"))
  
  #SGL
  beta <- grp_sgl[[2]]$coef
  sgl.coef <- data.frame(beta)
  rownames(sgl.coef) <- colnames(Z.test)
  sgl.coef <- abs(sgl.coef)
  #sgl.coef <- sort(sgl.coef, decreasing = TRUE)
  sgl.snps <- rownames(sgl.coef)[apply(sgl.coef, 1, function(u) any(u>0.5))]
  print(paste("SGL snps",length(sgl.snps)))
  write(sgl.snps, "sglsnps.txt")
  
  #Using training set#
  #Predictions
  sgl.probabilities.train <- grp_sgl[[2]]$sgl.probabilities.train
  
  #Precision
  # cm_sgl = as.matrix(table(Actual = y.train, Predicted = sgl.probabilities.train))
  # precision_sgl_train <- precision(cm_sgl)
  # print(paste(mean(precision_sgl_train), "precision_sgl_train"))
  # 
  # #Recall
  # recall_sgl_train <- recall(cm_sgl)
  # print(paste(mean(recall_sgl_train), "recall_sgl_train"))
  # 
  # #F1-score
  # print(paste(F1_Score(y.train, sgl.probabilities.train), "F1 train sgl"))
  # 
  # #AUC
  # auc_sgl_train <- auc((y.train), sgl.probabilities.train)
  # print(paste(auc_sgl_train, "auc_sgl_train"))
  
  # #RMSE
  # rmse_sgl_train <- RMSE(y.train, sgl.probabilities.train)
  # print(paste(rmse_sgl_train, "rmse_sgl_train"))
  # 
  # #R-squared
  # rsq_sgl_train <- RSQUARE(y.train, sgl.probabilities.train)
  # print(paste(rsq_sgl_train, "rsq_sgl_train"))
  
  #Accuracy
  #print(paste(sum(y.train == as.integer(sgl.probabilities.train)) / length(y.train), "Accuracy train sgl"))
  
  #Using testing set#
  #Predictions
  sgl.probabilities.test <- grp_sgl[[2]]$sgl.probabilities.test
  
  #Precision
  # cm_sgl = as.matrix(table(Actual = y.test, Predicted = sgl.probabilities.test))
  # precision_sgl_test <- precision(cm_sgl)
  # print(paste(mean(precision_sgl_test), "precision_sgl_test"))
  # 
  # #Recall
  # recall_sgl_test <- recall(cm_sgl)
  # print(paste(mean(recall_sgl_test), "recall_sgl_test"))
  # 
  # #F1-score
  # print(paste(F1_Score(y.test, sgl.probabilities.test), "F1 test sgl"))
  # 
  # #AUC
  # auc_sgl_test <- auc((y.test), sgl.probabilities.test)
  # print(paste(auc_sgl_test, "auc_sgl_test"))
  
  # #RMSE
  # rmse_sgl_test <- RMSE(y.test, sgl.probabilities.test)
  # print(paste(rmse_sgl_test, "rmse_sgl_test"))
  # 
  # #R-squared
  # rsq_sgl_test <- RSQUARE(y.test, sgl.probabilities.test)
  # print(paste(rsq_sgl_test, "rsq_sgl_test"))
  
  #Accuracy
  #print(paste(sum(y.test == as.integer(sgl.probabilities.test)) / length(y.test), "Accuracy test sgl"))
  ########
  
  #Final SNPs from the proposed algorithm
  final_snps <- c(grplasso.snps, sgl.snps)
  print(paste("Algorithm snps",length(final_snps)))
  write(final_snps, "algosnps.txt")
}

# create aggregated model by averaging predictions from individual models
agg_model_predictions_train <- rowMeans(do.call(cbind, list(ridge_lasso_elnet[[1]]$predictions_train_ridge, 
                                                            ridge_lasso_elnet[[2]]$predictions_train_lasso, 
                                                            ridge_lasso_elnet[[3]]$predictions_train_enet,
                                                            grp_sgl[[1]]$predictions_train_grplasso,
                                                            grp_sgl[[2]]$sgl.probabilities.train)))

agg_model_predictions_test <- rowMeans(do.call(cbind, list(ridge_lasso_elnet[[1]]$predictions_test_ridge, 
                                                           ridge_lasso_elnet[[2]]$predictions_test_lasso, 
                                                           ridge_lasso_elnet[[3]]$predictions_test_enet,
                                                           grp_sgl[[1]]$predictions_test_grplasso,
                                                           grp_sgl[[2]]$sgl.probabilities.test)))

# calculate performance metric of aggregated model
#Train set
#Precision
# cm_agg_model = as.matrix(table(Actual = y_train, Predicted = agg_model_predictions_train))
# precision_agg_model_train <- precision(cm_agg_model)
# print(paste(mean(precision_agg_model_train), "precision_agg_model_train"))
# 
#Recall
# recall_agg_model_train <- recall(cm_agg_model)
# print(paste(mean(recall_agg_model_train), "recall_agg_model_train"))
# 
# #F1-score
# print(paste(F1_Score(y_train, agg_model_predictions_train), "F1 train agg_model"))
# 
# #AUC
# auc_agg_model_train <- auc((y_train), agg_model_predictions_train)
# print(paste(auc_agg_model_train, "auc_agg_model_train"))

# #RMSE
# rmse_agg_model_train <- RMSE(y_train, agg_model_predictions_train)
# print(paste(rmse_agg_model_train, "rmse_agg_model_train"))
# 
# #R-squared
# rsq_agg_model_train <- RSQUARE(y_train, agg_model_predictions_train)
# print(paste(rsq_agg_model_train, "rsq_agg_model_train"))

#Accuracy
#print(paste(sum(y_train == as.integer(agg_model_predictions_train)) / length(y_train), "Accuracy train agg_model"))
########

#Test set
#Precision
# cm_agg_model = as.matrix(table(Actual = y_test, Predicted = agg_model_predictions_test))
# precision_agg_model_test <- precision(cm_agg_model)
# print(paste(mean(precision_agg_model_test), "precision_agg_model_test"))
# 
# #Recall
# recall_agg_model_test <- recall(cm_agg_model)
# print(paste(mean(recall_agg_model_test), "recall_agg_model_test"))
# 
# #F1-score
# print(paste(F1_Score(y_test, agg_model_predictions_test), "F1 test agg_model"))
# 
# #AUC
# auc_agg_model_test <- auc((y_test), agg_model_predictions_test)
# print(paste(auc_agg_model_test, "auc_agg_model_test"))

# #RMSE
# rmse_agg_model_test <- RMSE(y_test, agg_model_predictions_test)
# print(paste(rmse_agg_model_test, "rmse_agg_model_test"))
# 
# #R-squared
# rsq_agg_model_test <- RSQUARE(y_test, agg_model_predictions_test)
# print(paste(rsq_agg_model_test, "rsq_agg_model_test"))

#Accuracy
print(paste(sum(y_test == as.integer(agg_model_predictions_test)) / length(y_test), "Accuracy test agg_model"))
########
