##############
# Clear environment
rm(list=ls())
gc()
options(warn=-1)

##############
# Import libraries
library(dplyr)
library(ROCR)
library(logr)
library(matrixStats)
library(Matrix)
library(glmnet)
library(SGL)
library(pROC)
library(MLmetrics)
library(gglasso)
library(caTools)
library(caret)
library(tidyverse)
##############
# Import data
Geno <- read.pedfile("genotype.ped") 
char.pheno <- read.table("phenotypes.pheno", header = TRUE, stringsAsFactors = FALSE, sep = " ")

##############
# Data Pre-processing
set.seed(100)
# Re-code the data in ped file
Geno[Geno == 'A'] <- 0  # Converting A to 0
Geno[Geno == 'T'] <- 1  # Converting T to 1
Geno[Geno == 'G'] <- 2  # Converting G to 2
Geno[Geno == 'C'] <- 3  # Converting C to 3

#Convert the phenotype to matrix
y <- matrix(char.pheno$Anthocyanin_22) #Change the phenotype accordingly
rownames(y) <- char.pheno$IID
index <- !is.na(y)
y <- y[index, 1, drop = FALSE]

#Imputing SNP null values
for (j in 1:ncol(Geno)) {
  Geno[, j] <- ifelse(is.na(Geno[, j]), mean(Geno[, j], na.rm = TRUE), Geno[, 
                                                                            j])
}
Geno_y <- Geno[index, ] 
pheno_final <- data.frame(famid = rownames(y), y = y)

df <- merge(Geno_y, pheno_final, by = 'famid')
df_final <- df[, 7:204760] #comment here what is done
df_final <- sapply(df_final, as.numeric)
df_final <- df_final[sample(nrow(df_final)),]
n <- dim(df_final)[[1]] 
d <- dim(df_final)[[2]] 

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
#Perform 5 fold cross validation
for(i in unique(folds)){
  
  print(paste("Fold:",i))
  
  #Segement the data by fold using the which() function 
  testIndexes <- which(folds==i, arr.ind=TRUE)
  testData <- df_final[testIndexes, ]
  trainData <- df_final[-testIndexes, ]
  data.train <- as.data.frame(trainData)
  data.test <- as.data.frame(testData)
  x_test <- as.matrix(data.test[,1:(d-1)])
  y_test <- as.matrix(data.test[, d])
  x_train <- as.matrix(data.train[,1:(d-1)])
  y_train <- as.matrix(data.train[, d])
  
  # alpha = 0, Ridge
  ##############
  set.seed(100)
  ridge.cv <- cv.glmnet(x_train, y_train, type.measure="auc", alpha=0, nfolds = 5)
  ridge.fit <- glmnet(x_train, y_train, type.measure="auc"
                      , alpha=0, lambda = ridge.cv$lambda.1se)
  # Sort and Extract Coefficients and corresponding SNPs
  ridge.coef <- ridge.fit$beta[, which(ridge.fit$lambda==ridge.cv$lambda.1se)]
  ridge.coef <- abs(ridge.coef)
  ridge.coef <- sort(ridge.coef, decreasing = TRUE)
  ridge.snps <- names(ridge.coef)
  ridge.nsnps <- length(ridge.coef[ridge.coef>=mean(ridge.coef)])
  ridge.snps <- ridge.snps[1:ridge.nsnps]
  #Export the SNPs extracted from Ridge
  write(ridge.snps, "ridgesnps.txt")
  
  #USing training set#
  #Predictions
  predictions_train_ridge <- predict(ridge.fit, newx = x_train)

  #Precision
  cm_ridge = as.matrix(table(Actual = y_train, Predicted = predictions_train_ridge))
  precision_ridge_train <- precision(cm_ridge)
  print(paste(mean(precision_ridge_train), "precision_ridge_train"))
  
  #Recall
  recall_ridge_train <- recall(cm_ridge)
  print(paste(mean(recall_ridge_train), "recall_ridge_train"))
  
  #F1-score
  print(paste(F1_Score(predictions_train_ridge, y_train), "F1 train ridge"))

  #AUC
  auc_ridge_train <- auc(y_train, predictions_train_ridge)
  print(paste(auc_ridge_train, "auc_ridge_train"))
  
  #RMSE
  rmse_ridge_train <- RMSE(y_train, predictions_train_ridge)
  print(paste(rmse_ridge_train, "rmse_ridge_train"))
  
  #R-squared
  rsq_ridge_train <- RSQUARE(y_train, predictions_train_ridge)
  print(paste(rsq_ridge_train, "rsq_ridge_train"))
  
  #Accuracy
  print(paste(accuracy(predictions_train_ridge, y_train), "Accuracy train ridge"))
  
  #USing testing set#
  #Predictions
  predictions_test_ridge <- predict(ridge.fit, s = ridge.cv$lambda.1se, newx = x_test)

  #Precision
  cm_ridge = as.matrix(table(Actual = y_test, Predicted = predictions_test_ridge))
  precision_ridge_test <- precision(cm_ridge)
  print(paste(mean(precision_ridge_test), "precision_ridge_test"))
  
  #Recall
  recall_ridge_test <- recall(cm_ridge)
  print(paste(mean(recall_ridge_test), "recall_ridge_test"))
  
  #F1-score
  print(paste(F1_Score(predictions_test_ridge, y_test), "F1 test lasso"))

  #AUC
  auc_ridge_test <- auc(y_test, predictions_test_ridge)
  print(paste(auc_ridge_test, "auc_ridge_test"))
  
  #RMSE
  rmse_ridge_test <- RMSE(y_test, predictions_test_ridge)
  print(paste(rmse_ridge_test, "rmse_ridge_test"))
  
  #R-squared
  rsq_ridge_test <- RSQUARE(y_test, predictions_test_ridge)
  print(paste(rsq_ridge_test, "rsq_ridge_test"))
  
  #Accuracy
  print(paste(accuracy(predictions_test_ridge, y_test), "Accuracy test lasso"))
  
  ########

  # alpha = 1, Lasso
  ##############
  set.seed(100)
  lasso.cv <- cv.glmnet(x_train, y_train, type.measure="auc", alpha=1, nfolds = 5)
  lasso.fit <- glmnet(x_train, y_train, type.measure="auc"
                      , alpha=1, lambda = lasso.cv$lambda.1se)
  # Sort and Extract Coefficients and corresponding SNPs
  lasso.coef <- lasso.fit$beta[, which(lasso.fit$lambda==lasso.cv$lambda.1se)]
  lasso.coef <- abs(lasso.coef)
  lasso.coef <- sort(lasso.coef, decreasing = TRUE)
  lasso.snps <- names(lasso.coef)
  lasso.nsnps <- length(lasso.coef[lasso.coef>=mean(lasso.coef)])
  lasso.snps <- lasso.snps[1:lasso.nsnps]
  #Export the SNPs extracted from LASSO
  write(lasso.snps, "lassosnps.txt")
  
  #Using training set#
  #Predictions
  predictions_train_lasso <- predict(lasso.fit, s = lasso.cv$lambda.1se, newx = x_train)

  #Precision
  cm_lasso = as.matrix(table(Actual = y_train, Predicted = predictions_train_lasso))
  precision_lasso_train <- precision(cm_lasso)
  print(paste(mean(precision_lasso_train), "precision_lasso_train"))
  
  #Recall
  recall_lasso_train <- recall(cm_lasso)
  print(paste(mean(recall_lasso_train), "recall_lasso_train"))
  
  #F1-score
  print(paste(F1_Score(predictions_train_lasso, y_train), "F1 train lasso"))

  #AUC
  auc_lasso_train <- auc(y_train, predictions_train_lasso)
  print(paste(auc_lasso_train, "auc_lasso_train"))
  
  #RMSE
  rmse_lasso_train <- RMSE(y_train, predictions_train_lasso)
  print(paste(rmse_lasso_train, "rmse_lasso_train"))
  
  #R-squared
  rsq_lasso_train <- RSQUARE(y_train, predictions_train_lasso)
  print(paste(rsq_lasso_train, "rsq_lasso_train"))
  
  #Accuracy
  print(paste(accuracy(predictions_train_lasso, y_train), "Accuracy train lasso"))
  
  #Using testing set#
  #Predictions
  predictions_test_lasso <- predict(lasso.fit, s = lasso.cv$lambda.1se, newx = x_test)

  #Precision 
  cm_lasso = as.matrix(table(Actual = y_test, Predicted = predictions_test_lasso))
  precision_lasso_test <- precision(cm_lasso)
  print(paste(mean(precision_lasso_test), "precision_lasso_test"))
  
  #Recall
  recall_lasso_test <- recall(cm_lasso)
  print(paste(mean(recall_lasso_test), "recall_lasso_test"))
  
  #F1-score
  print(paste(F1_Score(predictions_test_lasso, y_test), "F1 test lasso"))

  #AUC
  auc_lasso_test <- auc(y_test, predictions_test_lasso)
  print(paste(auc_lasso_test, "auc_lasso_test"))
  
  #RMSE
  rmse_lasso_test <- RMSE(y_test, predictions_test_lasso)
  print(paste(rmse_lasso_test, "rmse_lasso_test"))
  
  #R-squared
  rsq_lasso_test <- RSQUARE(y_test, predictions_test_lasso)
  print(paste(rsq_lasso_test, "rsq_lasso_test"))
  
  #Accuracy
  print(paste(accuracy(predictions_test_lasso, y_test), "Accuracy test lasso"))
  ########

  # alpha = 0.5, Elasticnet
  ##############
  set.seed(100)
  elnet.cv <- cv.glmnet(x_train, y_train, type.measure="auc", alpha=0.5, nfolds = 5)
  elnet.fit <- glmnet(x_train, y_train, type.measure="auc"
                      , alpha=0.5, lambda = elnet.cv$lambda.1se)
  # Sort and Extract Coefficients and corresponding SNPs
  elnet.coef <- elnet.fit$beta[, which(elnet.fit$lambda==elnet.cv$lambda.1se)]
  elnet.coef <- abs(elnet.coef)
  elnet.coef <- sort(elnet.coef, decreasing = TRUE)
  elnet.snps <- names(elnet.coef)
  elnet.nsnps <- length(elnet.coef[elnet.coef>=mean(elnet.coef)])
  elnet.snps <- elnet.snps[1:elnet.nsnps]
  #Export the SNPs extracted from Elastic net
  write(elnet.snps, "elnetsnps.txt")
  
  #Using training set#
  #Predictions
  predictions_train_elnet <- predict(elnet.fit, s = elnet.cv$lambda.1se, newx = x_train)

  #Precision
  cm_elnet = as.matrix(table(Actual = y_train, Predicted = predictions_train_elnet))
  precision_elnet_train <- precision(cm_elnet)
  print(paste(mean(precision_elnet_train), "precision_elnet_train"))
  
  #Recall
  recall_elnet_train <- recall(cm_elnet)
  print(paste(mean(recall_elnet_train), "recall_elnet_train"))
  
  #F1-score
  print(paste(F1_Score(predictions_train_elnet, y_train), "F1 train elnet"))

  #AUC
  auc_elnet_train <- auc(y_train, predictions_train_elnet)
  print(paste(auc_elnet_train, "auc_elnet_train"))
  
  #RMSE
  rmse_elnet_train <- RMSE(y_train, predictions_train_elnet)
  print(paste(rmse_elnet_train, "rmse_elnet_train"))
  
  #R-squared
  rsq_elnet_train <- RSQUARE(y_train, predictions_train_elnet)
  print(paste(rsq_elnet_train, "rsq_elnet_train"))
  
  #Accuracy
  print(paste(accuracy(predictions_train_elnet, y_train), "Accuracy train elnet"))
  
  #Using testing set#
  #Predictions
  predictions_test_elnet <- predict(elnet.fit, s = elnet.cv$lambda.1se, newx = x_test)

  #Precision
  cm_elnet = as.matrix(table(Actual = y_test, Predicted = predictions_test_elnet))
  precision_elnet_test <- precision(cm_elnet)
  print(paste(mean(precision_elnet_test), "precision_elnet_test"))
  
  #Recall
  recall_elnet_test <- recall(cm_elnet)
  print(paste(mean(recall_elnet_test), "recall_elnet_test"))
  
  #F1-score
  print(paste(F1_Score(predictions_test_elnet, y_test), "F1 test elnet"))

  #AUC
  auc_elnet_test <- auc(y_test, predictions_test_elnet)
  print(paste(auc_elnet_test, "auc_elnet_test"))
  
  #RMSE
  rmse_elnet_test <- RMSE(y_test, predictions_test_elnet)
  print(paste(rmse_elnet_test, "rmse_elnet_test"))
  
  #R-squared
  rsq_elnet_test <- RSQUARE(y_test, predictions_test_elnet)
  print(paste(rsq_elnet_test, "rsq_elnet_test"))
  
  #Accuracy
  print(paste(accuracy(predictions_test_elnet, y_test), "Accuracy test elnet"))
  ########

  # Create SNP pool (Union of 3 methods)
  ##############
  snp_pool1 <- unique(c(ridge.snps,lasso.snps,elnet.snps))
  #Export SNPs extracted from SNP pool
  write(snp_pool1, "snp_pool1.txt")
  ########

  #Create data compatible with group LASSO and SGL
  ##############
  Z.train <- x_train[,snp_pool1]
  y.train <- y_train
  Z.test <- x_test[,snp_pool1]
  y.test <- (y_test)

  #Hierarchical Clustering for group formation
  snps_df_sc <- as.data.frame(scale(t(Z.train)))
  dist_mat <- dist(snps_df_sc, method = 'euclidean')
  hclust_avg <- hclust(dist_mat, method = 'average')
  cut_avg <- cutree(hclust_avg, k = 200)

  Z.train.sd <- colSds(Z.train)
  Z.train    <- Z.train[, which(Z.train.sd > 1e-7)]
  Z.test     <- Z.test[, which(Z.train.sd > 1e-7)]
  grp <- cut_avg[which(Z.train.sd > 1e-7)]
  Z.train.sd <- Z.train.sd[which(Z.train.sd > 1e-7)]
  ########

  # Group Lasso
  ##############
  #Data preparation for Group LASSO
  Z.train <- scale(Z.train, center = TRUE, scale = Z.train.sd)
  grpYTrain <- y.train
  grpYTest <- y.test
  grpXTrain <- Z.train
  grpXTest <- Z.test
  
  #Training the model
  grplasso.cv <- cv.gglasso(grpXTrain, grpYTrain, nfolds = 5,
                            group = grp)
  grplasso.fit <- gglasso(grpXTrain, grpYTrain,
                          lambda = grplasso.cv$lambda.1se, group = grp)

  # Sort and Extract Coefficients and corresponding SNPs
  grplasso.coef <- coef(grplasso.fit)
  names(grplasso.coef) <- colnames(grpXTrain)
  grplasso.coef <- abs(grplasso.coef)
  grplasso.coef <- sort(grplasso.coef, decreasing = TRUE)
  grplasso.snps <- names(grplasso.coef)
  grplasso.snps <- grplasso.snps[grplasso.snps != ""]
  grplasso.nsnps <- length(grplasso.coef[grplasso.coef >= mean(grplasso.coef)])
  if(grplasso.nsnps > ncol(x_train))
  {
    grplasso.nsnps <- ncol(x_train)
  }
  grplasso.snps <- grplasso.snps[1:grplasso.nsnps]
  #Export SNPs extracted from Group LASSO
  write(grplasso.snps, "grplassosnps.txt")
  
  #Using training set#
  #Predictions
  predictions_train_grplasso <- predict(grplasso.fit, s = grplasso.cv$lambda.1se,
                                        newx = grpXTrain)

  #Precision
  cm_grplasso = as.matrix(table(Actual = grpYTrain, Predicted = predictions_train_grplasso))
  precision_grplasso_train <- precision(cm_grplasso)
  print(paste(mean(precision_grplasso_train), "precision_grplasso_train"))
  
  #Recall
  recall_grplasso_train <- recall(cm_grplasso)
  print(paste(mean(recall_grplasso_train), "recall_grplasso_train"))
  
  #F1-score
  print(paste(F1_Score(predictions_train_grplasso, grpYTrain), "F1 train grplasso"))

  #AUC
  auc_grplasso_train <- auc(grpYTrain, predictions_train_grplasso)
  print(paste(auc_grplasso_train, "auc_grplasso_train"))
  
  #RMSE
  rmse_grplasso_train <- RMSE(grpYTrain, predictions_train_grplasso)
  print(paste(rmse_grplasso_train, "rmse_grplasso_train"))
  
  #R-squared
  rsq_grplasso_train <- RSQUARE(grpYTrain, predictions_train_grplasso)
  print(paste(rsq_grplasso_train, "rsq_grplasso_train"))
  
  #Accuracy
  print(paste(accuracy(predictions_train_grplasso, grpYTrain), "Accuracy train grplasso"))
  
  #Using testing set#
  #Predictions
  predictions_test_grplasso <- predict(grplasso.fit, s = grplasso.cv$lambda.1se,
                                       newx = grpXTest)

  #Precision 
  cm_grplasso = as.matrix(table(Actual = grpYTest, Predicted = predictions_test_grplasso))
  precision_grplasso_test <- precision(cm_grplasso)
  print(paste(mean(precision_grplasso_test), "precision_grplasso_test"))
  
  #Recall
  recall_grplasso_test <- recall(cm_grplasso)
  print(paste(mean(recall_grplasso_test), "recall_grplasso_test"))
  
  #F1-score
  print(paste(F1_Score(predictions_test_grplasso, grpYTest), "F1 test grplasso"))

  #AUC
  auc_grplasso_test <- auc(grpYTest, predictions_test_grplasso)
  print(paste(auc_grplasso_test, "auc_grplasso_test"))
  
  #RMSE
  rmse_grplasso_test <- RMSE(grpYTest, predictions_test_grplasso)
  print(paste(rmse_grplasso_test, "rmse_grplasso_test"))
  
  #R-squared
  rsq_grplasso_test <- RSQUARE(grpYTest, predictions_test_grplasso)
  print(paste(rsq_grplasso_test, "rsq_grplasso_test"))
  
  #Accuracy
  print(paste(accuracy(predictions_test_grplasso, grpYTest), "Accuracy test grplasso"))
  
  #######
  
  # SGL
  ##############
  #Data Preparations for SGL
  Z.train <- scale(Z.train, center = TRUE, scale = Z.train.sd)
  data.SGL <- list(x = Z.train, y = y.train)
  set.seed(100)
  
  #Training the model
  sgl.cv <- cvSGL(data = data.SGL, index = grp, nfold= 5)
  sgl.fit <- SGL(data = data.SGL, index = grp)
  # Sort and Extract Coefficients and corresponding SNPs
  beta <- sgl.fit$beta[,20]
  sgl.coef <- data.frame(beta)
  row.names(sgl.coef) <- colnames(Z.test)
  sgl.coef <- abs(sgl.coef)
  sgl.coef <- sort(sgl.coef, decreasing = TRUE)
  sgl.snps <- row.names(sgl.coef)[apply(sgl.coef, 1, function(u) any(u>0.5))]
  #Export the extracted SNPs from SGL
  write(sgl.snps, "sglsnps.txt")
  
  #Using training set#
  #Predictions
  sgl.probabilities.train <- predictSGL(sgl.fit, lam = which.min(sgl.cv$lambdas), newX = Z.train)

  #Precision
  cm_sgl = as.matrix(table(Actual = y.train, Predicted = sgl.probabilities.train))
  precision_sgl_train <- precision(cm_sgl)
  print(paste(mean(precision_sgl_train), "precision_sgl_train"))
  
  #Recall
  recall_sgl_train <- recall(cm_sgl)
  print(paste(mean(recall_sgl_train), "recall_sgl_train"))
  
  #F1-score
  print(paste(F1_Score(y.train, sgl.probabilities.train), "F1 train sgl"))

  #AUC
  auc_sgl_train <- auc((y.train), sgl.probabilities.train)
  print(paste(auc_sgl_train, "auc_sgl_train"))
  
  #RMSE
  rmse_sgl_train <- RMSE(y.train, sgl.probabilities.train)
  print(paste(rmse_sgl_train, "rmse_sgl_train"))
  
  #R-squared
  rsq_sgl_train <- RSQUARE(y.train, sgl.probabilities.train)
  print(paste(rsq_sgl_train, "rsq_sgl_train"))
  
  #Accuracy
  print(paste(accuracy(y.train, sgl.probabilities.train), "Accuracy train sgl"))
  
  #Using testing set#
  #Predictions
  sgl.probabilities.test <- predictSGL(sgl.fit, lam = which.min(sgl.cv$lambdas), newX = Z.test)
  
  #Precision
  cm_sgl = as.matrix(table(Actual = y.test, Predicted = sgl.probabilities.test))
  precision_sgl_test <- precision(cm_sgl)
  print(paste(mean(precision_sgl_test), "precision_sgl_test"))
  
  #Recall
  recall_sgl_test <- recall(cm_sgl)
  print(paste(mean(recall_sgl_test), "recall_sgl_test"))
  
  #F1-score
  print(paste(F1_Score(y.test, sgl.probabilities.test), "F1 test sgl"))
  
  #AUC
  auc_sgl_test <- auc((y.test), sgl.probabilities.test)
  print(paste(auc_sgl_test, "auc_sgl_test"))
  
  #RMSE
  rmse_sgl_test <- RMSE(y.test, sgl.probabilities.test)
  print(paste(rmse_sgl_test, "rmse_sgl_test"))
  
  #R-squareed
  rsq_sgl_test <- RSQUARE(y.test, sgl.probabilities.test)
  print(paste(rsq_sgl_test, "rsq_sgl_test"))
  
  #Accuracy
  print(paste(accuracy(y.test, sgl.probabilities.test), "Accuracy test sgl"))
  ########
  
  #Union of SNPs extracted from Group LASSO and SGL
  ###########
  snp_set <- c(grplasso.snps,sgl.snps)
}
