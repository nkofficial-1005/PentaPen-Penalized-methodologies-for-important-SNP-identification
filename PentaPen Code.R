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
library(foreach)
library(doParallel)
library(grpreg)
##############
# Import data
Geno <- read.pedfile("genotype.ped") 
char.pheno <- read.table("phenotypes.pheno", header = TRUE, stringsAsFactors = FALSE, sep = " ")

##############
# Data Pre-processing
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
df_final <- df[, 7:204760] #select the data set consisting of SNPs only
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
#Perform 5 fold cross validation and parallel computation
for(i in 1){
  
  print(paste("Fold:",i))
  #Segement your data by fold using the which() function 
  train_index <- sample(nrow(df_final), 0.7*nrow(df_final))
  train_data <- df_final[train_index,]
  val_data <- df_final[-train_index,]
  x_train <- train_data[,-ncol(train_data)]
  y_train  <- train_data[,ncol(train_data)]
  x_test <- val_data[,-ncol(val_data)]
  y_test <- val_data[,ncol(val_data)]
  nfolds <- 5
  
  #Parallel computing for Ridge, LASSO, Elastic net
  cl <- makeCluster(detectCores() - 1) # create a cluster with all but one core
  registerDoParallel(cl) # register the cluster with foreach
  ridge_lasso_elnet <- foreach(i=1:3, .packages="glmnet") %dopar% {
    if(i == 1) { # ridge 
      cv_ridge <- cv.glmnet(x=x_train,
                            y=y_train,
                            alpha=0, nfolds=nfolds)
      fit_ridge <- glmnet(x=x_train,
                          y=y_train,
                          alpha=0, lambda = cv_ridge$lambda.1se)
      return(list(coef=coef(cv_ridge, s="lambda.min"), # coefficients at lambda with minimum cross-validation error
                  cv_error=cv_ridge$cvm[cv_ridge$lambda == cv_ridge$lambda.1se], # cross-validation error at minimum lambda
                  lambda=cv_ridge$lambda,
                  lambdase=cv_ridge$lambda.1se,
                  predictions_train_ridge = predict(fit_ridge, s = cv_ridge$lambda.1se, newx = x_train),
                  predictions_test_ridge = predict(fit_ridge, s = cv_ridge$lambda.1se, newx = x_test),
                  selected_nsnps=which(coef(cv_ridge, s="lambda.min")!=0)))
    } else if(i == 2) { # lasso 
      cv_lasso <- cv.glmnet(x=x_train,
                            y=y_train,
                            alpha=1, nfolds=nfolds)
      fit_lasso <- glmnet(x_train, y_train, alpha=1, 
                          lambda = cv_lasso$lambda.1se)
      return(list(coef=coef(cv_lasso, s="lambda.min"), # coefficients at lambda with minimum cross-validation error
                  cv_error=cv_lasso$cvm[cv_lasso$lambda == cv_lasso$lambda.1se], # cross-validation error at minimum lambda
                  lambda=cv_lasso$lambda,
                  lambdase=cv_lasso$lambda.1se,
                  predictions_train_lasso = predict(fit_lasso, s = cv_lasso$lambda.1se, newx = x_train),
                  predictions_test_lasso = predict(fit_lasso, s = cv_lasso$lambda.1se, newx = x_test),
                  selected_nsnps=which(coef(cv_lasso, s="lambda.min")!=0)))
    } else { # elastic net 
      cv_enet <- cv.glmnet(x=train_data[,-ncol(train_data)],
                           y=train_data[,ncol(train_data)],
                           alpha=0.5, nfolds=nfolds)
      fit_enet <- glmnet(x_train, y_train, alpha=0.5, 
                          lambda = cv_enet$lambda.1se)
      return(list(coef=coef(cv_enet, s="lambda.min"), # coefficients at lambda with minimum cross-validation error
                  cv_error=cv_enet$cvm[cv_enet$lambda == cv_enet$lambda.1se], # cross-validation error at minimum lambda
                  lambda=cv_enet$lambda,
                  lambdase=cv_enet$lambda.1se,
                  predictions_train_enet = predict(fit_enet, s = cv_enet$lambda.1se, newx = x_train),
                  predictions_test_enet = predict(fit_enet, s = cv_enet$lambda.1se, newx = x_test),
                  selected_nsnps=which(coef(cv_enet, s="lambda.min")!=0)))
    }
  }
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
  grp_sgl <- foreach(i=1:2, .packages=c("gglasso","SGL")) %dopar% {
    if(i == 1) { # group lasso
      cv_glasso <- cv.gglasso(x=Z.train,
                              y=y.train,
                              group=grp, nfolds=nfolds)
      fit_glasso <- gglasso(x=Z.train,
                              y=y.train,
                              lambda = cv_glasso$lambda.1se, group = grp)
      return(list(coef=coef(cv_glasso, s="lambda.min"), # coefficients at lambda with minimum cross-validation error
                  cv_error=cv_glasso$cvm[cv_glasso$lambda == cv_glasso$lambda.1se], # cross-validation error at minimum lambda
                  predictions_train_grplasso=predict(fit_glasso, s = cv_glasso$lambda.1se,
                                                     newx = Z.train),
                  predictions_test_grplasso=predict(fit_glasso, s = cv_glasso$lambda.1se,
                                                       newx = Z.test),
                  selected_snps=which(coef(cv_glasso, s="lambda.min")!=0))) # selected SNPs at minimum lambda
    } else { # sparse group lasso
      cv_sgl <- cvSGL(data = data.SGL, index = grp, nfold=nfolds)
      sgl_fit <- SGL(data = data.SGL, index = grp)
      return(list(coef=cv_sgl$fit$beta[,which.min(cv_sgl$lambdas)],#coef(cv_sgl, s=which.min(cv_sgl$lambdas)), # coefficients at lambda with minimum cross-validation error
                  lambda=which.min(cv_sgl$lambdas),
                  lam=cv_sgl$lambdas,
                  sgl.probabilities.train=predictSGL(sgl_fit, lam = which.min(cv_sgl$lambdas), newX = Z.train_sgl),
                  sgl.probabilities.test =predictSGL(sgl_fit, lam = which.min(cv_sgl$lambdas), newX = Z.test_sgl)))
    }
  }
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
# #Recall
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
