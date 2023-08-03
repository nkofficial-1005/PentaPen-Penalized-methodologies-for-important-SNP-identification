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