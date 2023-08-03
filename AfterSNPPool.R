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