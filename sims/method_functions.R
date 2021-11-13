## @knitr methods

error_return <- list(beta = NA,
                     vnames = NA,
                     nonzero_coef = NA,
                     active = c(""),
                     not_active = NA,
                     yvalid_hat = NA,
                     msevalid = NA,
                     causal = NA,
                     not_causal = NA,
                     yvalid = NA)

# All pairs Method
# lasso back tracking split -----------------------------------------------

lassoBTsplit <- new_method("lassoBT", "LassoBT",
                           method = function(model, draw) {
                             
                             tryCatch({
                               
                               # fitBT <- cvLassoBT(x = draw[["xtrain_lasso"]],
                               #                    y = draw[["ytrain"]], iter_max=10, nperms=1, nfolds = 10, nlambda = 100)
                               fitBT <- LassoBT(x = draw[["xtrain_lasso"]],
                                                y = draw[["ytrain"]], iter_max=10, nlambda = 100)
                               
                               ytest_hat <- predict(fitBT, newx = draw[["xtest_lasso"]])
                               msetest <- colMeans((draw[["ytest"]] - ytest_hat)^2)
                               lambda.min.index <- which(msetest==min(msetest), arr.ind = TRUE)
                               lambda.min <- fitBT$lambda[lambda.min.index[1,1]]
                               iter.min <- lambda.min.index[1,2]
                               
                               coefBT <- as.matrix(predict(fitBT, type = "coef",
                                                           s = lambda.min, iter = iter.min))
                               nzcoef <- coefBT[sail:::nonzero(coefBT),,drop=F]
                               
                               ints <- grep(":",rownames(nzcoef)[-1], value=T)
                               if (length(ints) != 0) {
                                 inters <- sapply(stringr::str_split(ints, ":"), function(i) paste(colnames(draw[["xtrain_lasso"]])[as.numeric(as.character(i))], collapse = ":"))
                                 mains <- colnames(draw[["xtrain_lasso"]])[as.numeric(as.character(setdiff(rownames(nzcoef)[-1], ints)))]
                               } else {
                                 inters <- character(0)
                                 mains <- colnames(draw[["xtrain_lasso"]])[as.numeric(as.character(rownames(nzcoef)[-1]))]
                               }
                               
                               active <- if (length(c(mains, inters)) == 0) " " else c(mains, inters)
                               yvalid_hat <- tryCatch({
                                 as.matrix(predict(fitBT, newx = draw[["xvalid_lasso"]], s = lambda.min, iter = iter.min, type = "response"))
                               },
                               error = function(err) {
                                 return(matrix(0, nrow = nrow(draw[["xvalid_lasso"]]), ncol = 1))
                               } # return NULL on error
                               )
                               
                               msevalid <- mean((draw[["yvalid"]] - drop(yvalid_hat))^2)
                               
                               return(list(beta = coefBT,
                                           # cvfit = fitBT,
                                           vnames = draw[["vnames_lasso"]],
                                           nonzero_coef = nzcoef,
                                           active = active,
                                           not_active = setdiff(draw[["vnames"]], active),
                                           yvalid_hat = yvalid_hat,
                                           msevalid = msevalid,
                                           causal = draw[["causal"]],
                                           not_causal = draw[["not_causal"]],
                                           yvalid = draw[["yvalid"]]))
                             },
                             error = function(err) {
                               return(error_return)
                             }
                             )
                           })



# GxE Method



# sail split --------------------------------------------------------------



sailsplit <- new_method("sail_split", "Sail_split",
                        method = function(model, draw) {
                          tryCatch({
                            fit <- sail(x = draw[["xtrain"]], y = draw[["ytrain"]], e = draw[["etrain"]],
                                        basis = function(i) splines::bs(i, degree = 5))
                            
                            ytest_hat <- predict(fit, newx = draw[["xtest"]], newe = draw[["etest"]])
                            msetest <- colMeans((draw[["ytest"]] - ytest_hat)^2)
                            lambda.min.index <- as.numeric(which.min(msetest))
                            lambda.min <- fit$lambda[which.min(msetest)]
                            
                            yvalid_hat <- predict(fit, newx = draw[["xvalid"]], newe = draw[["evalid"]], s = lambda.min)
                            msevalid <- mean((draw[["yvalid"]] - drop(yvalid_hat))^2)
                            
                            nzcoef <- predict(fit, s = lambda.min, type = "nonzero")
                            
                            return(list(beta = coef(fit, s = lambda.min)[-1,,drop=F],
                                        # fit = fit,
                                        vnames = draw[["vnames"]],
                                        nonzero_coef = nzcoef,
                                        active = fit$active[[lambda.min.index]],
                                        not_active = setdiff(draw[["vnames"]], fit$active[[lambda.min.index]]),
                                        yvalid_hat = yvalid_hat,
                                        msevalid = msevalid,
                                        causal = draw[["causal"]],
                                        not_causal = draw[["not_causal"]],
                                        yvalid = draw[["yvalid"]]))
                          },
                          error = function(err) {
                            return(error_return)
                          }
                          )
                        })


# sail split linear -------------------------------------------------------


sailsplitlinear <- new_method("linearsail", "Linear Sail",
                              method = function(model, draw) {
                                tryCatch({
                                  fit <- sail(x = draw[["xtrain"]], y = draw[["ytrain"]], e = draw[["etrain"]],
                                              basis = function(i) i)
                                  
                                  ytest_hat <- predict(fit, newx = draw[["xtest"]], newe = draw[["etest"]])
                                  msetest <- colMeans((draw[["ytest"]] - ytest_hat)^2)
                                  lambda.min.index <- as.numeric(which.min(msetest))
                                  lambda.min <- fit$lambda[which.min(msetest)]
                                  
                                  yvalid_hat <- predict(fit, newx = draw[["xvalid"]], newe = draw[["evalid"]], s = lambda.min)
                                  msevalid <- mean((draw[["yvalid"]] - drop(yvalid_hat))^2)
                                  
                                  nzcoef <- predict(fit, s = lambda.min, type = "nonzero")
                                  
                                  return(list(beta = coef(fit, s = lambda.min)[-1,,drop=F],
                                              # fit = fit,
                                              vnames = draw[["vnames"]],
                                              nonzero_coef = nzcoef,
                                              active = fit$active[[lambda.min.index]],
                                              not_active = setdiff(draw[["vnames"]], fit$active[[lambda.min.index]]),
                                              yvalid_hat = yvalid_hat,
                                              msevalid = msevalid,
                                              causal = draw[["causal"]],
                                              not_causal = draw[["not_causal"]],
                                              yvalid = draw[["yvalid"]]))
                                },
                                error = function(err) {
                                  return(error_return)
                                }
                                )
                              })

# Gesso Split -------------------------------------------------------


gessosplit <- new_method("Gesso Split", "Gesso Split",
                              method = function(model, draw) {
                                tryCatch({
                                  gesso_fit <- gesso.fit(G = draw[["xtrain"]], E = draw[["etrain"]], Y = draw[["ytrain"]],
                                                   normalize=TRUE,family = "gaussian", min_working_set_size = 30)
                                  # "gaussian" for continuous outcome and "binomial" for binary
                                  coefficients = gesso.coef(fit=gesso_fit$fit, lambda=gesso_fit$lambda_min)
                                  # predicted_value <-predict
                                  beta_0 = coefficients$beta_0
                                  beta_e = coefficients$beta_e
                                  beta_g = coefficients$beta_g
                                  beta_gxe = coefficients$beta_gxe
                                  new_G = draw[["xtest"]]
                                  new_E = draw[["etest"]]
                                  ytest_hat = gesso.predict(beta_0, beta_e, beta_g, beta_gxe, new_G, new_E)
                                  # ytest_hat <- predict(fit, newx = draw[["xtest"]], newe = draw[["etest"]])
                                  msetest <- colMeans((draw[["ytest"]] - ytest_hat)^2)
                                  lambda.min.index <- as.numeric(which.min(msetest))
                                  # lambda.min <- fit$lambda[which.min(msetest)]

                                  new_G_valid = draw[["xvalid"]]
                                  new_E_valid = draw[["evalid"]]
                                  yvalid_hat = gesso.predict(beta_0, beta_e, beta_g, beta_gxe, new_G_valid, new_E_valid)

                                  msevalid <- mean((draw[["yvalid"]] - drop(yvalid_hat))^2)

                                  # nzcoef <- predict(gesso_fit, s = gesso_fit$lambda_min, type = "nonzero")
                                  nzcoef <- beta_0


                                  return(list(beta = coefficients,
                                              # fit = fit,
                                              vnames = draw[["vnames"]],
                                              nonzero_coef = nzcoef, # beta coeff
                                              active = fit$active[[lambda.min.index]], # index, name of index 
                                              not_active = setdiff(draw[["vnames"]], fit$active[[lambda.min.index]]),
                                              yvalid_hat = yvalid_hat,
                                              msevalid = msevalid,
                                              causal = draw[["causal"]],
                                              not_causal = draw[["not_causal"]],
                                              yvalid = draw[["yvalid"]]))
                                },
                                error = function(err) {
                                  return(error_return)
                                }
                                )
                              })


# 
# 
# # glinternet split --------------------------------------------------------
# 
# GLinternetsplit <- new_method("GLinternet", "GLinternet",
#                               method = function(model, draw) {
#                                 
#                                 tryCatch({
#                                   fitGL <- glinternet(X = draw[["xtrain_lasso"]], Y = draw[["ytrain"]],
#                                                       numLevels = rep(1, ncol(draw[["xtrain_lasso"]])),
#                                                       nLambda = 100, interactionCandidates = c(1),
#                                                       verbose = F)
#                                   
#                                   ytest_hat <- predict(fitGL, X = draw[["xtest_lasso"]])
#                                   msetest <- colMeans((draw[["ytest"]] - ytest_hat)^2)
#                                   lambda.min.index <- as.numeric(which.min(msetest))
#                                   lambda.min <- fitGL$lambda[which.min(msetest)]
#                                   
#                                   yvalid_hat <- predict(fitGL, X = draw[["xvalid_lasso"]], lambda = lambda.min)
#                                   msevalid <- mean((draw[["yvalid"]] - drop(yvalid_hat))^2)
#                                   
#                                   # nzcoef <- coef(fit, s = lambda.min)[nonzeroCoef(coef(fit, s = lambda.min)),,drop=F]
#                                   
#                                   tc <- coef(fitGL, lambdaIndex = lambda.min.index)
#                                   mains <- colnames(draw[["xtrain_lasso"]])[tc[[1]]$mainEffects$cont]
#                                   inters <- paste0(colnames(draw[["xtrain_lasso"]])[tc[[1]]$interactions$contcont[,2]],":E")
#                                   
#                                   return(list(beta = NULL,
#                                               vnames = draw[["vnames_lasso"]],
#                                               nonzero_coef = NULL,
#                                               active = c(mains, inters),
#                                               not_active = setdiff(draw[["vnames"]], c(mains, inters)),
#                                               yvalid_hat = yvalid_hat,
#                                               msevalid = msevalid,
#                                               causal = draw[["causal"]],
#                                               not_causal = draw[["not_causal"]],
#                                               yvalid = draw[["yvalid"]]))
#                                 },
#                                 error = function(err) {
#                                   return(error_return)
#                                 }
#                                 )
#                               })
# 
# 
