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



gessosplit <- new_method("gesso","Gesso Split",
                           method = function(model, draw) {
                             
                             tryCatch({
                               
                               fit_model <- gesso.fit(G = draw[["x_train"]], E = draw[["e_train"]], Y = draw[["y_train"]],
                                                      normalize=TRUE,family = "gaussian", min_working_set_size = 30)
                               
                               
                               msetest<-1000000
                               lambda1_index<- -1
                               lambda2_index<- -1
                               for(i in 1:(length(fit_model$lambda_1))){
                                 for(j in 1:(length(fit_model$lambda_2))){
                                   lambda_pair<-tibble(lambda_1 =fit_model$lambda_1[i],lambda_2 = fit_model$lambda_2[j])
                                   coeff<-gesso.coef(fit_model,lambda_pair)
                                   y_test_hat<-gesso.predict(coeff$beta_0,coeff$beta_e,coeff$beta_g,coeff$beta_gxe,draw[["x_test"]], draw[["e_test"]])
                                   temp_mse<-mean((y_test_hat-draw[["y_test"]])^2)
                                   
                                   if (temp_mse < msetest){
                                     msetest<-temp_mse
                                     lambda1_index<-i
                                     lambda2_index<-j
                                   }
                                 }
                               }
                               
                               min_lambdapair <- tibble(lambda_1 =fit_model$lambda_1[lambda1_index],lambda_2 = fit_model$lambda_2[lambda1_index])
                               coeff <- gesso.coef(fit_model,min_lambdapair)
                               y_valid_hat <- gesso.predict(coeff$beta_0,coeff$beta_e,coeff$beta_g,coeff$beta_gxe,draw[["x_valid"]], draw[["e_valid"]])
                               mse_valid <- mean((y_valid_hat-draw[["y_valid"]])^2)
                               
                               active_name_list<-c() # active variable name
                               beta_Matrix<-matrix(0,nrow=length(draw[["vnames"]])+1,ncol=1) # variable matrix with name and value
                               rowname_list<-append("(Intercept)",draw[["vnames"]])
                               rownames(beta_Matrix) <- rowname_list
                               colnames(beta_Matrix) <- "1"
                               beta_Matrix[1,]<-coeff$beta_0
                               
                               # loop for main effect
                               for (i in (1:(length(coeff$beta_g))) ){
                                 if ((coeff$beta_g)[i] !=0 ){
                                   vanme<-paste("X",i,sep="")
                                   active_name_list<-append(active_name_list,vanme)
                                   beta_Matrix[i+1,]<-(coeff$beta_g)[i]
                                 }
                               }
                               # betaE
                               if ((coeff$beta_e) !=0 ){
                                 active_name_list<-append(active_name_list,"E")
                                 beta_Matrix[length(coeff$beta_g)+2,]<-coeff$beta_e
                               }
                               # loop for interaction effect
                               for (i in (1:(length(coeff$beta_gxe))) ){
                                 if ((coeff$beta_gxe)[i] !=0 ){
                                   xname<-paste("X",i,sep="")
                                   vname<-paste(xname,":E",sep="")
                                   active_name_list<-append(active_name_list,vname)
                                   beta_Matrix[length(coeff$beta_g)+2+i,]<-coeff$beta_gxe[i]
                                 }
                               }
                               
                               return(list(beta = as.matrix(beta_Matrix[-1,1]),
                                           vnames = draw[["vnames"]],
                                           nonzero_coef = as.matrix(beta_Matrix[beta_Matrix[, 1] != 0, ]),
                                           active = active_name_list,
                                           not_active = setdiff(draw[["vnames"]], active_name_list),
                                           yvalid_hat = y_valid_hat,
                                           msevalid = mse_valid,
                                           causal = draw[["causal"]],
                                           not_causal = draw[["not_causal"]],
                                           yvalid = draw[["y_valid"]]))
                             },
                             error = function(err) {
                               return(error_return)
                             }
                             )
                           })


sailsplit <- new_method("sail", "Sail Split",
                         method = function(model, draw) {
                           
                           tryCatch({
                             
                             fit <- sail(x = draw[["x_train"]], y = draw[["y_train"]], e = draw[["e_train"]],
                                         basis = function(i) i)
                             
                             ytest_hat <- predict(fit, newx = draw[["x_test"]], newe = draw[["e_test"]])
                             msetest <- colMeans((draw[["y_test"]] - ytest_hat)^2)
                             lambda.min.index <- as.numeric(which.min(msetest))
                             lambda.min <- fit$lambda[which.min(msetest)]
                             
                             yvalid_hat <- predict(fit, newx = draw[["x_valid"]], newe = draw[["e_valid"]], s = lambda.min)
                             msevalid <- mean((draw[["y_valid"]] - drop(yvalid_hat))^2)
                             
                             nzcoef <- predict(fit, s = lambda.min, type = "nonzero")
                             
                             return(list(beta = coef(fit, s = lambda.min)[-1,,drop=F],
                                         vnames = draw[["vnames"]],
                                         nonzero_coef = nzcoef,
                                         active = fit$active[[lambda.min.index]],
                                         not_active = setdiff(draw[["vnames"]], fit$active[[lambda.min.index]]),
                                         yvalid_hat = yvalid_hat,
                                         msevalid = msevalid,
                                         causal = draw[["causal"]],
                                         not_causal = draw[["not_causal"]],
                                         yvalid = draw[["y_valid"]]))
                           },
                           error = function(err) {
                             return(error_return)
                           }
                           )
                         })


glinternetsplit_XE <- new_method("glinternetsplit_XE", "glinternetsplit XE",
                        method = function(model, draw) {
                          
                          tryCatch({
                            fitGL <- glinternet(X = draw[["x_train_lasso"]], Y = draw[["y_train"]],
                                                numLevels = rep(1, ncol(draw[["x_train_lasso"]])),
                                                nLambda = 100, interactionCandidates = c(1),
                                                verbose = F)
                            
                            y_test_hat <- predict(fitGL, X = draw[["x_test_lasso"]])
                            mse_test <- colMeans((draw[["y_test"]] - y_test_hat)^2)
                            lambda.min.index <- as.numeric(which.min(mse_test))
                            lambda.min <- fitGL$lambda[which.min(mse_test)]
                            
                            y_valid_hat <- predict(fitGL, X = draw[["x_valid_lasso"]], lambda = lambda.min)
                            mse_valid <- mean((draw[["y_valid"]] - drop(y_valid_hat))^2)
        
                            fit_coef <- coef(fitGL, lambdaIndex = lambda.min.index)
                            mains <- colnames(draw[["x_train_lasso"]])[fit_coef[[1]]$mainEffects$cont]
                            inters <- paste0(colnames(draw[["x_train_lasso"]])[fit_coef[[1]]$interactions$contcont[,2]],":E")
                            
                            beta_Matrix<-matrix(0,nrow=length(draw[["vnames_lasso"]]),ncol=1) # variable matrix with name and value
                            rownames(beta_Matrix) <- draw[["vnames_lasso"]]
  
                            for (i in (1:length(mains))){
                              beta_Matrix[mains[i],]<-fit_coef[[1]]$mainEffectsCoef$cont[[i]]
                            }
                            
                            for (i in (1:length(inters))){
                              beta_Matrix[inters[i],]<-fit_coef[[1]]$interactionsCoef$contcont[[i]]
                            }
                            
                            return(list(beta = beta_Matrix,
                                        vnames = draw[["vnames_lasso"]],
                                        nonzero_coef = as.matrix(beta_Matrix[beta_Matrix[, 1] != 0, ]),
                                        active = c(mains, inters),
                                        not_active = setdiff(draw[["vnames_lasso"]], c(mains, inters)),
                                        yvalid_hat = y_valid_hat,
                                        msevalid = mse_valid,
                                        causal = draw[["causal"]],
                                        not_causal = draw[["not_causal"]],
                                        yvalid = draw[["y_valid"]]))
                          },
                          error = function(err) {
                            return(error_return)
                          }
                          )
                        })


# 
# # All pairs Method
# # lasso back tracking split -----------------------------------------------
# 
# lassoBTsplit <- new_method("lassoBT", "LassoBT",
#                            method = function(model, draw) {
#                              
#                              tryCatch({
#                                
#                                # fitBT <- cvLassoBT(x = draw[["xtrain_lasso"]],
#                                #                    y = draw[["ytrain"]], iter_max=10, nperms=1, nfolds = 10, nlambda = 100)
#                                fitBT <- LassoBT(x = draw[["xtrain_lasso"]],
#                                                 y = draw[["ytrain"]], iter_max=10, nlambda = 100)
#                                
#                                ytest_hat <- predict(fitBT, newx = draw[["xtest_lasso"]])
#                                msetest <- colMeans((draw[["ytest"]] - ytest_hat)^2)
#                                lambda.min.index <- which(msetest==min(msetest), arr.ind = TRUE)
#                                lambda.min <- fitBT$lambda[lambda.min.index[1,1]]
#                                iter.min <- lambda.min.index[1,2]
#                                
#                                coefBT <- as.matrix(predict(fitBT, type = "coef",
#                                                            s = lambda.min, iter = iter.min))
#                                nzcoef <- coefBT[sail:::nonzero(coefBT),,drop=F]
#                                
#                                ints <- grep(":",rownames(nzcoef)[-1], value=T)
#                                if (length(ints) != 0) {
#                                  inters <- sapply(stringr::str_split(ints, ":"), function(i) paste(colnames(draw[["xtrain_lasso"]])[as.numeric(as.character(i))], collapse = ":"))
#                                  mains <- colnames(draw[["xtrain_lasso"]])[as.numeric(as.character(setdiff(rownames(nzcoef)[-1], ints)))]
#                                } else {
#                                  inters <- character(0)
#                                  mains <- colnames(draw[["xtrain_lasso"]])[as.numeric(as.character(rownames(nzcoef)[-1]))]
#                                }
#                                
#                                active <- if (length(c(mains, inters)) == 0) " " else c(mains, inters)
#                                yvalid_hat <- tryCatch({
#                                  as.matrix(predict(fitBT, newx = draw[["xvalid_lasso"]], s = lambda.min, iter = iter.min, type = "response"))
#                                },
#                                error = function(err) {
#                                  return(matrix(0, nrow = nrow(draw[["xvalid_lasso"]]), ncol = 1))
#                                } # return NULL on error
#                                )
#                                
#                                msevalid <- mean((draw[["yvalid"]] - drop(yvalid_hat))^2)
#                                
#                                return(list(beta = coefBT,
#                                            # cvfit = fitBT,
#                                            vnames = draw[["vnames_lasso"]],
#                                            nonzero_coef = nzcoef,
#                                            active = active,
#                                            not_active = setdiff(draw[["vnames"]], active),
#                                            yvalid_hat = yvalid_hat,
#                                            msevalid = msevalid,
#                                            causal = draw[["causal"]],
#                                            not_causal = draw[["not_causal"]],
#                                            yvalid = draw[["yvalid"]]))
#                              },
#                              error = function(err) {
#                                return(error_return)
#                              }
#                              )
#                            })
# 
# 
