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


# # Environmental Variable Method-----------
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
                            
                            beta_Matrix<-matrix(0,nrow=length(draw[["vnames_lasso"]])+1,ncol=1)
                            rowname_list<-append("(Intercept)",draw[["vnames"]]) # keep the vnames order the same as other methods
                            rownames(beta_Matrix) <- rowname_list
                            colnames(beta_Matrix) <- "1"
                            beta_Matrix[1,]<- 0
  
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

lassoBTsplit <- new_method("lassoBT", "LassoBT",
                           method = function(model, draw) {
                             # Remove try-catch since the method will throw warning messg even you give the argument value in a correct format/
                             tryCatch({ 
                               fitBT <- LassoBT(x = draw[["x_train"]],
                                                y = draw[["y_train"]], iter_max=10, nlambda = 100)
                               ytest_hat <- predict(fitBT, newx = draw[["x_test"]])
                               msetest <- colMeans((draw[["y_test"]] - ytest_hat)^2)
                               lambda.min.index <- which(msetest==min(msetest), arr.ind = TRUE)
                               lambda.min <- fitBT$lambda[lambda.min.index[1,1]]
                               iter.min <- lambda.min.index[1,2]

                               coefBT <- as.matrix(predict(fitBT, type = "coef",
                                                           s = lambda.min, iter = iter.min))

                               ints <- grep(":",rownames(coefBT)[-1], value=T)
                               if (length(ints) != 0) {
                                 mains <- colnames(draw[["x_train"]])[as.numeric(as.character(setdiff(rownames(coefBT)[-1], ints)))]
                                 for(i in (1:length(ints))){ # change the name of vnames in row
                                     each_vname<-stringr::str_split(ints, ":")[[i]]
                                     a<-each_vname[1]
                                     b<-each_vname[2]
                                     if(strtoi(a)>strtoi(b)){
                                       print(paste(b, a, sep=":"))
                                       ints[i]<-paste(b, a, sep=":")
                                     }
                                 }
                                 inters <- sapply(stringr::str_split(ints, ":"), function(i) paste(colnames(draw[["x_train"]])[as.numeric(as.character(i))], collapse = ":"))
                                 
                               } else {
                                 inters <- character(0)
                                 mains <- colnames(draw[["x_train"]])[as.numeric(as.character(rownames(nzcoef)[-1]))]
                               }
                               rowname_list<-append("(Intercept)",mains)
                               rownames(coefBT)<-append(rowname_list,inters)
                               nzcoef <- coefBT[sail:::nonzero(coefBT),,drop=F]
                               active<-if (length(rownames(nzcoef)[-1]) == 0) " " else rownames(nzcoef)[-1]
                               yvalid_hat <- tryCatch({
                                 as.matrix(predict(fitBT, newx = draw[["x_valid"]], s = lambda.min, iter = iter.min, type = "response"))
                               },
                               error = function(err) {
                                 return(matrix(0, nrow = nrow(draw[["x_valid"]]), ncol = 1))
                               } # return NULL on error
                               )

                               msevalid <- mean((draw[["y_valid"]] - drop(yvalid_hat))^2)

                               return(list(beta = coefBT,
                                           vnames = draw[["vnames"]],
                                           nonzero_coef = nzcoef,
                                           active = active,
                                           not_active = setdiff(draw[["vnames"]], active),
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

glinternetsplit_XX <- new_method("glinternetsplit_XX", "glinternetsplit XX",
                                 method = function(model, draw) {
                                   
                                   tryCatch({
                                     fitGL <- glinternet(X = draw[["x_train"]], Y = draw[["y_train"]],
                                                         numLevels = rep(1, ncol(draw[["x_train"]])),
                                                         nLambda = 100,
                                                         verbose = F)
                                     
                                     y_test_hat <- predict(fitGL, X = draw[["x_test"]])
                                     mse_test <- colMeans((draw[["y_test"]] - y_test_hat)^2)
                                     lambda.min.index <- as.numeric(which.min(mse_test))
                                     lambda.min <- fitGL$lambda[which.min(mse_test)]
                                     
                                     # y_valid_hat <- predict(fitGL, X = draw[["x_valid"]], lambda = lambda.min)
                                     # mse_valid <- mean((draw[["y_valid"]] - drop(y_valid_hat))^2)
                                     # 
                                     fit_coef <- coef(fitGL, lambdaIndex = lambda.min.index)
                                     
                                     # beta_Matrix<-matrix(0,nrow=length(draw[["vnames"]])+1,ncol=1)
                                     # rowname_list<-append("(Intercept)",draw[["vnames"]]) # keep the vnames order the same as other methods
                                     # rownames(beta_Matrix) <- rowname_list
                                     # colnames(beta_Matrix) <- "1"
                                     # beta_Matrix[1,]<- 0
                                     # 
                                     # if (length(tc[[1]]$mainEffects$cont) > 0){
                                     #   mains <- colnames(draw[["x_train"]])[fit_coef[[1]]$mainEffects$cont]
                                     #   
                                     #   for (i in (1:length(mains))){
                                     #     beta_Matrix[mains[i],]<-fit_coef[[1]]$mainEffectsCoef$cont[[i]]
                                     #   }
                                     # }
                                     # if (length(tc[[1]]$interactions$contcont[,1]) > 0){
                                     #   inters <- paste(colnames(draw[["x_train"]])[fit_coef[[1]]$interactions$contcont[,1]],
                                     #                   colnames(draw[["x_train"]])[fit_coef[[1]]$interactions$contcont[,2]], sep=":")
                                     #   for (i in (1:length(inters))){
                                     #     beta_Matrix[inters[i],]<-fit_coef[[1]]$interactionsCoef$contcont[[i]]
                                     #   }
                                     #   
                                     # }
                                     
                                     y_valid_hat <- predict(fitGL, X = draw[["x_valid"]], lambda = lambda.min)
                                     mse_valid <- mean((draw[["y_valid"]] - y_valid_hat)^2)
                                     
                                     
                                     mains <- colnames(draw[["x_train"]])[fit_coef[[1]]$mainEffects$cont]
                                     inters <- paste(colnames(draw[["x_train"]])[fit_coef[[1]]$interactions$contcont[,1]],
                                                     colnames(draw[["x_train"]])[fit_coef[[1]]$interactions$contcont[,2]], sep=":")
                                     
                                     beta_Matrix<-matrix(0,nrow=length(draw[["vnames"]])+1,ncol=1)
                                     rowname_list<-append("(Intercept)",draw[["vnames"]]) # keep the vnames order the same as other methods
                                     rownames(beta_Matrix) <- rowname_list
                                     colnames(beta_Matrix) <- "1"
                                     beta_Matrix[1,]<- 0
                                     
                                     for (i in (1:length(mains))){
                                       tryCatch({
                                         beta_Matrix[mains[i],]<-fit_coef[[1]]$mainEffectsCoef$cont[[i]]
                                       },
                                       error = function(err) {
                                         continue
                                       })
                                     }
                                     
                                     for (i in (1:length(inters))){
                                       tryCatch({
                                         beta_Matrix[inters[i],]<-fit_coef[[1]]$interactionsCoef$contcont[[i]]
                                       },
                                       error = function(err) {
                                         continue
                                       })
                                     }
                                     
                                     return(list(beta = beta_Matrix,
                                                 vnames = draw[["vnames"]],
                                                 nonzero_coef = as.matrix(beta_Matrix[beta_Matrix[, 1] != 0, ]),
                                                 active = c(mains, inters),
                                                 not_active = setdiff(draw[["vnames"]], c(mains, inters)),
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

hiernet_split <- new_method("hinernet", "hinernet",
                                 method = function(model, draw) {
                                   
                                   tryCatch({
                                     
                                     fit_hierNet_path<-hierNet.path(draw[["x_train"]], draw[["y_train"]],strong=TRUE)
                                     
                                     ytest_hat <- predict(fit_hierNet_path, newx = draw[["x_test"]])
                                     col_msetest<-c()
                                     
                                     for (i in 1:length(ytest_hat[1,])){
                                       msetest <- mean((draw[["y_test"]] - ytest_hat[,i])^2)
                                       col_msetest<-append(col_msetest,msetest)
                                     }
                                     lambda_min_index<-which.min(col_msetest)
                                     lambda_min<-fit_hierNet_path$lamlist[which.min(col_msetest)]
                                     
                                     # bp:p-vector of estimated "positive part" main effect (p=# features)
                                     # bn:p-vector of estimated "negative part" main effect (p=# features)
                                     # bp-bn:overall main effect estimated
                                     beta_main<-fit_hierNet_path$bp[,lambda_min_index]-fit_hierNet_path$bn[,lambda_min_index] #overall main effect estimated
                                     as.vector(beta_main)
                                     
                                     beta_Matrix<-matrix(0,nrow=length(draw[["vnames"]])+1,ncol=1)
                                     rowname_list<-append("(Intercept)",draw[["vnames"]])
                                     rownames(beta_Matrix) <- rowname_list
                                     colnames(beta_Matrix) <- "1"
                                     beta_Matrix[1,]<-0
                                     
                                     num_p<-dim(draw[["x_train"]])[2]
                                     for (i in (1:num_p)){
                                       beta_Matrix[i+1]<-beta_main[i]
                                     }
                                     
                                     inter<-c()
                                     
                                     for (i in (1:(num_p-1))){
                                       for (j in ((i+1):num_p)){
                                         inter<-append(inter, fit_hierNet_path$th[,,lambda_min_index][i,j])
                                         
                                       }
                                     }
                                     
                                     beta_Matrix[(num_p+2):(length(beta_Matrix))]<-inter
                                     
                                     nzcoef<-as.matrix(beta_Matrix[beta_Matrix[, 1] != 0, ])
                                     
                                     active<-rownames(nzcoef)
                                     
                                     y_valid_hat <- predict(fit_hierNet_path, newx = draw[["x_valid"]])[,lambda_min_index]
                                     mse_valid <- mean((draw[["y_valid"]] - y_valid_hat)^2)
                       
                                     return(list(beta = beta_Matrix,
                                                 vnames = draw[["vnames"]],
                                                 nonzero_coef = nzcoef,
                                                 active = active,
                                                 not_active = setdiff(draw[["vnames"]], active),
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

sprintr_split <- new_method("sprintr", "sprintr",
                            method = function(model, draw) {
                                
                              tryCatch({

                                fit_sprintr<-sprinter(draw[["x_train"]],draw[["y_train"]])
                                y_test_hat<-predict(fit_sprintr,draw[["x_test"]])
                               
                                msetest<- 10000000
                                lambda1_min<- -1
                                lambda3_min<- -1
                                for (i in (2:length(fit_sprintr$lambda1))){
                                  for (j in (2:(dim(fit_sprintr$lambda3)[1]))){
                                    temp<-mean((draw[["y_test"]] - y_test_hat[[i]][,j])^2)
                                    if(temp < msetest){
                                      msetest<-temp
                                      lambda1_min<-i
                                      lambda3_min<-j
                                    }
                                  }
                                } # loop to find the best parameter value for lambda1 and lambda3

                                beta_Matrix<-matrix(0,nrow=length(draw[["vnames"]])+1,ncol=1)
                                rowname_list<-append("(Intercept)",draw[["vnames"]])
                                rownames(beta_Matrix) <- rowname_list
                                colnames(beta_Matrix) <- "1"
                                beta_Matrix[1,]<-fit_sprintr$step1$a0[lambda1_min]

                                num_p<-dim(draw[["x_train"]])[2]
                                # main beta
                                for(i in (1:num_p)){
                                  beta_Matrix[(i+1),]<-fit_sprintr$step1$beta[,lambda1_min][i]
                                }
                                # 
                                # inter beta
                                inter_names<-c()
                                fit_sprintr$step2[[lambda1_min]]
                                num_inter<-nrow(fit_sprintr$step2[[lambda1_min]])
                                inter_value<-fit_sprintr$step3[[lambda1_min]]$coef[,lambda3_min][(num_p+1):(num_p+num_inter)]
                                for(i in (1:num_inter)){
                                  index1 <- fit_sprintr$step2[[lambda1_min]][i,1]
                                  index1_name<-paste0("X",index1)
                                  index2 <- fit_sprintr$step2[[lambda1_min]][i,2]
                                  index2_name<-paste0("X",index2)
                                  if (strtoi(index1) < strtoi(index2)){
                                    inter_names<-append(inter_names,paste(index1_name,index2_name, sep=":"))
                                  }else if (strtoi(index1) > strtoi(index2)){
                                    inter_names<-append(inter_names,paste(index2_name,index1_name, sep=":"))
                                  }else{# if index1==index2, X1:X1
                                    inter_value<-inter_value[-i]
                                  }
                                }

                                for(i in (1:(length(inter_value)))){
                                  name<-inter_names[i]
                                  beta_Matrix[name,]<-inter_value[i]
                                }# put value into the beta

                                nz_coef<-as.matrix(beta_Matrix[beta_Matrix[, 1] != 0, ])
                                active<-rownames(as.matrix(nz_coef[-1,]))

                                y_valid_hat <- predict(fit_sprintr,draw[["x_valid"]])[[lambda1_min]][,lambda3_min]
                                mse_valid<-mean((draw[["y_valid"]] - y_valid_hat)^2)

                                return(list(beta = beta_Matrix,
                                            vnames = draw[["vnames"]],
                                            nonzero_coef = nz_coef,
                                            active = active,
                                            not_active = setdiff(draw[["vnames"]], active),
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

xyz_split <- new_method("xyz", "xyz",
                            method = function(model, draw) {
                              
                               tryCatch({
                                xyz_fit <- xyz_regression(draw[["x_train"]],draw[["y_train"]],n_lambda=10,alpha=0.9,L=10)
                                n_p<-dim(draw[["x_train"]])[2]
                                
                                y_test_hat<-predict(xyz_fit,draw[["x_test"]],10,n_p) # 30 = value_p
                                msetest <- colMeans((draw[["y_test"]] - y_test_hat)^2)   
                                lambda.min.index <- as.numeric(which.min(msetest[2:10]))
                                # lambda.min.index ==1 means select nothing, so we don't consider the case
                                lambda.min.index <- lambda.min.index+1
                                # since lmabda.min.index select from range 2:10, the index position need +1 to correspond the true index for fit.
                                main_list<-xyz_fit[[1]][lambda.min.index][[1]]
                                main_name_list<-paste0("X",main_list)
                                
                                active_name_list <- c() # active variable name
                                beta_Matrix <- matrix(0,nrow=length(draw[["vnames"]])+1,ncol=1) # variable matrix with name and value
                                rowname_list <- append("(Intercept)",draw[["vnames"]])
                                rownames(beta_Matrix) <- rowname_list
                                colnames(beta_Matrix) <- "1"
                                beta_Matrix[1,] <- xyz_fit[[6]][lambda.min.index]
                                
                                for (i in 1:length(main_name_list)){
                                  beta_Matrix[main_name_list[i],] <- xyz_fit[[2]][lambda.min.index][[1]][i]
                                }
                                
                                n_inters<-length(xyz_fit[[4]][lambda.min.index][[1]])
                                inter_names<-c()
                                inter_values<-xyz_fit[[4]][lambda.min.index][[1]]
                                
                                for (i in 1:n_inters){
                                  index1<-xyz_fit[[3]][lambda.min.index][[1]][,i][1]
                                  index1_name<-paste0("X",index1)
                                  index2<-xyz_fit[[3]][lambda.min.index][[1]][,i][2]
                                  index2_name<-paste0("X",index2)
                                  if(index1 < index2){
                                    inter<-paste(index1_name,index2_name,sep=":")
                                    inter_names<-append(inter_names,inter)
                                  }else if (index1 > index2){
                                    inter<-paste(index2_name,index1_name,sep=":")
                                    inter_names<-append(inter_names,inter)
                                  }else { # if index1 == index2
                                    inter_values<-inter_values[-i]
                                  }
                                }
                                
                                for(i in (1: length(inter_values))){
                                  beta_Matrix[inter_names[i],]<-inter_values[i]
                                }
                                
                                nz_coef <- as.matrix(beta_Matrix[beta_Matrix[, 1] != 0, ])
                                active<-rownames(as.matrix(nz_coef[-1,]))
                                
                                y_valid_hat<-predict(xyz_fit,draw[["x_valid"]],10,n_p)[,lambda.min.index]
                                mse_valid <- mean((draw[["y_valid"]] - y_valid_hat)^2)
                                
                                return(list(beta = beta_Matrix,
                                            vnames = draw[["vnames"]],
                                            nonzero_coef = nz_coef,
                                            active = active,
                                            not_active = setdiff(draw[["vnames"]], active),
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

ramp_split <- new_method("ramp", "ramp",
                        method = function(model, draw) {
                          
                          tryCatch({
                            n_p<-ncol(draw[["x_train"]])
                            
                            fit_ramp<-RAMP(draw[["x_train"]],draw[["y_train"]],family = "gaussian", penalty = "LASSO", gamma = NULL,
                                           inter = TRUE, hier = "Strong", eps = 1e-15, tune = "EBIC",
                                           penalty.factor = rep(1, n_p), inter.penalty.factor = 1, 
                                           max.iter = 100, n.lambda = 100,
                                           ebic.gamma = 1, refit = TRUE, trace = FALSE)
                            y_test_hat<-predict(fit_ramp,draw[["x_test"]],type ="response",allpath=T)
                            
                            mse_min<-10000
                            min_index<- -1
                           
                            # the last one usually be NA since is the criteria for stop, so skip to loop on the last one.
                            for (i in (1:(length(y_test_hat[1,])-1))){
                              mse<-mean((y_test_hat[,i]-draw[["y_test"]])^2)
                              
                              if (mse < mse_min){
                                mse_min <- mse
                                min_index<-i
                              }
                            }
                            
                            temp<-fit_ramp$interInd.list[min_index][[1]]
                            inter_name<-c()
                            inter_values<-fit_ramp$beta.i.mat[[min_index]]
      
                            beta_Matrix<-matrix(0,nrow=length(draw[["vnames"]])+1,ncol=1)
                            rowname_list<-append("(Intercept)",draw[["vnames"]])
                            rownames(beta_Matrix) <- rowname_list
                            colnames(beta_Matrix) <- "1"
                            beta_Matrix[1,]<-fit_ramp$a0.list[min_index]
                            if (length(temp)>0){
                              for (i in 1:(length(temp))){
                                index1<-strsplit(temp[i] , "X")[[1]][2]
                                index1_name <- paste0("X",index1)
                                index2<-strsplit(temp[i] , "X")[[1]][3]
                                index2_name <- paste0("X",index2)
                                if(strtoi(index1) < strtoi(index2)){
                                  inter_name<-append(inter_name,paste(index1_name,index2_name, sep=":"))
                                }else if (strtoi(index1) > strtoi(index2)){
                                  inter_name<-append(inter_name,paste(index2_name,index1_name, sep=":"))
                                }else {
                                  inter_values <- inter_values[-i]
                                }
                              }
                            }
                            
                            #main
                            main_value<-fit_ramp$beta.m.mat[,min_index]
                            if (length(main_value)>0){
                              for(i in 1:n_p){
                                #a<-main_value[i]
                                beta_Matrix[(i+1),]<-main_value[i]
                              }
                            }
                            #inter
                            if (length(temp)>0 & (length(inter_values))>0 &(length(inter_name))>0 ){
                              for (i in 1:length(inter_values)){
                                if(inter_name[i] %in% rownames(beta_Matrix)){
                                  beta_Matrix[inter_name[i],]<-inter_values[i]
                                }     
                              }
                            }
                            nz_coef<-as.matrix(beta_Matrix[beta_Matrix[, 1] != 0, ])
                            active<-rownames(as.matrix(nz_coef[-1,]))
                            
                            y_valid_hat<-predict(fit_ramp,draw[["x_valid"]],type ="response",allpath=T)[,min_index]
                            mse_valid<-mean((y_valid_hat-draw[["y_valid"]])^2)
                      
                            return(list(beta = beta_Matrix,
                                        vnames = draw[["vnames"]],
                                        nonzero_coef = nz_coef,
                                        active = active,
                                        not_active = setdiff(draw[["vnames"]], active),
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

family_split <- new_method("family", "family",
                         method = function(model, draw) {
                           
                           tryCatch({
                           alphas<- c(0.01,0.5,0.99)
                           lambdas<- seq(0.1,0.3,length = 50)
                             fit_family<- FAMILY(draw[["x_train"]], draw[["x_train"]], draw[["y_train"]], lambdas,
                                                 alphas,quad = TRUE,iter=500, verbose = TRUE )
                             
                             yhat<- predict(fit_family, draw[["x_test"]], draw[["x_test"]])
                             mse_hSH <-apply(yhat,c(2,3), "-" ,draw[["y_test"]])
                             mse_hSH<- apply(mse_hSH^2,c(2,3),sum)
                             im<- which(mse_hSH==min(mse_hSH),TRUE)
                             
                             coef<-coef(fit_family,XequalZ = TRUE)[[im[2]]][[im[1]]]
                             
                             # opt_value_alpha<-alphas[1]
                             # opt_value_lambdas<-lambdas[1]
                             # opt_index_alpha<- 1
                             # opt_index_lambdas<- 1
                             # mse_min <- 99999999
                             # 
                             # for(i in (1:length(alphas))){
                             #   for (j in (1:length(lambdas))){
                             #     mse<-mean(yhat[,j,i]-draw[["y_test"]])^2
                             #     if (mse < mse_min){
                             #       mse_min<-mse
                             #       opt_value_alpha<-alphas[i]
                             #       opt_value_lambdas<-lambdas[j]
                             #       opt_index_alpha<-i
                             #       opt_index_lambdas<-j
                             #     }
                             #   }
                             # }
                             # coef<-coef(fit_family,XequalZ = TRUE)[[opt_index_alpha]][[opt_index_lambdas]]
                             
                             beta_Matrix<-matrix(0,nrow=length(draw[["vnames"]])+1,ncol=1)
                             rowname_list<-append("(Intercept)",draw[["vnames"]])
                             rownames(beta_Matrix) <- rowname_list
                             colnames(beta_Matrix) <- "1"
                             beta_Matrix[1,]<-coef$intercept
                             
                             n_main<-dim(coef$mains)[1]
                             n_inter<-dim(coef$interacts)[1]
                             
                             #main
                             if((!is.null(n_main)) & n_main>0){
                               for (i in (1:n_main)){
                                   each_name<-paste0("X",coef$mains[i,"X"])
                                   beta_Matrix[each_name,]<-coef$mains[i,"Coef. est"]
                                 
                               }
                             }
                             #inter
                             if(!is.null(n_inter)){
                               for (i in (1:n_inter)){
                                 index1<-coef$interacts[i,"X"]
                                 index1_name<-paste0("X",index1)
                                 index2<-coef$interacts[i,"Z"]
                                 index2_name<-paste0("X",index2)

                                 if(index1 < index2){
                                   inter_name<-paste(index1_name,index2_name, sep=":")
                                   beta_Matrix[each_name,]<-coef$interacts[i,"Coef. est"]
                                 }
                               }
                             }
                             
                             nz_coef <- as.matrix(beta_Matrix[beta_Matrix[, 1] != 0, ])
                             active<-rownames(as.matrix(nz_coef[-1,]))
                             
                             y_valid_hat<- predict(fit_family, draw[["x_valid"]], draw[["x_valid"]])[,im[1],im[2]]
                             mse_valid <- mean((draw[["y_valid"]] - y_valid_hat)^2)
                             
                             return(list(beta = beta_Matrix,
                                         vnames = draw[["vnames"]],
                                         nonzero_coef = nz_coef,
                                         active = active,
                                         not_active = setdiff(draw[["vnames"]], active),
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

# # All pairs Method For Binary response-------------------


family_binary_split <- new_method("family_binary", "family_binary",
                           method = function(model, draw) {
                             
                             tryCatch({
                               alphas<- c(0.01,0.5,0.99)
                               lambdas<- seq(0.1,0.3,length = 50)
                               fit_family<- FAMILY(draw[["x_train"]], draw[["x_train"]], draw[["y_train"]],family = "binomial", lambdas,
                                                   alphas,quad = TRUE,iter=500, verbose = TRUE )
                               
                               yhat<- predict(fit_family, draw[["x_test"]], draw[["x_test"]])
                               mse_hSH <-apply(yhat,c(2,3), "-" ,draw[["y_test"]])
                               mse_hSH<- apply(mse_hSH^2,c(2,3),sum)
                               im<- which(mse_hSH==min(mse_hSH),TRUE)
                               
                               coef<-coef(fit_family,XequalZ = TRUE)[[im[2]]][[im[1]]]
                               
                            
                               beta_Matrix<-matrix(0,nrow=length(draw[["vnames"]])+1,ncol=1)
                               rowname_list<-append("(Intercept)",draw[["vnames"]])
                               rownames(beta_Matrix) <- rowname_list
                               colnames(beta_Matrix) <- "1"
                               beta_Matrix[1,]<-coef$intercept
                               
                               n_main<-dim(coef$mains)[1]
                               n_inter<-dim(coef$interacts)[1]
                               
                               #main
                               if((!is.null(n_main)) & n_main>0){
                                 for (i in (1:n_main)){
                                   each_name<-paste0("X",coef$mains[i,"X"])
                                   beta_Matrix[each_name,]<-coef$mains[i,"Coef. est"]
                                   
                                 }
                               }
                               #inter
                               if(!is.null(n_inter)){
                                 for (i in (1:n_inter)){
                                   index1<-coef$interacts[i,"X"]
                                   index1_name<-paste0("X",index1)
                                   index2<-coef$interacts[i,"Z"]
                                   index2_name<-paste0("X",index2)
                                   
                                   if(index1 < index2){
                                     inter_name<-paste(index1_name,index2_name, sep=":")
                                     beta_Matrix[each_name,]<-coef$interacts[i,"Coef. est"]
                                   }
                                 }
                               }
                               
                               nz_coef <- as.matrix(beta_Matrix[beta_Matrix[, 1] != 0, ])
                               active<-rownames(as.matrix(nz_coef[-1,]))
                               
                               y_valid_hat<- predict(fit_family, draw[["x_valid"]], draw[["x_valid"]])[,im[1],im[2]]
                               y_valid_hat_b<- rbinom(length(y_valid_hat), size = 1, prob = y_valid_hat) # Convert prob to binary
                               mse_valid <- mean((draw[["y_valid"]] - y_valid_hat_b)^2)
                               
                               return(list(beta = beta_Matrix,
                                           vnames = draw[["vnames"]],
                                           nonzero_coef = nz_coef,
                                           active = active,
                                           not_active = setdiff(draw[["vnames"]], active),
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

ramp_binary_split <- new_method("ramp_binary", "ramp_binary",
                         method = function(model, draw) {
                           
                            tryCatch({
                             n_p<-ncol(draw[["x_train"]])
                             
                             fit_ramp<-RAMP(draw[["x_train"]],draw[["y_train"]],family = "binomial", penalty = "LASSO", gamma = NULL,
                                            inter = TRUE, hier = "Strong", eps = 1e-15, tune = "EBIC",
                                            penalty.factor = rep(1, n_p), inter.penalty.factor = 1, 
                                            max.iter = 100, n.lambda = 100,
                                            ebic.gamma = 1, refit = TRUE, trace = FALSE)
                             y_test_hat<-predict(fit_ramp,draw[["x_test"]],type ="response",allpath=T)
                             
                             mse_min<-10000
                             min_index<- -1
                             
                             # the last one usually be NA since is the criteria for stop, so skip to loop on the last one.
                             for (i in (1:(length(y_test_hat[1,])-1))){
                               mse<-mean((y_test_hat[,i]-draw[["y_test"]])^2)
                               if (is.na(mse)){
                                 print(i)
                               }
                               else if (mse < mse_min){
                                 mse_min <- mse
                                 min_index<-i
                               }
                             }
                             
                             temp<-fit_ramp$interInd.list[min_index][[1]]
                             inter_name<-c()
                             inter_values<-fit_ramp$beta.i.mat[[min_index]]
                             
                             beta_Matrix<-matrix(0,nrow=length(draw[["vnames"]])+1,ncol=1)
                             rowname_list<-append("(Intercept)",draw[["vnames"]])
                             rownames(beta_Matrix) <- rowname_list
                             colnames(beta_Matrix) <- "1"
                             beta_Matrix[1,]<-fit_ramp$a0.list[min_index]
                             if (length(temp)>0){
                               for (i in 1:(length(temp))){
                                 index1<-strsplit(temp[i] , "X")[[1]][2]
                                 index1_name <- paste0("X",index1)
                                 index2<-strsplit(temp[i] , "X")[[1]][3]
                                 index2_name <- paste0("X",index2)
                                 if(strtoi(index1) < strtoi(index2)){
                                   inter_name<-append(inter_name,paste(index1_name,index2_name, sep=":"))
                                 }else if (strtoi(index1) > strtoi(index2)){
                                   inter_name<-append(inter_name,paste(index2_name,index1_name, sep=":"))
                                 }else {
                                   inter_values <- inter_values[-i]
                                 }
                               }
                             }
                             
                             #main
                             main_value<-fit_ramp$beta.m.mat[,min_index]
                             if (length(main_value)>0){
                               for(i in 1:n_p){
                                 #a<-main_value[i]
                                 beta_Matrix[(i+1),]<-main_value[i]
                               }
                             }
                             #inter
                             if (length(temp)>0 & (length(inter_values))>0 &(length(inter_name))>0 ){
                               for (i in 1:length(inter_values)){
                                 if(inter_name[i] %in% rownames(beta_Matrix)){
                                   beta_Matrix[inter_name[i],]<-inter_values[i]
                                 }     
                               }
                             }
                             nz_coef<-as.matrix(beta_Matrix[beta_Matrix[, 1] != 0, ])
                             active<-rownames(as.matrix(nz_coef[-1,]))
                             
                             y_valid_hat<-predict(fit_ramp,draw[["x_valid"]],type ="class",allpath=T)[,min_index]
                             mse_valid<-mean((as.integer(y_valid_hat)-draw[["y_valid"]])^2)
                             
                             return(list(beta = beta_Matrix,
                                         vnames = draw[["vnames"]],
                                         nonzero_coef = nz_coef,
                                         active = active,
                                         not_active = setdiff(draw[["vnames"]], active),
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


glinternet_binary_split_XX <- new_method("glinternet_binary_XX", "glinternet binary XX",
                                 method = function(model, draw) {
                                   
                                   tryCatch({
                                     fitGL <- glinternet(X = draw[["x_train"]], Y = draw[["y_train"]],
                                                         numLevels = rep(1, ncol(draw[["x_train"]])),
                                                         nLambda = 100, family = "binomial",
                                                         verbose = F)
                                     
                                     y_test_hat <- predict(fitGL, X = draw[["x_test"]])
                                     mse_test <- colMeans((draw[["y_test"]] - y_test_hat)^2)
                                     lambda.min.index <- as.numeric(which.min(mse_test))
                                     lambda.min <- fitGL$lambda[which.min(mse_test)]
                                     
                                     # y_valid_hat <- predict(fitGL, X = draw[["x_valid"]], lambda = lambda.min)
                                     # mse_valid <- mean((draw[["y_valid"]] - drop(y_valid_hat))^2)
                                  
                                     # 
                                     fit_coef <- coef(fitGL, lambdaIndex = lambda.min.index)
                                    
                                     
                                     y_valid_hat <- predict(fitGL, X = draw[["x_valid"]], lambda = lambda.min)
                                     y_valid_hat_b<- rbinom(length(y_valid_hat), size = 1, prob = y_valid_hat) # Convert prob to binary
                                     mse_valid <- mean((draw[["y_valid"]] - y_valid_hat_b)^2)
                                     
                                     
                                     mains <- colnames(draw[["x_train"]])[fit_coef[[1]]$mainEffects$cont]
                                     inters <- paste(colnames(draw[["x_train"]])[fit_coef[[1]]$interactions$contcont[,1]],
                                                     colnames(draw[["x_train"]])[fit_coef[[1]]$interactions$contcont[,2]], sep=":")
                                     
                                     beta_Matrix<-matrix(0,nrow=length(draw[["vnames"]])+1,ncol=1)
                                     rowname_list<-append("(Intercept)",draw[["vnames"]]) # keep the vnames order the same as other methods
                                     rownames(beta_Matrix) <- rowname_list
                                     colnames(beta_Matrix) <- "1"
                                     beta_Matrix[1,]<- 0
                                     
                                     for (i in (1:length(mains))){
                                       tryCatch({
                                         beta_Matrix[mains[i],]<-fit_coef[[1]]$mainEffectsCoef$cont[[i]]
                                       },
                                       error = function(err) {
                                         #continue
                                       })
                                     }
                                     
                                     for (i in (1:length(inters))){
                                       tryCatch({
                                         beta_Matrix[inters[i],]<-fit_coef[[1]]$interactionsCoef$contcont[[i]]
                                       },
                                       error = function(err) {
                                         #continue
                                       })
                                     }
                                     
                                     return(list(beta = beta_Matrix,
                                                 vnames = draw[["vnames"]],
                                                 nonzero_coef = as.matrix(beta_Matrix[beta_Matrix[, 1] != 0, ]),
                                                 active = c(mains, inters),
                                                 not_active = setdiff(draw[["vnames"]], c(mains, inters)),
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

hiernet_binary_split <- new_method("hinernet_binary", "hinernet_binary",
                            method = function(model, draw) {
                              
                              tryCatch({
                                
                                fit_hierNet_path<-hierNet.logistic.path(draw[["x_train"]], draw[["y_train"]],strong=TRUE)
                                
                                ytest_hat <- predict(fit_hierNet_path, newx = draw[["x_test"]])
                                col_msetest<-c()
                                
                                for (i in 1:length(ytest_hat$prob[1,])){
                                  msetest <- mean((draw[["y_test"]] - ytest_hat$prob[,i])^2)
                                  col_msetest<-append(col_msetest,msetest)
                                }
                                lambda_min_index<-which.min(col_msetest)
                                lambda_min<-fit_hierNet_path$lamlist[which.min(col_msetest)]
                                
                                # bp:p-vector of estimated "positive part" main effect (p=# features)
                                # bn:p-vector of estimated "negative part" main effect (p=# features)
                                # bp-bn:overall main effect estimated
                                beta_main<-fit_hierNet_path$bp[,lambda_min_index]-fit_hierNet_path$bn[,lambda_min_index] #overall main effect estimated
                                as.vector(beta_main)
                                
                                beta_Matrix<-matrix(0,nrow=length(draw[["vnames"]])+1,ncol=1)
                                rowname_list<-append("(Intercept)",draw[["vnames"]])
                                rownames(beta_Matrix) <- rowname_list
                                colnames(beta_Matrix) <- "1"
                                beta_Matrix[1,]<-0
                                
                                num_p<-dim(draw[["x_train"]])[2]
                                for (i in (1:num_p)){
                                  beta_Matrix[i+1]<-beta_main[i]
                                }
                                
                                inter<-c()
                                
                                for (i in (1:(num_p-1))){
                                  for (j in ((i+1):num_p)){
                                    inter<-append(inter, fit_hierNet_path$th[,,lambda_min_index][i,j])
                                    
                                  }
                                }
                                
                                beta_Matrix[(num_p+2):(length(beta_Matrix))]<-inter
                                
                                nzcoef<-as.matrix(beta_Matrix[beta_Matrix[, 1] != 0, ])
                                
                                active<-rownames(nzcoef)
                                
                                y_valid_hat <- predict(fit_hierNet_path, newx = draw[["x_valid"]])$yhat[,lambda_min_index]
                                mse_valid <- mean((draw[["y_valid"]] - y_valid_hat)^2)
                                
                                return(list(beta = beta_Matrix,
                                            vnames = draw[["vnames"]],
                                            nonzero_coef = nzcoef,
                                            active = active,
                                            not_active = setdiff(draw[["vnames"]], active),
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



