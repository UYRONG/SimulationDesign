## @knitr models

# n should be the total of train and test. 2*n is chosen as the validation
make_gendata_Paper_data_split <- function(n, p, corr, betaE, SNR, lambda.type, parameterIndex) {
  
  main <- paste0("X", seq_len(p))
  vnames <- c(main, "E", paste0(main,":E"))
  # used for glmnet and lasso backtracking
  # f.identity <- function(i) i
  
  new_model(name = "gendata_thesis_split_v4",
            label = sprintf("n = %s, p = %s, corr = %s, betaE = %s, SNR = %s, index = %s, lambda = %s",
                            n, p, corr, betaE, SNR, parameterIndex, lambda.type),
            params = list(n = n, p = p, corr = corr, betaE = betaE, SNR = SNR,
                          lambda.type = lambda.type, parameterIndex = parameterIndex),
            simulate = function(n, p, corr, betaE, SNR, parameterIndex, lambda.type, nsim) {
              models <- list()
              for(i in seq(nsim)) {
                # test and training set
                DT <- sail::gendata(n = n, p = p, corr = corr, betaE = betaE,
                                    SNR = SNR, parameterIndex = parameterIndex)
                
                #validation set
                DT_val <- sail::gendata(n = 2*n, p = p, corr = corr, betaE = betaE,
                                        SNR = SNR, parameterIndex = parameterIndex)
                
                trainIndex <- drop(caret::createDataPartition(DT$y, p = 0.5, list = FALSE, times = 1))
                testIndex <- setdiff(seq_along(DT$y), trainIndex)
                
                xtrain <- DT$x[trainIndex,,drop=FALSE]
                xtest <- DT$x[testIndex,,drop=FALSE]
                xvalid <- DT_val$x
                
                etrain <- DT$e[trainIndex]
                etest <- DT$e[testIndex]
                evalid <- DT_val$e
                
                ytrain <- DT$y[trainIndex]
                ytest <- DT$y[testIndex]
                yvalid <- DT_val$y
                
                main <- colnames(DT$x)
                vnames <- c(main, "E", paste0(main,":E"))
                vnames_lasso <- c("E", main) # needs to be in this order for glinternet
                
                # as is done in pliable lasso, only feed (E,X) to glmnet
                X_main <- cbind(E = etrain, xtrain)
                
                # test set
                X_main_test <- cbind(E = etest, xtest)
                
                # validation set
                X_main_valid <- cbind(E = evalid, xvalid)
                
                models[[i]] <- list(xtrain = xtrain, etrain = etrain, ytrain = ytrain, xtrain_lasso = X_main,
                                    xtest = xtest, etest = etest, ytest = ytest, xtest_lasso = X_main_test,
                                    xvalid = xvalid, evalid = evalid, yvalid = yvalid, xvalid_lasso = X_main_valid,
                                    causal = DT$causal, not_causal = DT$not_causal,
                                    vnames = vnames, vnames_lasso = vnames_lasso)
              }
              return(models)
            })
  
}

# Simulation Case 1 : All pairs, Continuous Outcome
generate_data_split_case_1 <- function(n, p, corr, SNR, lambda.type, parameterIndex) {
  
  main <- paste0("X", seq_len(p))
  vnames <- c(main, "E", paste0(main,":E"))
  
  new_model(name = "gendata_thesis_split_v4",
            label = sprintf("n = %s, p = %s, corr = %s, SNR = %s, index = %s, lambda = %s",
                            n, p, corr, SNR, parameterIndex, lambda.type),
            params = list(n = n, p = p, corr = corr, SNR = SNR,
                          lambda.type = lambda.type, parameterIndex = parameterIndex),
            simulate = function(n, p, corr, SNR, parameterIndex, lambda.type, nsim) {
              models <- list()
              for(i in seq(nsim)) {
                # test and training set
                DT <- sail::gendata(n = n, p = p, corr = corr,
                                    SNR = SNR, parameterIndex = parameterIndex)
                
                #validation set
                DT_val <- sail::gendata(n = 2*n, p = p, corr = corr, 
                                        SNR = SNR, parameterIndex = parameterIndex)
                
                trainIndex <- drop(caret::createDataPartition(DT$y, p = 0.5, list = FALSE, times = 1))
                testIndex <- setdiff(seq_along(DT$y), trainIndex)
                
                xtrain <- DT$x[trainIndex,,drop=FALSE]
                xtest <- DT$x[testIndex,,drop=FALSE]
                xvalid <- DT_val$x
                
                etrain <- DT$e[trainIndex]
                etest <- DT$e[testIndex]
                evalid <- DT_val$e
                
                ytrain <- DT$y[trainIndex]
                ytest <- DT$y[testIndex]
                yvalid <- DT_val$y
                
                main <- colnames(DT$x)
                vnames <- c(main, "E", paste0(main,":E"))
                vnames_lasso <- c("E", main) # needs to be in this order for glinternet
                
                # as is done in pliable lasso, only feed (E,X) to glmnet
                X_main <- cbind(E = etrain, xtrain)
                
                # test set
                X_main_test <- cbind(E = etest, xtest)
                
                # validation set
                X_main_valid <- cbind(E = evalid, xvalid)
                
                models[[i]] <- list(xtrain = xtrain, etrain = etrain, ytrain = ytrain, xtrain_lasso = X_main,
                                    xtest = xtest, etest = etest, ytest = ytest, xtest_lasso = X_main_test,
                                    xvalid = xvalid, evalid = evalid, yvalid = yvalid, xvalid_lasso = X_main_valid,
                                    causal = DT$causal, not_causal = DT$not_causal,
                                    vnames = vnames, vnames_lasso = vnames_lasso)
              }
              return(models)
            })
  
}

# Simulation Case 2 : All pairs, Binary Outcome
generate_data_split_case_2 <- function(n, p, corr, betaE, SNR, lambda.type, parameterIndex) {
  
  main <- paste0("X", seq_len(p))
  vnames <- c(main, "E", paste0(main,":E"))
  # used for glmnet and lasso backtracking
  # f.identity <- function(i) i
  
  new_model(name = "gendata_thesis_split_v4",
            label = sprintf("n = %s, p = %s, corr = %s, betaE = %s, SNR = %s, index = %s, lambda = %s",
                            n, p, corr, betaE, SNR, parameterIndex, lambda.type),
            params = list(n = n, p = p, corr = corr, betaE = betaE, SNR = SNR,
                          lambda.type = lambda.type, parameterIndex = parameterIndex),
            simulate = function(n, p, corr, betaE, SNR, parameterIndex, lambda.type, nsim) {
              models <- list()
              for(i in seq(nsim)) {
                # test and training set
                DT <- sail::gendata(n = n, p = p, corr = corr, betaE = betaE,
                                    SNR = SNR, parameterIndex = parameterIndex)
                
                #validation set
                DT_val <- sail::gendata(n = 2*n, p = p, corr = corr, betaE = betaE,
                                        SNR = SNR, parameterIndex = parameterIndex)
                
                trainIndex <- drop(caret::createDataPartition(DT$y, p = 0.5, list = FALSE, times = 1))
                testIndex <- setdiff(seq_along(DT$y), trainIndex)
                
                xtrain <- DT$x[trainIndex,,drop=FALSE]
                xtest <- DT$x[testIndex,,drop=FALSE]
                xvalid <- DT_val$x
                
                etrain <- DT$e[trainIndex]
                etest <- DT$e[testIndex]
                evalid <- DT_val$e
                
                ytrain <- DT$y[trainIndex]
                ytest <- DT$y[testIndex]
                yvalid <- DT_val$y
                
                main <- colnames(DT$x)
                vnames <- c(main, "E", paste0(main,":E"))
                vnames_lasso <- c("E", main) # needs to be in this order for glinternet
                
                # as is done in pliable lasso, only feed (E,X) to glmnet
                X_main <- cbind(E = etrain, xtrain)
                
                # test set
                X_main_test <- cbind(E = etest, xtest)
                
                # validation set
                X_main_valid <- cbind(E = evalid, xvalid)
                
                models[[i]] <- list(xtrain = xtrain, ytrain = ytrain, xtrain_lasso = X_main,
                                    xtest = xtest, ytest = ytest, xtest_lasso = X_main_test,
                                    xvalid = xvalid, yvalid = yvalid, xvalid_lasso = X_main_valid,
                                    causal = DT$causal, not_causal = DT$not_causal,
                                    vnames = vnames, vnames_lasso = vnames_lasso)
              }
              return(models)
            })
  
}
# Simulation Case 3 : GXE, Continuous Outcome, n=100, p=30
generate_data_split_case_3 <- function(n, p, corr, betaE, SNR, lambda.type, parameterIndex) {
  
  main <- paste0("X", seq_len(p))
  vnames <- c(main, "E", paste0(main,":E"))
  # used for glmnet and lasso backtracking
  # f.identity <- function(i) i
  
  new_model(name = "gendata_thesis_split_v4",
            label = sprintf("n = %s, p = %s, corr = %s, betaE = %s, SNR = %s, index = %s, lambda = %s",
                            n, p, corr, betaE, SNR, parameterIndex, lambda.type),
            params = list(n = n, p = p, corr = corr, betaE = betaE, SNR = SNR,
                          lambda.type = lambda.type, parameterIndex = parameterIndex),
            simulate = function(n, p, corr, betaE, SNR, parameterIndex, lambda.type, nsim) {
              models <- list()
              for(i in seq(nsim)) {
                # test and training set
                DT <- sail::gendata(n = n, p = p, corr = corr, betaE = betaE,
                                    SNR = SNR, parameterIndex = parameterIndex)
                
                #validation set
                DT_val <- sail::gendata(n = 2*n, p = p, corr = corr, betaE = betaE,
                                        SNR = SNR, parameterIndex = parameterIndex)
                
                trainIndex <- drop(caret::createDataPartition(DT$y, p = 0.5, list = FALSE, times = 1))
                testIndex <- setdiff(seq_along(DT$y), trainIndex)
                
                xtrain <- DT$x[trainIndex,,drop=FALSE]
                xtest <- DT$x[testIndex,,drop=FALSE]
                xvalid <- DT_val$x
                
                etrain <- DT$e[trainIndex]
                etest <- DT$e[testIndex]
                evalid <- DT_val$e
                
                ytrain <- DT$y[trainIndex]
                ytest <- DT$y[testIndex]
                yvalid <- DT_val$y
                
                main <- colnames(DT$x)
                vnames <- c(main, "E", paste0(main,":E"))
                vnames_lasso <- c("E", main) # needs to be in this order for glinternet
                
                # as is done in pliable lasso, only feed (E,X) to glmnet
                X_main <- cbind(E = etrain, xtrain)
                
                # test set
                X_main_test <- cbind(E = etest, xtest)
                
                # validation set
                X_main_valid <- cbind(E = evalid, xvalid)
                
                models[[i]] <- list(xtrain = xtrain, etrain = etrain, ytrain = ytrain, xtrain_lasso = X_main,
                                    xtest = xtest, etest = etest, ytest = ytest, xtest_lasso = X_main_test,
                                    xvalid = xvalid, evalid = evalid, yvalid = yvalid, xvalid_lasso = X_main_valid,
                                    causal = DT$causal, not_causal = DT$not_causal,
                                    vnames = vnames, vnames_lasso = vnames_lasso)
              }
              return(models)
            })
  
}

