## @knitr models
pacman::p_load(truncnorm)

# DEFUALT hierarchy = STRONG
gendata_XX <- function(n, p, corr = 0,
                         SNR = 2, hierarchy = c("strong", "weak", "none"),
                         nonlinear = TRUE, interactions = TRUE, causal, not_causal) {
  # this is modified from "VARIABLE SELECTION IN NONPARAMETRIC ADDITIVE MODEL" huang et al, Ann Stat.

  if (!requireNamespace("truncnorm", quietly = TRUE)) {
    stop("Package \"truncnorm\" needed for this function to simulate data. Please install it.",
         call. = FALSE
    )
  }
  
  hierarchy <- match.arg(hierarchy)
  
  # covariates
  W <- replicate(n = p, truncnorm::rtruncnorm(n, a = 0, b = 1))
  U <- truncnorm::rtruncnorm(n, a = 0, b = 1)
  V <- truncnorm::rtruncnorm(n, a = 0, b = 1)
  
  X1 <- (W[, 1] + corr * U) / (1 + corr)
  X2 <- (W[, 2] + corr * U) / (1 + corr)
  X3 <- (W[, 3] + corr * U) / (1 + corr)
  X4 <- (W[, 4] + corr * U) / (1 + corr)
  
  X <- (W[, 5:p] + corr * V) / (1 + corr)
  
  Xall <- cbind(X1, X2, X3, X4, X)
  
  colnames(Xall) <- paste0("X", seq_len(p))
  
  # see "Variable Selection in NonParametric Addditive Model" Huang Horowitz and Wei

    f1 <- function(x)  5 * x
    f2 <- function(x)  3 * (x + 1)
    f3 <- function(x)  4 * x
    f4 <- function(x)  6 * (x - 2)

  # error
  error <- stats::rnorm(n)
  
  Y.star <- f1(X1) +
    f2(X2) +
    f3(X3) +
    f4(X4) 

  scenario <- "new"
  
  k <- sqrt(stats::var(Y.star) / (SNR * stats::var(error)))
  
  Y <- Y.star + as.vector(k) * error
  
  return(list(
    x = Xall, y = Y, Y.star = Y.star, f1 = f1(X1),
    f2 = f2(X2), f3 = f3(X3), f4 = f4(X4), 
    f1.f = f1, f2.f = f2, f3.f = f3, f4.f = f4,
    X1 = X1, X2 = X2, X3 = X3, X4 = X4, scenario = scenario,
    causal = causal, not_causal = not_causal
  ))
}

generate_data <- function(n, p, corr, SNR) {
  if (!requireNamespace("truncnorm", quietly = TRUE)) {
    stop("Package \"truncnorm\" needed for this function to simulate data. Please install it.",
         call. = FALSE
    )
  }
  
  main <- paste0("X", seq_len(p))
  vnames <- c(main, "E", paste0(main, ":E"))
  
  hierarchy <- "strong"
  nonlinear <- FALSE
  interactions <- FALSE
  causal <- c("X1", "X2", "X3", "X4")
  
  not_causal <- setdiff(vnames, causal)
  
  DT <- gendata_XX(
    n = n, p = p, corr = corr, SNR = SNR,
    hierarchy = hierarchy, nonlinear = nonlinear, interactions = interactions,
    causal = causal, not_causal = not_causal
  )
  return(DT)
}





# n should be the total of train and test. 2*n is chosen as the validation------------------------------------------------------------------------
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

# Simulation Case 1 : All pairs, Continuous Outcome------------------------------------------------------------------------
generate_data_split_case_1 <- function(n, p, corr, SNR, lambda.type, parameterIndex = 1) {
  
  main <- paste0("X", seq_len(p))
  vnames <- c(main, "E", paste0(main,":E"))
  
  new_model(name = "generate_data_split_case_1",
            label = sprintf("n = %s, p = %s, corr = %s, SNR = %s, index = %s, lambda = %s",
                            n, p, corr, SNR, parameterIndex, lambda.type),
            params = list(n = n, p = p, corr = corr, SNR = SNR,
                          lambda.type = lambda.type, parameterIndex = parameterIndex),
            simulate = function(n, p, corr, SNR, parameterIndex, lambda.type, nsim = 1) {
              models <- list()
              for(i in seq(nsim)) {
                # test and training set
                DT <- generate_data(n = n, p = p, corr = corr,
                                    SNR = SNR)
                
                #validation set
                DT_val <- generate_data(n = 2*n, p = p, corr = corr, 
                                        SNR = SNR)
                
                trainIndex <- drop(caret::createDataPartition(DT$y, p = 0.5, list = FALSE, times = 1))
                testIndex <- setdiff(seq_along(DT$y), trainIndex)
                
                xtrain <- DT$x[trainIndex,,drop=FALSE]
                xtest <- DT$x[testIndex,,drop=FALSE]
                xvalid <- DT_val$x
                
                ytrain <- DT$y[trainIndex]
                ytest <- DT$y[testIndex]
                yvalid <- DT_val$y
                
                main <- colnames(DT$x)
                vnames <- c(main, "E", paste0(main,":E"))
                vnames_lasso <- c("E", main) # needs to be in this order for glinternet
                
                # as is done in pliable lasso, only feed (E,X) to glmnet
                X_main <- cbind(xtrain)
                
                # test set
                X_main_test <- cbind(xtest)
                
                # validation set
                X_main_valid <- cbind( xvalid)
                
                models[[i]] <- list(xtrain = xtrain, ytrain = ytrain, xtrain_lasso = X_main,
                                    xtest = xtest,ytest = ytest, xtest_lasso = X_main_test,
                                    xvalid = xvalid,  yvalid = yvalid, xvalid_lasso = X_main_valid,
                                    causal = DT$causal, not_causal = DT$not_causal,
                                    vnames = vnames, vnames_lasso = vnames_lasso)
              }
              return(models)
            })
  
}

try <-generate_data_split_case_1 (100, 30, 0, 2, lambda.type = "lambda.min", parameterIndex=1)



sprintf("n = %s, p = %s, corr = %s, SNR = %s, index = %s, lambda = %s",
        100, 30, 0, 2, "1", "lambda.min")





# Simulation Case 2 : All pairs, Binary Outcome-------------------------------------------------------------------------------
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

