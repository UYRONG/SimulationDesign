## @knitr models
pacman::p_load(truncnorm)

# generate data for X_E case
y_function_XE<-function(n, p, corr=0, betaE = 2, SNR = 2, case){
  # covariates
  W <- replicate(n = p, truncnorm::rtruncnorm(n, a = 0, b = 1))
  U <- truncnorm::rtruncnorm(n, a = 0, b = 1)
  V <- truncnorm::rtruncnorm(n, a = 0, b = 1)
  
  E <- stats::rnorm(n)
  
  Y <- 0
  
  X1 <- (W[, 1] + corr * U) / (1 + corr)
  X2 <- (W[, 2] + corr * U) / (1 + corr)
  X3 <- (W[, 3] + corr * U) / (1 + corr)
  X4 <- (W[, 4] + corr * U) / (1 + corr)
  X5 <- (W[, 5] + corr * U) / (1 + corr)
  
  
  if (case==1){
    X <- (W[, 6:p] + corr * V) / (1 + corr)
    
    Xall <- cbind(X1, X2, X3, X4, X5, X)
    
    colnames(Xall) <- paste0("X", seq_len(p))
    
    f1 <- function(x) 3 * x
    f2 <- function(x) -3 * x
    f3 <- function(x) -3 * x
    f4 <- function(x) 3 * x
    f5 <- function(x) 3 * x
    f4.inter = function(x, e) 1.5 * e * x
    f5.inter = function(x, e) -1.5 * e * x
    
    # 5 main effects + 2 interaction
    Y <- f1(X1) + f2(X2) + f3(X3) + f4(X4) + f5(X5) + betaE*E + f4.inter(X4,E) + f5.inter(X5,E) + rnorm(n)
  }
  else{
    X6 <- (W[, 6] + corr * U) / (1 + corr)
    X7 <- (W[, 7] + corr * U) / (1 + corr)
    X8 <- (W[, 8] + corr * U) / (1 + corr)
    X9 <- (W[, 9] + corr * U) / (1 + corr)
    X10 <- (W[, 10] + corr * U) / (1 + corr)
    
    X <- (W[, 11:p] + corr * V) / (1 + corr)
    
    Xall <- cbind(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X)
    
    colnames(Xall) <- paste0("X", seq_len(p))
    
    f1 <- function(x) 3 * x
    f2 <- function(x) -3 * x
    f3 <- function(x) -3 * x
    f4 <- function(x) 3 * x
    f5 <- function(x) 3 * x
    f6 <- function(x) 3 * x
    f7 <- function(x) -3 * x
    f8 <- function(x) -3 * x
    f9 <- function(x) 3 * x
    f10 <- function(x) 3 * x
    
    f4.inter = function(x, e) 1.5 * e * x
    f5.inter = function(x, e) -1.5 * e * x
    f6.inter = function(x, e) 1.5 * e * x
    f7.inter = function(x, e) -1.5 * e * x
    f8.inter = function(x, e) 1.5 * e * x
    
    Y <- f1(X1) + f2(X2) + f3(X3) + f4(X4) + f5(X5) +
      f6(X6) + f7(X7) + f8(X8) + f9(X9) + f10(X10)+ 
      betaE*E + f4.inter(X4,E) + f5.inter(X5,E) +
      f6.inter(X6,E) + f7.inter(X7,E) +
      f8.inter(X8,E) + rnorm(n)
  }
  return(list(x = Xall, y = Y ,e = E))
}

generate_data_case_XE<-function(n, p, corr, betaE, SNR, case){
  # test and training set
  # X <- matrix(rnorm(n*p), nrow=n)
  # X <- scale(X,TRUE,TRUE)
  gendata <- y_function_XE(n, p, corr, betaE, SNR, case)
  
  split_value <- n*0.8
  m <- 2*n
  gendata_valid <- y_function_XE(m, p, corr, betaE, SNR, case)
  
  main <- paste0("X", seq_len(p))
  
  vnames_lasso <- c("E", main) # needs to be in this order for glinternet
  
  vnames<-c(main)
  vnames<-append(vnames, "E") # others
  
  for (i in (1:(length(main)))){
    conc<-paste(main[i], "E", sep = ":")
    vnames<-append(vnames, conc)
    vnames_lasso<-append(vnames_lasso, conc)
  }
  
  true_beta<-matrix(0,nrow=(length(vnames)),ncol=1)
  true_beta_list<-vnames
  rownames(true_beta) <- true_beta_list
  colnames(true_beta) <- "true value"
  
  if (case ==1){
    causal <- c("X1","X2","X3","X4","X5","E","X4:E","X5:E")
    true_beta[1:5,1]<-c(3,-3,-3,3,3)
    true_beta[(p+5):(p+6),1]<-c(1.5,-1.5)
    true_beta["E",1]<-betaE
  }else if (case ==2){
    causal <- c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","E","X4:E","X5:E","X6:E","X7:E","X8:E")
    true_beta[1:10,1]<-c(3,-3,-3,3,3,3,-3,-3,3,3)
    true_beta[(p+5):(p+9),]<-c(1.5,-1.5,1.5,-1.5,1.5)
    true_beta["E",1]<-betaE
  }else if (case ==3){
    causal <- c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","E","X4:E","X5:E","X6:E","X7:E","X8:E")
    true_beta[1:10,1]<-c(3,-3,-3,3,3,3,-3,-3,3,3)
    true_beta[(p+5):(p+9),]<-c(1.5,-1.5,1.5,-1.5,1.5)
    true_beta["E",1]<-betaE
  }
  not_causal <- setdiff(vnames, causal)
  
  xtrain <- gendata$x[1:split_value,] 
  ytrain = gendata$y[1:split_value]
  etrain = gendata$e[1:split_value]
  xtest = gendata$x[(split_value+1):n,]
  ytest = gendata$y[(split_value+1):n]
  etest = gendata$e[(split_value+1):n]
  xvalid = gendata_valid$x
  yvalid = gendata_valid$y
  evalid = gendata_valid$e
  
  # as is done in pliable lasso, only feed (E,X) to glmnet
  x_main <- cbind(E = etrain, xtrain)
  
  # test set
  x_main_test <- cbind(E = etest, xtest)
  
  # validation set
  x_main_valid <- cbind(E = evalid, xvalid)
  
  result<-list("x_train" = xtrain , "y_train" = ytrain, "e_train" = etrain, "x_train_lasso" = x_main,
               "x_test" = xtest, "y_test" = ytest,"e_test" = etest, "x_test_lasso" = x_main_test,
               "x_valid" = xvalid, "y_valid" = yvalid, "e_valid" = evalid, "x_valid_lasso" = x_main_valid,
               "vnames" = vnames, "vnames_lasso" = vnames_lasso, "betaE"= betaE, "true_beta" = true_beta,
               "causal" = causal , "not_causal" = not_causal
  )
}

make_XE_data_split<-function(n, p, corr, betaE, SNR, case, lambda.type){
  
  new_model(name = "gendata_XE_split",
            label = sprintf("n = %s, p = %s, corr = %s, betaE = %s, SNR = %s, case = %s, lambda = %s",
                            n, p, corr, betaE, SNR, case, lambda.type),
            params = list(n = n, p = p, corr = corr, betaE = betaE, SNR = SNR, case = case,
                          lambda.type = lambda.type),
            simulate = function(n, p, corr, betaE, SNR, case, lambda.type, nsim) {
              models <- list()
              for(i in seq(nsim)) {
                
                models[[i]] <- generate_data_case_XE(n, p, corr, betaE, SNR, case)
              }
              return(models)
            })
}

# generate data for X_X case
y_function_XX<-function(n, p, corr=0, SNR = 2, case){
  # covariates
  W <- replicate(n = p, truncnorm::rtruncnorm(n, a = 0, b = 1))
  U <- truncnorm::rtruncnorm(n, a = 0, b = 1)
  V <- truncnorm::rtruncnorm(n, a = 0, b = 1)
  
  Y <- 0
  
  X1 <- (W[, 1] + corr * U) / (1 + corr)
  X2 <- (W[, 2] + corr * U) / (1 + corr)
  X3 <- (W[, 3] + corr * U) / (1 + corr)
  X4 <- (W[, 4] + corr * U) / (1 + corr)
  X5 <- (W[, 5] + corr * U) / (1 + corr)
  
  if (case==1){
    X <- (W[, 6:p] + corr * V) / (1 + corr)
    
    Xall <- cbind(X1, X2, X3, X4, X5, X)
    
    colnames(Xall) <- paste0("X", seq_len(p))
    
    f1 <- function(x) 2 * x
    f2 <- function(x) -1.5 * x
    f3 <- function(x) -2 * x
    f4 <- function(x) 1 * x
    f5 <- function(x) 1 * x
    
    # 5 main effects + 2 interaction
    Y <- f1(X1) + f2(X2) + f3(X3) + f4(X4) + f5(X5) + X1*X2 + X2*X3 + rnorm(n)
  }
  else{
    X6 <- (W[, 6] + corr * U) / (1 + corr)
    X7 <- (W[, 7] + corr * U) / (1 + corr)
    X8 <- (W[, 8] + corr * U) / (1 + corr)
    X9 <- (W[, 9] + corr * U) / (1 + corr)
    X10 <- (W[, 10] + corr * U) / (1 + corr)
    
    X <- (W[, 11:p] + corr * V) / (1 + corr)
    
    Xall <- cbind(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X)
    
    colnames(Xall) <- paste0("X", seq_len(p))
    
    f1 <- function(x) 2 * x
    f2 <- function(x) -1.5 * x
    f3 <- function(x) 1.25 * x
    f4 <- function(x) -1 * x
    f5 <- function(x) -2 * x
    f6 <- function(x) -1 * x
    f7 <- function(x) 1 * x
    f8 <- function(x) 1 * x
    f9 <- function(x) 1 * x
    f10 <- function(x) 1 * x
    # case-2
    # 10 main effects + 3 interaction
    if(case == 2){
      
      Y <- f1(X1) + f2(X2) + f3(X3) + f4(X4) + f5(X5) +
        f6(X6) + f7(X7) + f8(X8) + f9(X9) + f10(X10)+ 
        X1*X2 + X3*X4 + X5*X6 + rnorm(n)
      
    }
    # case-3
    else{
      # 10 main effects + 5 interaction
      Y <- f1(X1) + f2(X2) + f3(X3) + f4(X4) + f5(X5) +
        f6(X6) + f7(X7) + f8(X8) + f9(X9) + f10(X10)+ 
        X1*X2 + X3*X4 + X5*X6 + X1*X3 + X1*X6 + rnorm(n)
    }
  }
  return(list(x = Xall, y = Y))
}

generate_data_case_XX<-function(n, p, corr, SNR, case){
  gendata <- y_function_XX(n, p, corr, SNR, case)
  split_value <- n*0.8
  m <- 2*n
  gendata_valid <- y_function_XX(m, p, corr, SNR, case)
  
  main <- paste0("X", seq_len(p))
  vnames<-c(main)
  
  for (i in (1:(length(main)-1))){
    conc<-paste(main[i], main[(i+1):length(main)], sep = ":")
    vnames<-append(vnames, conc)
  }
  
  true_beta<-matrix(0,nrow=(length(vnames)),ncol=1)
  true_beta_list<-vnames
  rownames(true_beta) <- true_beta_list
  colnames(true_beta) <- "true value"
  
  if (case ==1){
    causal <- c("X1","X2","X3","X4","X5","X1:X2","X2:X3")
    true_beta[1:5,1]<-c(2, -1.5, -2, 1, 1)
    true_beta[p+1,1]<- 1
    true_beta[p+p,1]<- 1
    
  }else if (case ==2){
    causal <- c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X1:X2","X3:X4","X5:X6")
    true_beta[1:10,1]<-c(2, -1.5, 1.25, -1, -2, -1, 1, 1, 1, 1)
    true_beta[p+1,1]<- 1
    true_beta[p+p-1+p-2+1,1]<- 1 #  p+(p-1)+(p-2)+1 = 3*p-2
    true_beta[5*p-9,1]<- 1 # p+(p-1)+(p-2)+(p-3)+(p-4)+1 = 5*p-9
  }else if (case ==3){
    causal <- c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X1:X2","X3:X4","X5:X6","X1:X3","X1:X6")
    true_beta[1:10,1]<-c(2, -1.5, 1.25, -1, -2, -1, 1, 1, 1, 1)
    true_beta[p+1,1]<- 1
    true_beta[p+2,1]<- 1
    true_beta[p+5,1]<- 1
    true_beta[p+p-1+p-2+1,1]<- 1 #  p+(p-1)+(p-2)+1 = 3*p-2
    true_beta[5*p-9,1]<- 1 # p+(p-1)+(p-2)+(p-3)+(p-4)+1 = 5*p-9
  }
  not_causal <- setdiff(vnames, causal)
  
  xtrain <- gendata$x[1:split_value,] 
  ytrain = gendata$y[1:split_value]
  xtest = gendata$x[(split_value+1):n,]
  ytest = gendata$y[(split_value+1):n]
  xvalid = gendata_valid$x
  yvalid = gendata_valid$y
  
  result<-list("x_train" = xtrain , "y_train" = ytrain,
               "x_test" = xtest, "y_test" = ytest,
               "x_valid" = xvalid, "y_valid" = yvalid, 
               "vnames" = vnames, "true_beta" = true_beta,
               "causal" = causal , "not_causal" = not_causal
  )
  return (result)
}

make_XX_data_split<-function(n, p, corr, SNR, case, lambda.type){
  
  new_model(name = "gendata_XX_split",
            label = sprintf("n = %s, p = %s, corr = %s, SNR = %s, case = %s, lambda = %s",
                            n, p, corr, SNR, case, lambda.type),
            params = list(n = n, p = p, corr = corr, SNR = SNR, case = case,
                          lambda.type = lambda.type),
            simulate = function(n, p, corr, SNR, case, lambda.type, nsim) {
              models <- list()
              for(i in seq(nsim)) {
                
                models[[i]] <- generate_data_case_XX(n, p, corr, SNR, case)
              }
              return(models)
            })
}

# generate data for X_X case Binary
y_function_XX_binary<-function(n, p, corr=0, SNR = 2, case){
  # covariates
  W <- replicate(n = p, truncnorm::rtruncnorm(n, a = 0, b = 1))
  U <- truncnorm::rtruncnorm(n, a = 0, b = 1)
  V <- truncnorm::rtruncnorm(n, a = 0, b = 1)
  
  Y <- 0
  
  X1 <- (W[, 1] + corr * U) / (1 + corr)
  X2 <- (W[, 2] + corr * U) / (1 + corr)
  X3 <- (W[, 3] + corr * U) / (1 + corr)
  X4 <- (W[, 4] + corr * U) / (1 + corr)
  X5 <- (W[, 5] + corr * U) / (1 + corr)
  
  if (case==1){
    X <- (W[, 6:p] + corr * V) / (1 + corr)
    
    Xall <- cbind(X1, X2, X3, X4, X5, X)
    
    colnames(Xall) <- paste0("X", seq_len(p))
    
    f1 <- function(x) 2 * x
    f2 <- function(x) -1.5 * x
    f3 <- function(x) -2 * x
    f4 <- function(x) 1 * x
    f5 <- function(x) 1 * x
    
    # 5 main effects + 2 interaction
    Y <- f1(X1) + f2(X2) + f3(X3) + f4(X4) + f5(X5) + X1*X2 + X2*X3 + rnorm(n)
    Y_prob<- 1/(1+exp(-Y))
    Y<- rbinom(n, size = 1, prob = Y_prob)
  }
  else{
    X6 <- (W[, 6] + corr * U) / (1 + corr)
    X7 <- (W[, 7] + corr * U) / (1 + corr)
    X8 <- (W[, 8] + corr * U) / (1 + corr)
    X9 <- (W[, 9] + corr * U) / (1 + corr)
    X10 <- (W[, 10] + corr * U) / (1 + corr)
    
    X <- (W[, 11:p] + corr * V) / (1 + corr)
    
    Xall <- cbind(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X)
    
    colnames(Xall) <- paste0("X", seq_len(p))
    
    f1 <- function(x) 2 * x
    f2 <- function(x) -1.5 * x
    f3 <- function(x) 1.25 * x
    f4 <- function(x) -1 * x
    f5 <- function(x) -2 * x
    f6 <- function(x) -1 * x
    f7 <- function(x) 1 * x
    f8 <- function(x) 1 * x
    f9 <- function(x) 1 * x
    f10 <- function(x) 1 * x
    # case-2
    # 10 main effects + 3 interaction
    if(case == 2){
      
      Y <- f1(X1) + f2(X2) + f3(X3) + f4(X4) + f5(X5) +
        f6(X6) + f7(X7) + f8(X8) + f9(X9) + f10(X10)+ 
        X1*X2 + X3*X4 + X5*X6 + rnorm(n)
      Y_prob<- 1/(1+exp(-Y))
      Y<- rbinom(n, size = 1, prob = Y_prob)
      
    }
    # case-3
    else{
      # 10 main effects + 5 interaction
      Y <- f1(X1) + f2(X2) + f3(X3) + f4(X4) + f5(X5) +
        f6(X6) + f7(X7) + f8(X8) + f9(X9) + f10(X10)+ 
        X1*X2 + X3*X4 + X5*X6 + X1*X3 + X1*X6 + rnorm(n)
      Y_prob<- 1/(1+exp(-Y))
      Y<- rbinom(n, size = 1, prob = Y_prob)
    }
  }
  return(list(x = Xall, y = Y))
}

generate_data_case_XX_binary<-function(n, p, corr, SNR, case){
  gendata <- y_function_XX_binary(n, p, corr, SNR, case)
  split_value <- n*0.8
  m <- 2*n
  gendata_valid <- y_function_XX_binary(m, p, corr, SNR, case)
  
  main <- paste0("X", seq_len(p))
  vnames<-c(main)
  
  for (i in (1:(length(main)-1))){
    conc<-paste(main[i], main[(i+1):length(main)], sep = ":")
    vnames<-append(vnames, conc)
  }
  
  true_beta<-matrix(0,nrow=(length(vnames)),ncol=1)
  true_beta_list<-vnames
  rownames(true_beta) <- true_beta_list
  colnames(true_beta) <- "true value"
  
  if (case ==1){
    causal <- c("X1","X2","X3","X4","X5","X1:X2","X2:X3")
    true_beta[1:5,1]<-c(2, -1.5, -2, 1, 1)
    true_beta[p+1,1]<- 1
    true_beta[p+p,1]<- 1
    
  }else if (case ==2){
    causal <- c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X1:X2","X3:X4","X5:X6")
    true_beta[1:10,1]<-c(2, -1.5, 1.25, -1, -2, -1, 1, 1, 1, 1)
    true_beta[p+1,1]<- 1
    true_beta[p+p-1+p-2+1,1]<- 1 #  p+(p-1)+(p-2)+1 = 3*p-2
    true_beta[5*p-9,1]<- 1 # p+(p-1)+(p-2)+(p-3)+(p-4)+1 = 5*p-9
  }else if (case ==3){
    causal <- c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X1:X2","X3:X4","X5:X6","X1:X3","X1:X6")
    true_beta[1:10,1]<-c(2, -1.5, 1.25, -1, -2, -1, 1, 1, 1, 1)
    true_beta[p+1,1]<- 1
    true_beta[p+2,1]<- 1
    true_beta[p+5,1]<- 1
    true_beta[p+p-1+p-2+1,1]<- 1 #  p+(p-1)+(p-2)+1 = 3*p-2
    true_beta[5*p-9,1]<- 1 # p+(p-1)+(p-2)+(p-3)+(p-4)+1 = 5*p-9
  }
  not_causal <- setdiff(vnames, causal)
  
  xtrain <- gendata$x[1:split_value,] 
  ytrain = gendata$y[1:split_value]
  xtest = gendata$x[(split_value+1):n,]
  ytest = gendata$y[(split_value+1):n]
  xvalid = gendata_valid$x
  yvalid = gendata_valid$y
  
  result<-list("x_train" = xtrain , "y_train" = ytrain,
               "x_test" = xtest, "y_test" = ytest,
               "x_valid" = xvalid, "y_valid" = yvalid, 
               "vnames" = vnames, "true_beta" = true_beta,
               "causal" = causal , "not_causal" = not_causal
  )
  return (result)
}

make_XX_data_split_binary<-function(n, p, corr, SNR, case, lambda.type){
  
  new_model(name = "gendata_XX_split_binary",
            label = sprintf("n = %s, p = %s, corr = %s, SNR = %s, case = %s, lambda = %s",
                            n, p, corr, SNR, case, lambda.type),
            params = list(n = n, p = p, corr = corr, SNR = SNR, case = case,
                          lambda.type = lambda.type),
            simulate = function(n, p, corr, SNR, case, lambda.type, nsim) {
              models <- list()
              for(i in seq(nsim)) {
                
                models[[i]] <- generate_data_case_XX_binary(n, p, corr, SNR, case)
              }
              return(models)
            })
}


