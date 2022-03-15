install.packages("sail")
library(sail)


y_function<-function(n, p, corr=0, betaE = 2, SNR = 2, case){
  # covariates
  W <- replicate(n = p, truncnorm::rtruncnorm(n, a = 0, b = 1))
  U <- truncnorm::rtruncnorm(n, a = 0, b = 1)
  V <- truncnorm::rtruncnorm(n, a = 0, b = 1)
  
  E <- stats::rnorm(n)
  
  Y<-0
  
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
generate_data_case2<-function(n,p,betaE,case){
  # test and training set
  # X <- matrix(rnorm(n*p), nrow=n)
  # X <- scale(X,TRUE,TRUE)
  gendata <- y_function(n,p,0,0,betaE,case)
  
  # Y <- y_function(n,p,0,0,2,case)
  split_value <- n*0.8
  m <- 2*n
  gendata_valid <- y_function(m,p,0,0,betaE,case)
  # X_valid <- matrix(rnorm(m*p), nrow=m)
  # Y_valid <-y_function(n,p,0,0,2,case)
  
  main <- paste0("X", seq_len(p))
  vnames<-c(main)
  
  for (i in (1:(length(main)))){
    conc<-paste(main[i], "E", sep = ":")
    vnames<-append(vnames, conc)
  }
  if (case ==1){
    causal <- c("X1","X2","X3","X4","X5","X4:E","X5:E")
  }else if (case ==2){
    causal <- c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X4:E","X5:E","X6:E","X7:E","X8:E")
  }else if (case ==3){
    causal <- c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X4:E","X5:E","X6:E","X7:E","X8:E")
  }
  not_causal <- setdiff(vnames, causal)
  result<-list("x_train" = gendata$x[1:split_value,] , "y_train" = gendata$y[1:split_value], "e_train" = gendata$e[1:split_value] ,
               "x_test" = gendata$x[(split_value+1):n,], "y_test" = gendata$y[(split_value+1):n],"e_test" = gendata$e[(split_value+1):n],
               "x_valid" = gendata_valid$x, "y_valid" = gendata_valid$y, "e_valid" = gendata_valid$e,
               "vnames" = vnames, "betaE"= betaE,
               "causal" = causal , "not_causal" = not_causal
  )
  return (result)
}
betaE <- 2
# Third parameter indicate the case1
data<-generate_data_case2(100,30,betaE,1)

fit <- sail(x = data$x_train, y = data$y_train, e = data$e_train,
            basis = function(i) i)

ytest_hat <- predict(fit, newx = data$x_test, newe = data$e_test)
ytest_hat
msetest <- colMeans((data$y_test - ytest_hat)^2)
lambda.min.index <- as.numeric(which.min(msetest))
lambda.min <- fit$lambda[which.min(msetest)]

yvalid_hat <- predict(fit, newx = data$x_valid, newe = data$e_valid, s = lambda.min)
yvalid_hat
msevalid <- mean((data$y_valid - drop(yvalid_hat))^2)

nzcoef <- predict(fit, s = lambda.min, type = "nonzero")

as.matrix(coef(fit, s = lambda.min))[-1,,drop=F]

beta = coef(fit, s = lambda.min)
dim(data$x_train)

nzcoef
fit$active
fit$active[[lambda.min.index]]

fit$active[[lambda.min.index]]

coef(fit, s = lambda.min)[-1,,drop=F]
nzcoef
active
# return(list(beta = coef(fit, s = lambda.min)[-1,,drop=F],
#             # fit = fit,
#             vnames = draw[["vnames"]],
#             nonzero_coef = nzcoef,
#             active = fit$active[[lambda.min.index]],
#             not_active = setdiff(draw[["vnames"]], fit$active[[lambda.min.index]]),
#             yvalid_hat = yvalid_hat,
#             msevalid = msevalid,
#             causal = draw[["causal"]],
#             not_causal = draw[["not_causal"]],
#             yvalid = draw[["yvalid"]]))
# },