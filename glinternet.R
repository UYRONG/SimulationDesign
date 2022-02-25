y_function<-function(x,n,case){
  y<-0
  if (case==1){
    f1 <- function(x) 2 * x
    f2 <- function(x) -1.5 * x
    f3 <- function(x) -2 * x
    f4 <- function(x) 1 * x
    f5 <- function(x) 1 * x
    x1 <-x[, 1]
    x2 <-x[, 2]
    x3 <-x[, 3]
    x4 <-x[, 4]
    x5 <-x[, 5]
    # 5 main effects + 2 interaction
    y <- f1(x1) + f2(x2) + f3(x3) + f4(x4) + f5(x5) + x1*x2 + x2*x3 + rnorm(n)
  }
  else{
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
    x1 <-x[, 1]
    x2 <-x[, 2]
    x3 <-x[, 3]
    x4 <-x[, 4]
    x5 <-x[, 5]
    x6 <-x[, 6]
    x7 <-x[, 7]
    x8 <-x[, 8]
    x9 <-x[, 9]
    x10 <-x[, 10]
    # case-2
    # 10 main effects + 3 interaction
    if(case == 2){
      y <- f1(x1) + f2(x2) + f3(x3) + f4(x4) + f5(x5) +
        f6(x6) + f7(x7) + f8(x8) + f9(x9) + f10(x10)+ 
        x1*x2 + x3*x4 + x5*x6 +rnorm(n)
    }
    # case-3
    else{
      # 10 main effects + 5 interaction
      y <- f1(x1) + f2(x2) + f3(x3) + f4(x4) + f5(x5) +
        f6(x6) + f7(x7) + f8(x8) + f9(x9) + f10(x10)+ 
        x1*x2 + x3*x4 + x5*x6 + x1*x3 + x1*x6 + rnorm(n)
    }
  }
  return(y)
}
generate_data_case1<-function(n,p,case){
  # test and training set
  X <- matrix(rnorm(n*p), nrow=n)
  X <- scale(X,TRUE,TRUE)
  Y <-y_function(X,n,case)
  split_value <- n*0.8
  m <- 2*n
  X_valid <- matrix(rnorm(m*p), nrow=m)
  Y_valid <-y_function(X,m,case)
  
  main <- paste0("X", seq_len(p))
  vnames<-c(main)
  
  for (i in (1:(length(main)-1))){
    conc<-paste(main[i], main[(i+1):length(main)], sep = ":")
    vnames<-append(vnames, conc)
  }
  if (case ==1){
    causal <- c("X1","X2","X3","X4","X5","X1:X2","X2:X3")
  }else if (case ==2){
    causal <- c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X1:X2","X3:X4","X5:X6")
  }else if (case ==3){
    causal <- c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X1:X2","X3:X4","X5:X6","X1:X3","X1:X6")
  }
  not_causal <- setdiff(vnames, causal)
  result<-list("x_train" = X[1:split_value,] , "y_train" = Y[1:split_value] ,
               "x_test" = X[(split_value+1):n,], "y_test" = Y[(split_value+1):n],
               "x_valid" = X_valid, "y_valid" = Y_valid,
               "vnames" = vnames,
               "causal" = causal , "not_causal" = not_causal
              )
  return (result)
}
# Third parameter indicate the case1
data<-generate_data_case1(100,30,1)
# Third parameter indicate the case2
# data<-generate_data_case1(150,300,2)
# data<-generate_data_case1(500,2000,3)


####
#          glinternet ---------------------------------------
####
install.packages("glinternet")
library(glinternet)
fitGL <- glinternet(X = data$x_train, Y = data$y_train,
                    numLevels = rep(1, ncol(data$x_train)),
                    nLambda = 100,
                    verbose = F)

ytest_hat <- predict(fitGL, X = data$x_test)
msetest <- colMeans((data$y_test - ytest_hat)^2)
lambda.min.index <- as.numeric(which.min(msetest))
lambda.min <- fitGL$lambda[which.min(msetest)]

yvalid_hat <- predict(fitGL, X = data$x_valid, lambda = lambda.min)
msevalid <- mean((data$y_valid - drop(yvalid_hat))^2)

# nzcoef <- coef(fit, s = lambda.min)[nonzeroCoef(coef(fit, s = lambda.min)),,drop=F]

tc <- coef(fitGL, lambdaIndex = lambda.min.index)
mains <- colnames(draw[["xtrain_lasso"]])[tc[[1]]$mainEffects$cont]
inters <- paste0(colnames(draw[["xtrain_lasso"]])[tc[[1]]$interactions$contcont[,2]],":E")

return_list<-list(beta = NULL,
            vnames = draw[["vnames_lasso"]],
            nonzero_coef = NULL,
            active = c(mains, inters),
            not_active = setdiff(draw[["vnames"]], c(mains, inters)),
            yvalid_hat = yvalid_hat,
            msevalid = msevalid,
            causal = draw[["causal"]],
            not_causal = draw[["not_causal"]],
            yvalid = draw[["yvalid"]])








# Number of levels for each variable, of length nvars. Set to 1 for continuous variables.
# numLevels <- rep(1, 30)
# fit_glin <- glinternet(data$x_train, data$y_train, numLevels, nLambda = 100,
#                        verbose = F)
# ytest_hat<-predict(fit_glin,data$x_test,type="response" )
# 
# 
# class(fit_glin)
# 
# 
# msetest <- colMeans((data$y_test - ytest_hat)^2)
# msetest
# lambda.min.index <- as.numeric(which.min(msetest))
# lambda.min <- fit_glin$lambda[which.min(msetest)]
# lambda.min
# 
# yvalid_hat <- predict(fit_glin, data$x_valid, lambda = lambda.min)
# dim(yvalid_hat)
# msevalid <- mean((data$y_valid - drop(yvalid_hat))^2)
# msevalid
# 
# coeffs <- coef(fit_glin,lambdaIndex = lambda.min.index)
# coeffs

####
#          xyz ---------------------------------------
####

# library(devtools)
# install_github("gathanei/xyz")
# library(xyz)
# 
# # #construct a data matrix X
# # X<-matrix(sample(c(-1,1),replace=T,n*p),n,p)
# 
# # xyz_regression(xyz_regression(X, Y, lambdas = NULL, n_lambda = 10, alpha = 0.9, L = 10,
# # standardize = TRUE, standardize_response = TRUE))
# fit_xyz<-xyz_regression(data$x_train,data$y_train,n_lambda=10,alpha=0.9,L=10)
# fit_xyz
# 
# # without specified predict func
# 
# 
# 
# # # xyz_search(X, Y, L = 10, N = 100, binary = TRUE, negative = TRUE)
# # # L: An integer indicating how many projection steps are performed.
# # # N: A integer, controlling the number of pairs that will be returned in the end.
# # # binary: A logical indicating if X is binary or continuous
# # # negative A logical indicating if also negative interactions should be searched for.
# # 
# # result<-xyz_search(X,Y,L=10,N=10,binary=T,negative=T)
# # 
# # xyz
# # 
# # # gaussian response, continuous features
# # Y = rnorm(100)
# # X = matrix(rnorm(100*10), nrow=100)
# # numLevels = rep(1, 10)
# # fit = glinternet(X, Y, numLevels)
# # 
# # #binary response, continuous features
# # Y = rbinom(100, 1, 0.5)
# # fit = glinternet(X, Y, numLevels, family="binomial")
# # 
# # #binary response, categorical variables
# # X = matrix(sample(0:2, 100*10, replace=TRUE), nrow=100)
# # numLevels = rep(3, 10)
# # fit = glinternet(X, Y, numLevels, family="binomial")
# # 
# # Y = rnorm(100)
# # numLevels = sample(1:5, 10, replace=TRUE)
# # X = sapply(numLevels, function(x) if (x==1)
# #   rnorm(100) else sample(0:(x-1), 100, replace=TRUE))
# # fit = glinternet(X, Y, numLevels)
# # max(abs(fit$fitted - predict(fit, X)))
