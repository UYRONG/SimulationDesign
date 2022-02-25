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

##
#     backtrack----------
#

# library(LassoBacktracking)
# lassoBT_fit <- LassoBT(data$x_train, data$y_train, iter_max=100)
# # predict must have the sma number of variable as during the training.
# # glinternet:predict(lassoBT_fit, data$x_test)
# ytest_hat<-predict(lassoBT_fit, data$x_test, type = "response")
# dim(ytest_hat)
# msetest <-10000000
# lamb_index <-10000000
# temp_lambda<-1000000
# 
# for (i in 1:dim(ytest_hat)[3]){
#   temp_msetest <- colMeans((data$y_test - ytest_hat[,,i])^2)
#   lambda.min.index <- as.numeric(which.min(temp_msetest))
#   min_lambda<- lassoBT_fit$lambda[which.min(temp_msetest)]
#   print(min_lambda)
#   if (min_lambda<temp_lambda){
#     msetest<- temp_msetest
#     lamb_index<-lambda.min.index
#     temp_lambda<-min_lambda
#   }
# }
# 
# msetest
# lamb_index
# temp_lambda
# 
# yvalid_hat<-predict(lassoBT_fit, data$x_valid, s=temp_lambda)
# yvalid_hat
# msevalid <- mean((data$y_valid - drop(yvalid_hat))^2)
# msevalid
# 
# getNamespaceExports("LassoBacktracking")
# 
# coef<-coef(lassoBT_fit, s=temp_lambda)
# dim(coef)

# --backtrack version
library(LassoBacktracking)
fitBT <- LassoBT(x = data$x_train,
                 y = data$y_train, iter_max=10, nlambda = 100)

ytest_hat <- predict(fitBT, newx = data$x_test)
msetest <- colMeans((data$y_test - ytest_hat)^2)
lambda.min.index <- which(msetest==min(msetest), arr.ind = TRUE)
lambda.min <- fitBT$lambda[lambda.min.index[1,1]]
iter.min <- lambda.min.index[1,2]

coefBT <- as.matrix(predict(fitBT, type = "coef",
                            s = lambda.min, iter = iter.min))
coefBT
nzcoef <- coefBT[sail:::nonzero(coefBT),,drop=F]

ints <- grep(":",rownames(nzcoef)[-1], value=T)
##
# what is the colnames for x_train
##
if (length(ints) != 0) {
  inters <- sapply(stringr::str_split(ints, ":"), function(i) paste(colnames(data$x_train)[as.numeric(as.character(i))], collapse = ":"))
  mains <- colnames(data$x_train)[as.numeric(as.character(setdiff(rownames(nzcoef)[-1], ints)))]
} else {
  inters <- character(0)
  mains <- colnames(data$x_train)[as.numeric(as.character(rownames(nzcoef)[-1]))]
}
inters
active <- if (length(c(mains, inters)) == 0) " " else c(mains, inters)
yvalid_hat <- tryCatch({
  as.matrix(predict(fitBT, newx = data$x_valid, s = lambda.min, iter = iter.min, type = "response"))
},
error = function(err) {
  return(matrix(0, nrow = nrow(data$x_valid), ncol = 1))
} # return NULL on error
)

msevalid <- mean((data$y_valid - drop(yvalid_hat))^2)

returnlist<-list(beta = coefBT,
            # cvfit = fitBT,
            vnames = data$vnames,
            nonzero_coef = nzcoef,
            active = active,
            not_active = setdiff(data$vnames, active),
            yvalid_hat = yvalid_hat,
            msevalid = msevalid,
            causal = data$causal,
            not_causal = data$not_causal,
            yvalid = data$y_valid)
returnlist
