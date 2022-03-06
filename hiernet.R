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
#          hierNet ---------------------------------------
####

install.packages("hierNet")
library(hierNet)
getNamespaceExports("hierNet")

# fit_hierNet<-hierNet(data$x_train, data$y_train, lam=50, strong=TRUE)

# The reason for using hierNet.path instead of hierNet
# hierNet_path can provide a list of lambda value, whie the hierNet you have to specify the value for lambda before using it
fit_hierNet_path<-hierNet.path(data$x_train, data$y_train,strong=TRUE)
# print(fit_hierNet_path)
# fit_hierNet_path$lamlist
ytest_hat <- predict(fit_hierNet_path, newx = data$x_test)
# ytest_hat
# dim(ytest_hat[1])
# length(ytest_hat[1,])
col_msetest<-c()
for (i in 1:length(ytest_hat[1,])){
  msetest <- mean((data$y_test - ytest_hat[,i])^2)
  col_msetest<-append(col_msetest,msetest)
}
lambda_min_index<-which.min(col_msetest)
lambda_min<-fit_hierNet_path$lamlist[which.min(col_msetest)]

yvalid_hat <- predict(fit_hierNet_path, newx = data$x_valid)[,lambda_min_index]
msevalid <- mean((data$y_valid - yvalid_hat)^2)
length(yvalid_hat)
msevalid

# bp:p-vector of estimated "positive part" main effect (p=# features)
# bn:p-vector of estimated "positive part" main effect (p=# features)
# bp-bn:overall main effect estimated
beta_main<-fit_hierNet_path$bp[,lambda_min_index]-fit_hierNet_path$bn[,lambda_min_index]#overall main effect estimated
as.vector(beta_main)
beta_main
fit_hierNet_path$th[,,lambda_min_index]

data$vnames








