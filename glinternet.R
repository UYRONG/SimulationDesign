library(LassoBacktracking)

####
#          Lasso Backtracking ---------------------------------------
####

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
        x1*x2 + x1*x3 + x1*x4 + x1*x5 + x1*x6 + rnorm(n)
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
  
  result<-list("x_train" = X[1:split_value,] , "y_train" = Y[1:split_value] ,
               "x_test" = X[(split_value+1):n,], "y_test" = Y[(split_value+1):n],
               "x_valid" = X_valid, "y_valid" = Y_valid
)
  return (result)
}
# Third parameter indicate the case1
data<-generate_data_case1(100,30,1)
# Third parameter indicate the case2
# data<-generate_data_case1(150,300,2)
# data<-generate_data_case1(500,2000,3)


lassoBT_fit <- LassoBacktracking::LassoBT(data$x_train, data$y_train, iter_max=100)
# predict must have the sma number of variable as during the training.
# glinternet:predict(lassoBT_fit, data$x_test)
ytest_hat<-predict(lassoBT_fit, data$x_test, type = "response")
dim(ytest_hat)
msetest <-10000000
lamb_index <-10000000
temp_lambda<-1000000

for (i in 1:dim(ytest_hat)[3]){
  temp_msetest <- colMeans((data$y_test - ytest_hat[,,i])^2)
  lambda.min.index <- as.numeric(which.min(temp_msetest))
  min_lambda<- lassoBT_fit$lambda[which.min(temp_msetest)]
  print(min_lambda)
  if (min_lambda<temp_lambda){
    msetest<- temp_msetest
    lamb_index<-lambda.min.index
    temp_lambda<-min_lambda
  }
}

msetest
lamb_index
temp_lambda

yvalid_hat<-predict(lassoBT_fit, data$x_valid, s=temp_lambda)
yvalid_hat
msevalid <- mean((data$y_valid - drop(yvalid_hat))^2)
msevalid

getNamespaceExports("LassoBacktracking")

coef<-coef(lassoBT_fit, s=temp_lambda)
dim(coef)

ls("package:LassoBacktracking")
ls("package:glinternet")
getNamespaceExports("LassoBacktracking")


####
#          hierNet ---------------------------------------
####

install.packages("hierNet")
getNamespaceExports("hierNet")

fit_hierNet<-hierNet(data$x_train, data$y_train, lam=50, strong=TRUE)
print(fit_hierNet)
fit_hierNet_path<-hierNet.path(data$x_train, data$y_train, lam=50, strong=TRUE)
print(fit_hierNet_path)
ytest_hat <- predict(fit_hierNet, newx = data$x_test)
dim(ytest_hat)
ytest_hat_col<-ytest_hat[,1]

msetest <- mean((data$y_test - ytest_hat_col)^2)
msetest
# can not define the lambda
# still do valid?


####
#          glinternet ---------------------------------------
####

# Number of levels for each variable, of length nvars. Set to 1 for continuous variables.
numLevels <- rep(1, 30)
fit_glin <- glinternet(data$x_train, data$y_train, numLevels)
ytest_hat<-predict(fit_glin,data$x_test,type="response" )


class(fit_glin)


msetest <- colMeans((data$y_test - ytest_hat)^2)
msetest
lambda.min.index <- as.numeric(which.min(msetest))
lambda.min <- fit_glin$lambda[which.min(msetest)]
lambda.min

yvalid_hat <- predict(fit_glin, data$x_valid, lambda = lambda.min)
dim(yvalid_hat)
msevalid <- mean((data$y_valid - drop(yvalid_hat))^2)
msevalid

coeffs <- coef(fit_glin,lambdaIndex = lambda.min.index)
coeffs

####
#          xyz ---------------------------------------
####

library(devtools)
install_github("gathanei/xyz")
library(xyz)

# #construct a data matrix X
# X<-matrix(sample(c(-1,1),replace=T,n*p),n,p)

# xyz_regression(xyz_regression(X, Y, lambdas = NULL, n_lambda = 10, alpha = 0.9, L = 10,
# standardize = TRUE, standardize_response = TRUE))
fit_xyz<-xyz_regression(data$x_train,data$y_train,n_lambda=10,alpha=0.9,L=10)
fit_xyz

# without specified predict func



# xyz_search(X, Y, L = 10, N = 100, binary = TRUE, negative = TRUE)
# L: An integer indicating how many projection steps are performed.
# N: A integer, controlling the number of pairs that will be returned in the end.
# binary: A logical indicating if X is binary or continuous
# negative A logical indicating if also negative interactions should be searched for.

result<-xyz_search(X,Y,L=10,N=10,binary=T,negative=T)

xyz

# gaussian response, continuous features
Y = rnorm(100)
X = matrix(rnorm(100*10), nrow=100)
numLevels = rep(1, 10)
fit = glinternet(X, Y, numLevels)

#binary response, continuous features
Y = rbinom(100, 1, 0.5)
fit = glinternet(X, Y, numLevels, family="binomial")

#binary response, categorical variables
X = matrix(sample(0:2, 100*10, replace=TRUE), nrow=100)
numLevels = rep(3, 10)
fit = glinternet(X, Y, numLevels, family="binomial")

Y = rnorm(100)
numLevels = sample(1:5, 10, replace=TRUE)
X = sapply(numLevels, function(x) if (x==1)
  rnorm(100) else sample(0:(x-1), 100, replace=TRUE))
fit = glinternet(X, Y, numLevels)
max(abs(fit$fitted - predict(fit, X)))
