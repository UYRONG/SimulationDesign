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
  }else{
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
  
  result<-list("x_train" = xtrain , "y_train" = ytrain, "e_train" = etrain, "xtrain_lasso" = x_main,
               "x_test" = xtest, "y_test" = ytest,"e_test" = etest, "xtest_lasso" = x_main_test,
               "x_valid" = xvalid, "y_valid" = yvalid, "e_valid" = evalid, "xvalid_lasso" = x_main_valid,
               "vnames" = vnames, "vnames_lasso"=vnames_lasso,"betaE"= betaE,"true_beta" = true_beta,
               "causal" = causal , "not_causal" = not_causal
  )
  return (result)
}
betaE <- 2
# Third parameter indicate the case1
data<-generate_data_case2(100,30,betaE,1)
# Third parameter indicate the case2
# data<-generate_data_case1(150,300,2)
# data<-generate_data_case1(500,2000,3)


####
#          glinternet ---------------------------------------
####
# install.packages("glinternet")
library(glinternet)
fitGL <- glinternet(X = data$xtrain_lasso, Y = data$y_train,
                    numLevels = rep(1, ncol(data$xtrain_lasso)),
                    nLambda = 100,interactionCandidates = c(1),
                    verbose = F)

ytest_hat <- predict(fitGL, X = data$xtest_lasso)
msetest <- colMeans((data$y_test - ytest_hat)^2)
lambda.min.index <- as.numeric(which.min(msetest))
lambda.min <- fitGL$lambda[which.min(msetest)]

yvalid_hat <- predict(fitGL, X = data$xvalid_lasso, lambda = lambda.min)
yvalid_hat
msevalid <- mean((data$y_valid - drop(yvalid_hat))^2)

# nzcoef <- coef(fit, s = lambda.min)[nonzeroCoef(coef(fit, s = lambda.min)),,drop=F]
colnames(data$xtrain_lasso)
tc <- coef(fitGL, lambdaIndex = lambda.min.index)
mains <- colnames(data$xtrain_lasso)[tc[[1]]$mainEffects$cont]

inters <- paste0(colnames(data$xtrain_lasso)[tc[[1]]$interactions$contcont[,2]],":E")
inters
c(mains, inters)
beta_Matrix<-matrix(0,nrow=length(data$vnames_lasso),ncol=1) # variable matrix with name and value
rownames(beta_Matrix) <- data$vnames_lasso

for (i in (1:length(mains))){
  beta_Matrix[mains[i],]<-tc[[1]]$mainEffectsCoef$cont[[i]]
}

for (i in (1:length(inters))){
  beta_Matrix[inters[i],]<-tc[[1]]$interactionsCoef$contcont[[i]]
}

as.matrix(beta_Matrix[beta_Matrix[, 1] != 0, ])

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
