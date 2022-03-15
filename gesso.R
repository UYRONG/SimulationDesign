library(gesso)
library(dplyr)

# data = data.gen(sample_size=100, p=30, 
#                 n_g_non_zero=5, n_gxe_non_zero=3, 
#                 family="gaussian", mode="strong_hierarchical")

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
  vnames<-append(vnames, "E")
  
  for (i in (1:(length(main)))){
    conc<-paste(main[i], "E", sep = ":")
    vnames<-append(vnames, conc)
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
  result<-list("x_train" = gendata$x[1:split_value,] , "y_train" = gendata$y[1:split_value], "e_train" = gendata$e[1:split_value] ,
               "x_test" = gendata$x[(split_value+1):n,], "y_test" = gendata$y[(split_value+1):n],"e_test" = gendata$e[(split_value+1):n],
               "x_valid" = gendata_valid$x, "y_valid" = gendata_valid$y, "e_valid" = gendata_valid$e,
               "vnames" = vnames, "betaE"= betaE,"true_beta" = true_beta,
               "causal" = causal , "not_causal" = not_causal
  )
  return (result)
}
betaE <- 2
# Third parameter indicate the case1
data<-generate_data_case2(100,30,betaE,1)
data$true_beta


fit_model <- gesso.fit(G = data$x_train, E = data$e_train, Y = data$y_train,
                   normalize=TRUE,family = "gaussian", min_working_set_size = 30)


msetest<-1000000
lambda1_index<- -1
lambda2_index<- -1
for(i in 1:(length(fit_model$lambda_1))){
  for(j in 1:(length(fit_model$lambda_2))){
    lambda_pair<-tibble(lambda_1 =fit_model$lambda_1[i],lambda_2 = fit_model$lambda_2[j])
    coeff<-gesso.coef(fit_model,lambda_pair)
    y_test_hat<-gesso.predict(coeff$beta_0,coeff$beta_e,coeff$beta_g,coeff$beta_gxe,data$x_test, data$e_test)
    temp_mse<-mean((y_test_hat-data$y_test)^2)
    
    if (temp_mse < msetest){
      msetest<-temp_mse
      lambda1_index<-i
      lambda2_index<-j
    }
  }
}
msetest
lambda1_index
lambda2_index
min_lambdapair <- tibble(lambda_1 = fit_model$lambda_1[lambda1_index],lambda_2 = fit_model$lambda_2[lambda1_index])
coeff <- gesso.coef(fit_model,min_lambdapair)
y_valid_hat <- gesso.predict(coeff$beta_0,coeff$beta_e,coeff$beta_g,coeff$beta_gxe,data$x_valid, data$e_valid)
mse_valid <- mean((y_valid_hat-data$y_valid)^2)

coeff
data$vnames

active_name_list<-c()
beta_Matrix<-matrix(0,nrow=length(data$vnames)+1,ncol=1)
rowname_list<-append("(Intercept)",data$vnames)
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
    vanme<-paste(xname,":E",sep="")
    active_name_list<-append(active_name_list,vanme)
    beta_Matrix[length(coeff$beta_g)+2+i,]<-coeff$beta_gxe[i]
  }
}


active_name_list
beta_Matrix

beta_Matrix

setdiff(data$vnames, active_name_list)
class(matrix(0,nrow=length(data$vnames),ncol=1))
beta_Matrix<-matrix(0,nrow=length(data$vnames)+1,ncol=1)
rownames(beta_Matrix) <- append(data$vnames,"(Intercept)")
colnames(beta_Matrix) <- 1
beta_Matrix
typeof(nzcoef)
class(nzcoef)

nzout2 <- filter(beta_Matrix, "1" != 0)
nzout2
active_name_list
nactive_name_list<-setdiff(data$vnames,active_name_list)
nactive_name_list
nz<-as.matrix(beta_Matrix[beta_Matrix[, 1] != 0, ])
nz

beta<-as.matrix(beta_Matrix[-1,1])
beta

true_beta <- matrix(rep(0, ncol(design), ncol = 1))
dimnames(true_beta)[[1]] <- rownames(design)
# the first 5 main effects and the first 2 interactions are active
true_beta[c(1:(5*df),(p * df + 2):(p * df + 1 + 2 * df) ),1] <- rnorm(n = 7 * df)
true_beta["X_E",] <- betaE

true_beta<-matrix(0,nrow=20,ncol=1)
true_beta_list<-append("(Intercept)",data$vnames)
rownames(true_beta) <- true_beta_list
colnames(true_beta) <- "true value"
true_beta[2:6,1]<-c(3,-3,-3,3,3)
true_beta[p+2:p+3,]<-c(1.5,-1.5)

