# install.packages("devtools",dependencies = TRUE)
#library(devtools)
#devtools::install_github("NataliaZemlianskaia/gesso")

#install.packages("simulator")

library(gesso)
library(dplyr)
source("/Users/uyrong/Desktop/Simulation/SimulationDesign/sims/gen_data.R")
source("/Users/uyrong/Desktop/Simulation/SimulationDesign/sims/model_functions.R")
#---------------------------------try----------------------------------------------------------------------------------
data = data.gen()
tune_model = gesso.cv(data$G_train, data$E_train, data$Y_train)
coefficients = gesso.coef(tune_model$fit, tune_model$lambda_min)
beta_0 = coefficients$beta_0; beta_e = coefficients$beta_e
beta_g = coefficients$beta_g; beta_gxe = coefficients$beta_gxe
new_G = data$G_test; new_E = data$E_test
new_Y = gesso.predict(beta_0, beta_e, beta_g, beta_gxe, new_G, new_E)
cor(new_Y, data$Y_test)^2
length(beta_0)

data = data.gen(sample_size=100, p=30, 
                n_g_non_zero=5, n_gxe_non_zero=3, 
                family="gaussian", mode="strong_hierarchical")

fit_model <- gesso.fit(G = data$G_train, E = data$E_train, Y = data$Y_train,
                   normalize=TRUE,family = "gaussian", min_working_set_size = 30)

for (value in fit_model$beta_0){
ytest_hat <- gesso.predict(fit_model$beta_0, fit_model$beta_e, fit_model$beta_g,
                           fit_model$beta_gxe, data$G_test, data$E_test)

}
length(fit_model$beta_0) 
length(data$Y_test)
msetest <- (data$Y_test - ytest_hat)^2

msetest

min_mse_index <- as.numeric(which.min(msetest))

min_mse_index

#problem fit model has two lambda
fit_model$lambda_1[which.min(msetest)]
fit_model$lambda_2[which.min(msetest)]

min(fit_model$lambda_1)
min(fit_model$lambda_2)

min_lambda <-tibble(lambda_1=min(model$lambda_1), lambda_2=min(model$lambda_2))

new_model <- gesso.coef(fit=model, lambda=min_lambda)

Y_valid_hat <- gesso.predict(new_model$beta_0,new_model$beta_e, new_model$beta_g,
                        new_model$beta_gxe, data$G_valid, data$E_valid)
Y_valid_hat

msevalid <- mean((data$Y_valid - drop(Y_valid_hat))^2)

msevalid

# question about calculating in the part!!!

length(fit_model$lambda_1)

min_value <- 100
min_lambda1 <-0
min_lambda2 <-0

model <- gesso.fit(G = data$G_train, E = data$E_train, Y = data$Y_train,
                   normalize=TRUE,family = "gaussian", min_working_set_size = 30)

for (value in fit_model$lambda_1){
  print("check")
  for (value2 in fit_model$lambda_2){
    min_l <-tibble(lambda_1=value, lambda_2=value2)
    new_model <- gesso.coef(fit=model, lambda=min_l)
    Y_valid_hat <- gesso.predict(new_model$beta_0,new_model$beta_e, new_model$beta_g,
                                 new_model$beta_gxe, data$G_valid, data$E_valid)
    msevalid <- mean((data$Y_valid - drop(Y_valid_hat))^2)
    if (msevalid < min_value){
      min_value <- msevalid
      min_lambda1 <- value
      min_lambda2 <- value2
    }
  }
}
min_value #  0.438307
min_lambda1 # 0.0975839
min_lambda2 # 0.01788694



# is nzcoef a number or a list?
#parameter for new_model is beta_0/beta_e/beta_gxe/beta_c

new_model

#---------------------------------try----------------------------------------------------------------------------------


gendataPaper <- function(n, p, corr = 0,
                         E = truncnorm::rtruncnorm(n, a = -1, b = 1),
                         # E = rbinom(n,1,0.5),
                         betaE = 2, SNR = 2, hierarchy = c("strong", "weak", "none"),
                         nonlinear = TRUE, interactions = TRUE, causal, not_causal) {
  # this is modified from "VARIABLE SELECTION IN NONPARAMETRIC ADDITIVE MODEL" huang et al, Ann Stat.
  # n = 200
  # p = 10
  # corr = 1
  
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
  
  # W <- replicate(n = p, rnorm(n))
  # U <- rnorm(n)
  # V <- rnorm(n)
  
  X1 <- (W[, 1] + corr * U) / (1 + corr)
  X2 <- (W[, 2] + corr * U) / (1 + corr)
  X3 <- (W[, 3] + corr * U) / (1 + corr)
  X4 <- (W[, 4] + corr * U) / (1 + corr)
  
  X <- (W[, 5:p] + corr * V) / (1 + corr)
  
  Xall <- cbind(X1, X2, X3, X4, X)
  
  colnames(Xall) <- paste0("X", seq_len(p))
  
  # see "Variable Selection in NonParametric Addditive Model" Huang Horowitz and Wei
  if (nonlinear) {
    
    f1 <- function(x) 5 * x
    f2 <- function(x) 3 * (2 * x - 1)^2
    f3 <- function(x) 4 * sin(2 * pi * x) / (2 - sin(2 * pi * x))
    f4 <- function(x) 6 * (0.1 * sin(2 * pi * x) + 0.2 * cos(2 * pi * x) +
                             0.3 * sin(2 * pi * x)^2 + 0.4 * cos(2 * pi * x)^3 +
                             0.5 * sin(2 * pi * x)^3)
    f3.inter = function(x, e) e * f3(x)
    f4.inter = function(x, e) e * f4(x)
    
  } else {
    f1 <- function(x)  5 * x
    f2 <- function(x)  3 * (x + 1)
    f3 <- function(x)  4 * x
    f4 <- function(x)  6 * (x - 2)
    f3.inter <- function(x, e) e * f3(x)
    f4.inter <- function(x, e) e * f4(x)
    
  }
  # error
  error <- stats::rnorm(n)
  
  if (!nonlinear) {
    
    Y.star <- f1(X1) +
      f2(X2) +
      f3(X3) +
      f4(X4) +
      betaE * E +
      f3.inter(X3,E) +
      f4.inter(X4,E)
    
    scenario <- "2"
    
  } else {
    if (!interactions) {
      # main effects only; non-linear Scenario 3
      Y.star <- f1(X1) +
        f2(X2) +
        f3(X3) +
        f4(X4) +
        betaE * E
      scenario <- "3"
    } else {
      if (hierarchy == "none" & interactions) {
        # interactions only; non-linear
        Y.star <- E * f3(X3) +
          E * f4(X4)
        scenario <- "1c"
      } else if (hierarchy == "strong" & interactions) {
        # strong hierarchy; non-linear
        Y.star <- f1(X1) +
          f2(X2) +
          f3(X3) +
          f4(X4) +
          betaE * E +
          E * f3(X3) +
          E * f4(X4)
        scenario <- "1a"
      } else if (hierarchy == "weak" & interactions) {
        # weak hierarchy; linear
        Y.star <- f1(X1) +
          f2(X2) +
          # f3(X3) +
          # f4(X4) +
          betaE * E +
          E * f3(X3) +
          E * f4(X4)
        scenario <- "1b"
      }
    }
  }
  
  k <- sqrt(stats::var(Y.star) / (SNR * stats::var(error)))
  
  Y <- Y.star + as.vector(k) * error
  
  return(list(
    x = Xall, y = Y, e = E, Y.star = Y.star, f1 = f1(X1),
    f2 = f2(X2), f3 = f3(X3), f4 = f4(X4), betaE = betaE,
    f1.f = f1, f2.f = f2, f3.f = f3, f4.f = f4,
    X1 = X1, X2 = X2, X3 = X3, X4 = X4, scenario = scenario,
    causal = causal, not_causal = not_causal
  ))
}


gendata <- function(n, p, corr, E = truncnorm::rtruncnorm(n, a = -1, b = 1),
                    betaE, SNR, parameterIndex) {
  if (!requireNamespace("truncnorm", quietly = TRUE)) {
    stop("Package \"truncnorm\" needed for this function to simulate data. Please install it.",
         call. = FALSE
    )
  }
  
  main <- paste0("X", seq_len(p))
  vnames <- c(main, "E", paste0(main, ":E"))
  
  if (parameterIndex == 1) { # 1a
    hierarchy <- "strong"
    nonlinear <- TRUE
    interactions <- TRUE
    causal <- c("X1", "X2", "X3", "X4", "E", "X3:E", "X4:E")
  } else if (parameterIndex == 2) { # 1b
    hierarchy <- "weak"
    nonlinear <- TRUE
    interactions <- TRUE
    causal <- c("X1", "X2", "E", "X3:E", "X4:E")
  } else if (parameterIndex == 3) { # 1c
    hierarchy <- "none"
    nonlinear <- TRUE
    interactions <- TRUE
    causal <- c("X3:E", "X4:E")
  } else if (parameterIndex %in% c(4,6)) { # 2
    hierarchy <- "strong"
    nonlinear <- FALSE
    interactions <- TRUE
    causal <- c("X1", "X2", "X3", "X4", "E", "X3:E", "X4:E")
  } else if (parameterIndex == 5) { # 3
    hierarchy <- "strong"
    nonlinear <- TRUE
    interactions <- FALSE
    causal <- c("X1", "X2", "X3", "X4", "E")
  }
  
  not_causal <- setdiff(vnames, causal)
  
  DT <- gendataPaper(
    n = n, p = p, corr = corr,
    E = E,
    betaE = betaE, SNR = SNR,
    hierarchy = hierarchy, nonlinear = nonlinear, interactions = interactions,
    causal = causal, not_causal = not_causal
  )
  return(DT)
}
try <-gendata (100, 30, 0, 2, 2, lambda.type = "lambda.min", parameterIndex=1)

## tune the model hyperparameters   
# tune_model = gesso.cv(data$G_train, data$E_train, data$Y_train,
#                       grid_size=20, tolerance=1e-4,
#                       parallel=TRUE, nfold=4,
#                       normalize=TRUE, normalize_response=TRUE,
#                       seed=1)
# tune_model$lambda_min
# new_model = gesso.coef(fit=tune_model$fit, lambda=tune_model$lambda_min)

model <- gesso.fit(G = data$G_train, E = data$E_train, Y = data$Y_train,
          normalize=TRUE,family = "gaussian", min_working_set_size = 30)
model$beta_gxe
ytest_hat <- predict(model, newx = data$G_test, newe = data$E_test)
ytest_hat

msetest <- colMeans((data$Y_test - ytest_hat)^2)
msetest

coef_model = gesso.coef(fit=model, lambda=min_lambda)

new_Y = gesso.predict(model$beta_0, model$beta_e, model$beta_g, model$beta_gxe, data$G_test, data$E_test)
new_Y

colMeans((data$Y_test - new_Y)^2)
length(data$Y_test)
li = (data$Y_test - new_Y)^2
li

lambda.min.index <- as.numeric(which.min(li))
lambda.min.index

lambda.min1 <- model$lambda_2[which.min(li)]
lambda.min1

lambda.min2 <- model$lambda_2[which.min(li)]
lambda.min2

min(model$lambda_1)


tune_model$lambda_min

min_lambda <-tibble(lambda_1=min(model$lambda_1), lambda_2=min(model$lambda_2))
min_lambda
newnew = gesso.coef(fit=model, lambda=min_lambda)
newnew
new_Y1 = gesso.predict(new_model$beta_0,new_model$beta_e, new_model$beta_g, new_model$beta_gxe, data$G_test, data$E_test)
new_Y1





coefficients = gesso.coef(fit=model$fit, lambda=model$lambda_min)
model$lambda_min
coefficients
# predicted_value <-predict
beta_0 = coefficients$beta_0
beta_e = coefficients$beta_e
beta_g = coefficients$beta_g
beta_gxe = coefficients$beta_gxe
new_G = draw[["xtest"]]
new_E = draw[["etest"]]
ytest_hat = gesso.predict(beta_0, beta_e, beta_g, beta_gxe, new_G, new_E)
ytest_hat

## obtain interaction and main effect coefficietns corresponding to the best model
coefficients = gesso.coef(fit=tune_model$fit, lambda=tune_model$lambda_min)
gxe_coefficients = coefficients$beta_gxe                      
g_coefficients = coefficients$beta_g    

## check if all non-zero features were recovered by the model
features<-cbind(data$Beta_GxE[data$Beta_GxE != 0], gxe_coefficients[data$Beta_GxE != 0])
print(features)

##      [,1]       [,2]
## [1,] -1.5 -0.9711450
## [2,] -1.5 -1.9493914
## [3,] -1.5 -0.5148486
## [4,] -1.5 -1.4827039
## [5,]  1.5  1.8539925

## check if the largest estimated interaction effects correspond to the true non-zero coefficients
estimated_interaction<-(data$Beta_GxE[order(abs(gxe_coefficients), decreasing=TRUE)])[1:10]
## [1] -1.5  1.5 -1.5  1.5 -1.5  0.0  0.0  0.0  0.0  0.0
print(estimated_interaction)
## calculate principal selection metrics
selection = selection.metrics(true_b_g=data$Beta_G, 
                              true_b_gxe=data$Beta_GxE, 
                              estimated_b_g=g_coefficients,
                              estimated_b_gxe=gxe_coefficients)
summary<-cbind(selection)
print(summary)

##                 selection 
## b_g_non_zero    77        
## b_gxe_non_zero  63        
## mse_b_g         0.2775684 
## mse_b_gxe       0.5617425 
## sensitivity_g   1         
## specificity_g   0.8282051 
## precision_g     0.1298701 
## sensitivity_gxe 1         
## specificity_gxe 0.8531646 
## precision_gxe   0.07936508
## auc_g           0.999999  
## auc_gxe         0.999998  


