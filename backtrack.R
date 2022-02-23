
install.packages("LassoBacktracking")
library(LassoBacktracking)
x <- matrix(rnorm(100*250), 100, 250)
y <- x[, 1] + x[, 2] - x[, 1]*x[, 2] + x[, 3] + rnorm(100)
dim(x)

length(y)
out <- LassoBT(x, y, iter_max=10)
predict(out, newx=x[1:2, ])

x

# gendataPaper <- function(n, p, corr = 0,
#                          # E = truncnorm::rtruncnorm(n, a = -1, b = 1),
#                          # E = rbinom(n,1,0.5),
#                          # betaE = 2, SNR = 2, hierarchy = c("strong", "weak", "none"),
#                          # nonlinear = TRUE, 
#                          interactions = TRUE, causal, not_causal) {
#   # this is modified from "VARIABLE SELECTION IN NONPARAMETRIC ADDITIVE MODEL" huang et al, Ann Stat.
#   # n = 200
#   # p = 10
#   # corr = 1
#   
#   if (!requireNamespace("truncnorm", quietly = TRUE)) {
#     stop("Package \"truncnorm\" needed for this function to simulate data. Please install it.",
#          call. = FALSE
#     )
#   }
#   
#   hierarchy <- match.arg(hierarchy)
#   
#   # covariates
#   W <- replicate(n = p, truncnorm::rtruncnorm(n, a = 0, b = 1))
#   U <- truncnorm::rtruncnorm(n, a = 0, b = 1)
#   V <- truncnorm::rtruncnorm(n, a = 0, b = 1)
#   
#   # W <- replicate(n = p, rnorm(n))
#   # U <- rnorm(n)
#   # V <- rnorm(n)
#   
#   X1 <- (W[, 1] + corr * U) / (1 + corr)
#   X2 <- (W[, 2] + corr * U) / (1 + corr)
#   X3 <- (W[, 3] + corr * U) / (1 + corr)
#   X4 <- (W[, 4] + corr * U) / (1 + corr)
#   
#   X <- (W[, 5:p] + corr * V) / (1 + corr)
#   
#   Xall <- cbind(X1, X2, X3, X4, X)
#   
#   colnames(Xall) <- paste0("X", seq_len(p))
#   
#   # see "Variable Selection in NonParametric Addditive Model" Huang Horowitz and Wei
#  
#     
#     f1 <- function(x) 5 * x
#     f2 <- function(x)  3 * (x + 1)
#     f3 <- function(x)  4 * x
#     f4 <- function(x)  6 * (x - 2)
#     f3.inter <- function(x, e) e * f3(x)
#     f4.inter <- function(x, e) e * f4(x)
#     
#     Y.star <- f1(X1) +
#       f2(X2) +
#       f3(X3) +
#       f4(X4) +
# 
#   # error
#   error <- stats::rnorm(n)
#   
#   k <- sqrt(stats::var(Y.star) / (SNR * stats::var(error)))
#   
#   Y <- Y.star + as.vector(k) * error
#   
#   return(list(
#     x = Xall, y = Y, Y.star = Y.star, f1 = f1(X1),
#     f2 = f2(X2), f3 = f3(X3), f4 = f4(X4),
#     f1.f = f1, f2.f = f2, f3.f = f3, f4.f = f4,
#     X1 = X1, X2 = X2, X3 = X3, X4 = X4, 
#     causal = causal, not_causal = not_causal
#   ))
# }
# 
# gendata <- function(n, p, corr, SNR) {
#   if (!requireNamespace("truncnorm", quietly = TRUE)) {
#     stop("Package \"truncnorm\" needed for this function to simulate data. Please install it.",
#          call. = FALSE
#     )
#   }
#   
#   main <- paste0("X", seq_len(p))
#   vnames <- c(main, paste0(main))
#   
#   if (parameterIndex == 1) { # 1a
#     interactions <- TRUE
#     causal <- c("X1", "X2", "X3", "X4")
#   } else if (parameterIndex == 2) { # 1b
#     interactions <- TRUE
#     causal <- c("X1", "X2", "X3", "X4")
#   } else if (parameterIndex == 3) { # 1c
#     interactions <- TRUE
#     causal <- c("X1", "X2", "X3", "X4")
#   } 
#   
#   not_causal <- setdiff(vnames, causal)
#   
#   DT <- gendataPaper(
#     n = n, p = p, corr = corr, SNR = SNR,
#     interactions = interactions,
#     causal = causal, not_causal = not_causal
#   )
#   return(DT)
# }
# 
# make_gendata_Paper_data_split <- function(n, p, corr, SNR, lambda.type, parameterIndex) {
#   
#   main <- paste0("X", seq_len(p))
#   vnames <- c(main, paste0(main))
#   # used for glmnet and lasso backtracking
#   # f.identity <- function(i) i
#   
#   new_model(name = "gendata_thesis_split_v4",
#             label = sprintf("n = %s, p = %s, corr = %s, betaE = %s, SNR = %s, index = %s, lambda = %s",
#                             n, p, corr, SNR, parameterIndex, lambda.type),
#             params = list(n = n, p = p, corr = corr, SNR = SNR,
#                           lambda.type = lambda.type, parameterIndex = parameterIndex),
#             simulate = function(n, p, corr, SNR, parameterIndex, lambda.type, nsim) {
#               models <- list()
#               for(i in seq(nsim)) {
#                 # test and training set
#                 DT <- gendata(n = n, p = p, corr = corr, 
#                                     SNR = SNR, parameterIndex = parameterIndex)
#                 
#                 #validation set
#                 DT_val <- gendata(n = 2*n, p = p, corr = corr,
#                                         SNR = SNR, parameterIndex = parameterIndex)
#                 
#                 trainIndex <- drop(caret::createDataPartition(DT$y, p = 0.5, list = FALSE, times = 1))
#                 testIndex <- setdiff(seq_along(DT$y), trainIndex)
#                 
#                 xtrain <- DT$x[trainIndex,,drop=FALSE]
#                 xtest <- DT$x[testIndex,,drop=FALSE]
#                 xvalid <- DT_val$x
#                 
#                 ytrain <- DT$y[trainIndex]
#                 ytest <- DT$y[testIndex]
#                 yvalid <- DT_val$y
#                 
#                 main <- colnames(DT$x)
#                 vnames <- c(main, paste0(main))
#                 # vnames_lasso <- c("E", main) # needs to be in this order for glinternet
#                 
#                 # as is done in pliable lasso, only feed (E,X) to glmnet
#                 X_main <- cbind(E = xtrain)
#                 
#                 # test set
#                 X_main_test <- cbind(E = etest, xtest)
#                 
#                 # validation set
#                 X_main_valid <- cbind(E = evalid, xvalid)
#                 
#                 models[[i]] <- list(xtrain = xtrain, etrain = etrain, ytrain = ytrain, xtrain_lasso = X_main,
#                                     xtest = xtest, etest = etest, ytest = ytest, xtest_lasso = X_main_test,
#                                     xvalid = xvalid, evalid = evalid, yvalid = yvalid, xvalid_lasso = X_main_valid,
#                                     causal = DT$causal, not_causal = DT$not_causal,
#                                     vnames = vnames, vnames_lasso = vnames_lasso)
#               }
#               return(models)
#             })
#   
# }
# 
# data<-make_gendata_Paper_data_split()