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