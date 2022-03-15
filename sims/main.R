# This is the main simulator file

#library(simulator) # this file was created under simulator version 0.2.3
# If install packages simulator caused warning of the package is not available
# for the current version of R, try to upgrade the R version or download the package from GitHub.

# Regular ways of download the package
# install.packages("simulator") 

# Download the package from Github
devtools::install_github("jacobbien/simulator")
library(simulator)

install.packages('BiocManager')
pacman::p_load(simulator)
pacman::p_load(gesso)
pacman::p_load(sail)
pacman::p_load(glinternet)
pacman::p_load(magrittr)
pacman::p_load(dplyr)
pacman::p_load(tidyr)
pacman::p_load(ggplot2)

source("/Users/rongyu/Desktop/Simulation/sims/eval_functions.R")
source("/Users/rongyu/Desktop/Simulation/sims/method_functions.R")
source("/Users/rongyu/Desktop/Simulation/sims/model_functions.R")

## @knitr init

name_of_simulation <- "normal-mean-estimation-with-contamination"

## @knitr main
# # Simulation Case 7 : GXE, Continuous Outcome, n=100, p=30------------------------
sim <- new_simulation(name = "2022",
                      label = "2022_try",
                      dir = ".") %>%
  generate_model(make_XE_data_split,seed = 1234,
                 n = 100, p=30, corr= 0, betaE = 2, SNR = 2,
                 lambda.type = "lambda.min",parameterIndex = list(1),vary_along = "parameterIndex") %>%

  simulate_from_model(nsim = 1, index = 1:2) %>%
  run_method(list(gessosplit,sailsplit,glinternetsplit_XE),
             parallel = list(socket_names = 35,
                             libraries = c("splines",
                                           "magrittr","gesso","sail","glinternet","simulator","parallel")))
simulator::save_simulation(sim)

sim <- sim %>% evaluate(list(msevalid, tpr, fpr, nactive, r2))

sim %>% plot_eval(metric_name = "mse")
as.data.frame(evals(sim))

# # --------------------------------------------------------------------------------------------

# simulator::save_simulation(sim)
# sim
# ## Try 
# sim <- new_simulation(name = "2021",
#                       label = "2021_try",
#                       dir = ".") %>%
#   generate_model(make_gendata_Paper_data_split, seed = 1234,
#                  n = 400, p = 1000, corr = 0, betaE = 2, SNR = 2, lambda.type = "lambda.min",
#                  parameterIndex = list(1,2,3,4,5),
#                  vary_along = "parameterIndex") %>%
#   simulate_from_model(nsim = 1, index = 1:35) %>%
#   run_method(list(sailsplitlinear, lassoBTsplit),
#              parallel = list(socket_names = 35,
#                              libraries = c("LassoBacktracking","splines",
#                                            "magrittr","sail","simulator", "parallel")))
# ## @knitr plots
# 
# # plot_eval_by(sim, "hisloss", varying = "prob")
# 
# ## @knitr tables
# 
# # tabulate_eval(sim, "herloss", output_type = "markdown",format_args = list(digits = 1))
# 
# # Simulation Case 1 : All pairs, Continuous Outcome, n=100, p=30 ---------------------------------------------------
# # Simulation Case 2 : All pairs, Continuous Outcome, n=150, p=300 ---------------------------------------------------
# # Simulation Case 3 : All pairs, Continuous Outcome, n=500, p=2000 ---------------------------------------------------
# # Simulation Case 4 : All pairs, Binary Outcome, n=100, p=30 ---------------------------------------------------
# # Simulation Case 5 : All pairs, Binary Outcome, n=150, p=300 ---------------------------------------------------
# # Simulation Case 6 : All pairs, Binary Outcome, n=500, p=2000 ---------------------------------------------------
# # Simulation Case 7 : GXE, Continuous Outcome, n=100, p=30 ---------------------------------------------------
# sim <- new_simulation(name = "Simulation_Case_7",
#                       label = "Simulation_Case_7",
#                       dir = ".") %>%
#   generate_model(make_gendata_Paper_data_split, seed = 1234,
#                  n = 100, p = 30, corr = 0, betaE = 2, SNR = 2, lambda.type = "lambda.min",
#                  parameterIndex = list(1,2,3,4,5),
#                  vary_along = "parameterIndex") %>%
#   simulate_from_model(nsim = 1, index = 1:35) %>%
#   run_method(list(sailsplitlinear),
#              parallel = list(socket_names = 35,
#                              libraries = c("gesso","glinternet","magrittr","sail","simulator", "parallel")))
# 
# # Simulation Case 8 : GXE, Continuous Outcome, n=150, p=300 ---------------------------------------------------
# sim <- new_simulation(name = "Simulation Case 7",
#                       label = " n=150, p=300",
#                       dir = ".") %>%
#   generate_model(make_gendata_Paper_data_split, seed = 1234,
#                  n = 150, p = 300, corr = 0, betaE = 2, SNR = 2, lambda.type = "lambda.min",
#                  parameterIndex = list(1,2,3,4,5),
#                  vary_along = "parameterIndex") %>%
#   simulate_from_model(nsim = 1, index = 1:35) %>%
#   run_method(list(sailsplitlinear,gessosplit,GLinternetsplit),
#              parallel = list(socket_names = 35,
#                              libraries = c("gesso","glinternet","magrittr","sail","simulator", "parallel")))
# 
# # Simulation Case 9 : GXE, Continuous Outcome, n=500, p=2000 ---------------------------------------------------
# sim <- new_simulation(name = "Simulation Case 7",
#                       label = " n=500, p=2000",
#                       dir = ".") %>%
#   generate_model(make_gendata_Paper_data_split, seed = 1234,
#                  n = 500, p = 2000, corr = 0, betaE = 2, SNR = 2, lambda.type = "lambda.min",
#                  parameterIndex = list(1,2,3,4,5),
#                  vary_along = "parameterIndex") %>%
#   simulate_from_model(nsim = 1, index = 1:35) %>%
#   run_method(list(sailsplitlinear,gessosplit,GLinternetsplit),
#              parallel = list(socket_names = 35,
#                              libraries = c("gesso","glinternet","magrittr","sail","simulator", "parallel")))
# 
