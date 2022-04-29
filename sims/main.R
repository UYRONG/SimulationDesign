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
devtools::install_github("hugogogo/sprintr", build_vignettes = TRUE)
install_github("gathanei/xyz")
install.packages("RAMP")
pacman::p_load(simulator)
pacman::p_load(gesso)
pacman::p_load(sail)
pacman::p_load(glinternet)
pacman::p_load(hierNet)
pacman::p_load(sprintr)
pacman::p_load(LassoBacktracking)
pacman::p_load(xyz)
pacman::p_load(RAMP)
pacman::p_load(FAMILY)
pacman::p_load(magrittr)
pacman::p_load(stringr)
pacman::p_load(dplyr)
pacman::p_load(tidyr)
pacman::p_load(ggplot2)

source("/Users/rongyu/Desktop/Simulation/sims/eval_functions.R")
source("/Users/rongyu/Desktop/Simulation/sims/method_functions.R")
source("/Users/rongyu/Desktop/Simulation/sims/model_functions.R")
library(data.table)

## @knitr init

name_of_simulation <- "normal-mean-estimation-with-contamination"

# # Simulation Case 1 : All pairs, Continuous Outcome, n=100, p=30 ---------------------------------------------------
sim1 <- new_simulation(name = "2022_sim1",
                       label = "2022_sim1_XX_CASE1",
                       dir = ".") %>%
  generate_model(make_XX_data_split,seed = 1234,
                 n = 100, p=30, corr= 0, SNR = 2, case = 1,
                 lambda.type = "lambda.min") %>%
  
  simulate_from_model(nsim = 4, index = 1:5) %>%
  run_method(list(glinternetsplit_XX),
                  #ramp_split,xyz_split),
             #lassoBTsplit,glinternetsplit_XX,hiernet_split,sprintr_split，xyz_split
             
             parallel = list(socket_names = 35,
                             libraries = c("splines","stringr",
                                           "magrittr","LassoBacktracking","glinternet",
                                           "hierNet","sprintr","xyz","RAMP","FAMILY","simulator","parallel")))
simulator::save_simulation(sim1)

sim1 <- sim1 %>% evaluate(list(msevalid, tpr, fpr, nactive, r2))

sim1 %>% plot_eval(metric_name = "mse")
as.data.frame(evals(sim1))

df <- as.data.frame(evals(sim1))

saveRDS(df, file = "/Users/rongyu/Desktop/Simulation/sims/files/sim-2022_sim1.rds")

df <- readRDS("/Users/rongyu/Desktop/Simulation/sims/files/sim-2022_sim1.rds")
df <- df %>% separate(Model, into = c("simnames","corr","lambda.type","n","p","parameterIndex","SNR_2"),
                      sep = "/")
DT <- as.data.table(df)
DT

# # Simulation Case 2 : All pairs, Continuous Outcome, n=150, p=300 ---------------------------------------------------
# # Simulation Case 3 : All pairs, Continuous Outcome, n=500, p=2000 ---------------------------------------------------


# # Simulation Case 4 : All pairs, Binary Outcome, n=100, p=30 ---------------------------------------------------
sim4 <- new_simulation(name = "2022_sim4",
                       label = "2022_sim4_XX_CASE2",
                       dir = ".") %>%
  generate_model(make_XX_data_split_binary,seed = 1234,
                 n = 100, p=30, corr= 0, SNR = 2, case = 1,
                 lambda.type = "lambda.min") %>%
  
  simulate_from_model(nsim = 6, index = 1:35) %>%
  run_method(list(family_binary_split,ramp_binary_split,glinternet_binary_split_XX,hiernet_binary_split),
             
             parallel = list(socket_names = 35,
                             libraries = c("splines","stringr",
                                           "magrittr","LassoBacktracking","glinternet",
                                           "hierNet","sprintr","xyz","RAMP","FAMILY","simulator","parallel")))
simulator::save_simulation(sim4)

sim4 <- sim4 %>% evaluate(list(msevalid, tpr, fpr, nactive, r2))

sim4 %>% plot_eval(metric_name = "mse")
as.data.frame(evals(sim4))

df <- as.data.frame(evals(sim4))
saveRDS(df, file = "/Users/rongyu/Desktop/Simulation/sims/files/sim-2022_sim4.rds")

# # Simulation Case 5 : All pairs, Binary Outcome, n=150, p=300 ---------------------------------------------------
sim5 <- new_simulation(name = "2022_sim5",
                       label = "2022_sim5_XX_CASE2",
                       dir = ".") %>%
  generate_model(make_XX_data_split_binary,seed = 1234,
                 n = 150, p=300, corr= 0, SNR = 2, case = 2,
                 lambda.type = "lambda.min") %>%
  
  simulate_from_model(nsim = 3, index = 1:20) %>%
  run_method(list(family_binary_split,ramp_binary_split,glinternet_binary_split_XX,hiernet_binary_split),
             
             parallel = list(socket_names = 35,
                             libraries = c("splines","stringr",
                                           "magrittr","LassoBacktracking","glinternet",
                                           "hierNet","sprintr","xyz","RAMP","FAMILY","simulator","parallel")))
simulator::save_simulation(sim5)

sim5 <- sim5 %>% evaluate(list(msevalid, tpr, fpr, nactive, r2))

sim5 %>% plot_eval(metric_name = "mse")
as.data.frame(evals(sim5))

df <- as.data.frame(evals(sim5))
saveRDS(df, file = "/Users/rongyu/Desktop/Simulation/sims/files/sim-2022_sim5.rds")




# # Simulation Case 6 : All pairs, Binary Outcome, n=500, p=2000 ---------------------------------------------------

## @knitr main
# # Simulation Case 7 : GXE, Continuous Outcome, n=100, p=30------------------------
sim7 <- new_simulation(name = "2022_sim7",
                      label = "2022_sim7_XE_CASE1",
                      dir = ".") %>%
  generate_model(make_XE_data_split,seed = 1234,
                 n = 100, p=30, corr= 0, betaE = 2, SNR = 2, case = 1,
                 lambda.type = "lambda.min") %>%

  simulate_from_model(nsim = 1, index = 1:2) %>%
  run_method(list(gessosplit,sailsplit,glinternetsplit_XE),
             parallel = list(socket_names = 35,
                             libraries = c("splines",
                                           "magrittr","gesso","sail","glinternet","simulator","parallel")))
simulator::save_simulation(sim7)

sim7 <- sim7 %>% evaluate(list(msevalid, tpr, fpr, nactive, r2))

sim7 %>% plot_eval(metric_name = "mse")
as.data.frame(evals(sim7))

df <- as.data.frame(evals(sim7))

saveRDS(df, file = "/Users/rongyu/Desktop/Simulation/sims/files/sim-2022_sim7.rds")

df <- readRDS("/Users/rongyu/Desktop/Simulation/sims/files/sim-2022_sim7.rds")
df <- df %>% separate(Model, into = c("simnames","betaE","corr","lambda.type","n","p","parameterIndex","SNR_2"),
                      sep = "/")
DT <- as.data.table(df)
DT

# # --------------------------------------------------------------------------------------------


# # Simulation Case 8 : GXE, Continuous Outcome, n=150, p=300------------------------

sim8 <- new_simulation(name = "2022_sim8",
                       label = "2022_sim8_XE_CASE2",
                       dir = ".") %>%
  generate_model(make_XE_data_split,seed = 1234,
                 n = 150, p=300, corr= 0, betaE = 2, SNR = 2, case = 2,
                 lambda.type = "lambda.min", parameterIndex=list(1), vary_along = "parameterIndex") %>%
  
  simulate_from_model(nsim = 6, index = 1:10) %>%
  run_method(list(gessosplit,sailsplit,glinternetsplit_XE),
             parallel = list(socket_names = 35,
                             libraries = c("splines",
                                           "magrittr","gesso","sail","glinternet","simulator","parallel")))
simulator::save_simulation(sim8)

sim8 <- sim8 %>% evaluate(list(msevalid, tpr, fpr, nactive, r2))

sim8 %>% plot_eval(metric_name = "mse")
as.data.frame(evals(sim8))

df <- as.data.frame(evals(sim8))

saveRDS(df, file = "/Users/rongyu/Desktop/Simulation/sims/files/sim-2022_sim8.rds")

df <- readRDS("/Users/rongyu/Desktop/Simulation/sims/files/sim-2022_sim8.rds")
df <- df %>% separate(Model, into = c("simnames","betaE","corr","lambda.type","n","p","parameterIndex","SNR_2"),
                      sep = "/")
DT <- as.data.table(df)
DT


# # --------------------------------------------------------------------------------------------

# # Simulation Case 9 : GXE, Continuous Outcome, n=500, p=2000------------------------

sim9 <- new_simulation(name = "2022_sim9",
                       label = "2022_sim9_XE_CASE3",
                       dir = ".") %>%
  generate_model(make_XE_data_split,seed = 1234,
                 n = 500, p=2000, corr= 0, betaE = 2, SNR = 2, case = 3,
                 lambda.type = "lambda.min", parameterIndex=list(1), vary_along = "parameterIndex") %>%
  
  simulate_from_model(nsim = 6, index = 1:10) %>%
  run_method(list(gessosplit,sailsplit,glinternetsplit_XE),
             parallel = list(socket_names = 35,
                             libraries = c("splines",
                                           "magrittr","gesso","sail","glinternet","simulator","parallel")))
simulator::save_simulation(sim9)

sim9 <- sim9 %>% evaluate(list(msevalid, tpr, fpr, nactive, r2))

sim9 %>% plot_eval(metric_name = "mse")
as.data.frame(evals(sim9))

df <- as.data.frame(evals(sim9))

saveRDS(df, file = "/Users/rongyu/Desktop/Simulation/sims/files/sim-2022_sim9.rds")

df <- readRDS("/Users/rongyu/Desktop/Simulation/sims/files/sim-2022_sim9.rds")
df <- df %>% separate(Model, into = c("simnames","betaE","corr","lambda.type","n","p","parameterIndex","SNR_2"),
                      sep = "/")
DT <- as.data.table(df)
DT


