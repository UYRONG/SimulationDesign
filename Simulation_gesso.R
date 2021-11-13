# install.packages("devtools",dependencies = TRUE)
#library(devtools)
#devtools::install_github("NataliaZemlianskaia/gesso")

#install.packages("simulator")

library(gesso)
library(dplyr)

## generate the data: 400 main effects and 400 interaction effects 
## with 10 non-zero main effects and 5 non-zero interaction effects, sample size equal to 150
data = data.gen(sample_size=150, p=400, 
                n_g_non_zero=10, n_gxe_non_zero=5, 
                family="gaussian", mode="strong_hierarchical")

## tune the model hyperparameters   
# tune_model = gesso.cv(data$G_train, data$E_train, data$Y_train, 
#                       grid_size=20, tolerance=1e-4,
#                       parallel=TRUE, nfold=4,
#                       normalize=TRUE, normalize_response=TRUE,
#                       seed=1)
model <- gesso.fit(G = data$G_train, E = data$E_train, Y = data$Y_train,
          normalize=TRUE,family = "gaussian", min_working_set_size = 30)
model
coefficients = gesso.coef(fit=gesso_fit$fit, lambda=gesso_fit$lambda_min)
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


