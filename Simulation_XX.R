#devtools::install_github("hugogogo/sprintr", build_vignettes = TRUE)

# X*X Simulation Case-1 for sprintr
library(sprintr)
set.seed(123)
n<-100
p<-100
x<-matrix(data=rnorm(n*p), nrow=n, ncol=p)
# X~N(0,I) main effect beta1=1, beta2=-2, and betaj=0 for j>2.
# Two important interactions are X1*X3 with coeff=3 and X4*X5 = -4
# n=100,p=100
# 
y<-x[,1]-2*x[,2]+3*x[,1]*x[,3]-4*x[,4]*x[,5]+rnorm(100)

fit <-sprinter(x = x, y = y,square = FALSE)

fit$step2[[1]]

estimate<-fit$step3[[4]]$coef[,30]

fit_csvstep1<-sprinter(x=x,y=y,square=FALSE,cv_step1 = TRUE)
#print(fit,which=2)
#plot(fit,which=3)
fit_cv<-cv.sprinter(x=x,y=y,square=FALSE)
fit_cv$compact

print(fit_cv)
plot(fit_cv) #not work

# Prediction
# Defined for both object that computes the prediction for a new data matrix of main effects
newdata<-matrix(rnorm(20*p),nrow=20,ncol=p)
pred<-predict(fit,newdata=newdata)
print(pred)
pred_cv<-predict(fit_cv,newdata=newdata)



                    