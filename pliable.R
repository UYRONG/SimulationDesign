install.packages("/Users/rongyu/Desktop/pliable", repos = NULL, type="source")
library(pliable)









install.packages("devtools")
library(devtools)
install_github("cran/pliable")

n = 20 ; p = 3 ;nz=3
x = matrix(rnorm(n*p), n, p)
mx=colMeans(x)
sx=sqrt(apply(x,2,var))
x=scale(x,mx,sx)
z =matrix(rnorm(n*nz),n,nz)
mz=colMeans(z)
sz=sqrt(apply(z,2,var))
z=scale(z,mz,sz)
y =4*x[,1] +5*x[,1]*z[,3]+ 3*rnorm(n)
fit = pliable(x,z,y)

if (!requireNamespace("pliable", quietly = TRUE))
  install.packages("/Users/rongyu/Desktop/pliable", repos = NULL, type="source")
pliable::install(version = "3.11")
