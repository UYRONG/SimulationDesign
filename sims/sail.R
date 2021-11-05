if (!require("pacman")) install.packages("pacman")
pacman::p_install_gh("sahirbhatnagar/sail")

f.basis <- function(x) splines::bs(x, degree = 5)

print(head(f.basis))