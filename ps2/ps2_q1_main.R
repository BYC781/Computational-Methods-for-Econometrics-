library(tidyverse)
source("./ps2/functions.R")
data <- readxl::read_xlsx("./ps2/CilibertoTamerEconometrica.xlsx")
data <- cbind(data[ ,1:7], rowSums(data %>% select(2:7)), data[ ,8:27])
names(data)[8] <- 'total.N'
n.mkt <- dim(data)[1]
true.n <- cbind(data[,2:8])
get_XX_and_Z.mat_and_Y()

A <- c(1,1)
rho <- 0.5
B <- c(0, 1, 2)
d <- 5
init <- c(1,1,1,1)
fit <- optim(fn = like_oprobit, par = init, method = "Nelder-Mead")

# Print the results
cat("Maximum likelihood estimates:\n")
cat("beta0 =", fit$par[1], "\n")
cat("beta1 =", fit$par[2], "\n")
cat("beta2 =", fit$par[3], "\n")
cat("delta =", fit$par[4], "\n")

# MSM
T <- 10
# get (n.mkt * 7) matrix, each element is pred. number
# col: firm1, firm2, ..., firm6, total
n_hat <- get.n_hat(A,B,d,rho)
g_n <- as.matrix(g_n.fn()) # 21*1

set.seed(2048)
init.param = c(0.5, 0.5, 0.9, 0.1, 0.5, 1.9, 0.6)
msm.fit <- optim(fn = obj.fn, par = init.param, method = "BFGS")
