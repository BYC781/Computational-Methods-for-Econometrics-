library(tidyverse)
library(kable)
library(kableExtra)
source("./ps2/ps2_q1_fn.R")
flight.df <- readxl::read_xlsx("./ps2/CilibertoTamerEconometrica.xlsx")
dat2 <- (flight.df %>% 
            mutate(N = rowSums(across(airlineAA:airlineWN)), const = 1, .before = marketdistance))

firm.Names <- flight.df %>% select(starts_with("airline")) %>% colnames() %>% substring(8)

market.regressors <- c('const',"marketdistance","marketsize" )
firm.regressors <- c("marketpresence","mindistancefromhub")

market.regressors.i <- c(9, 10, 15)
firm.regressors.i <- c(18, 24)

n.mkt <- nrow(dat2)
firm.n <- length(firm.Names)

length.of.each.param <- c(length(firm.regressors),length(market.regressors),1,1 )
get_XX_and_Z.mat_and_Y()



init <- c(1,1,1,1)
fit <- optim(fn = like_oprobit, par = init, method = "BFGS")

# Print the results
cat("Maximum likelihood estimates:\n")
cat("beta0 =", fit$par[1], "\n")
cat("beta1 =", fit$par[2], "\n")
cat("beta2 =", fit$par[3], "\n")
cat("delta =", fit$par[4], "\n")

# MSM
T <- 100
set.seed(2048)
init.param = c(0.5, 0.5, 0.9, 0.1, 0.5, 1.9, 0.6)
A <- init.param[1:2] %>% as.matrix()
B <- init.param[3:5] %>% as.matrix()
d <- init.param[6]
rho <- init.param[7]

n_hat <- get.n_hat(A,B,d,rho)
v <- pred.error(n_hat)
g_n <- g_n.fn(v)
mom <- t(g_n) %*% g_n

msm.fit <- optim(fn = obj.fn, par = init.param, method = "BFGS")
kable(msm.fit[["par"]])
