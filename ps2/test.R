library(foreach)
library(doParallel)
library(tidyverse)
library(kableExtra)

source("./ps2/funs.R")
flight.df <- readxl::read_xlsx("./ps2/CilibertoTamerEconometrica.xlsx")
dat2 <- (flight.df %>% 
             mutate(N = rowSums(across(airlineAA:airlineWN)), const = 1, .before = marketdistance))

market.regressors.i <- c(9, 10, 15)
firm.regressors.i <- c(18, 24)
n.mkt <- nrow(dat2)
get_XX_and_Z.mat_and_Y()




n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
clusterEvalQ(cl=my.cluster, source("./ps2/funs.R"))
clusterEvalQ(cl=my.cluster, library(magrittr))
clusterExport(cl = my.cluster, c("n.mkt","Z1", "Z2", "Z3","Z4", "Z5", "Z6", "Z.mat", "XX", "true.n"))


T <- 100
init.param = c(0.5, 0.5, 0.9, 0.1, 0.5, 1.9, 0.6)
msm.fit <- optim(fn = obj.fn, par = init.param, method = "BFGS")


stopCluster(my.cluster)
