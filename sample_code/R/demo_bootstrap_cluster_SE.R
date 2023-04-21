 install.packages('multiwayvcov')
 install.packages('lmtest')
 library(multiwayvcov)
 library(lmtest) 
# Artificial balanced panel data set from Petersen (2009) 
# https://www.kellogg.northwestern.edu/faculty/petersen/htm/papers/standarderror.html
# for illustrating and benchmarking clustered standard errors. 
 data(petersen)
 m1 <- lm(y ~ x, data = petersen)
 summary(m1)
 
# Cluster by firm
 boot_firm <- cluster.boot(m1, petersen$firmid)
 coeftest(m1, boot_firm)

# Cluster by year
 boot_year <- cluster.boot(m1, petersen$year)
 coeftest(m1, boot_year)
 
# Double cluster by firm and year
 boot_both <- cluster.boot(m1, cbind(petersen$firmid, petersen$year))
 coeftest(m1, boot_both)
 

# Go multicore using the parallel package
 require(parallel)
 cl <- makeCluster(4)
 options(boot.ncpus = 4)
 boot_both <- cluster.boot(m1, cbind(petersen$firm, petersen$year), parallel = cl)
 stopCluster(cl)
 coeftest(m1, boot_both)
 
# Go multicore using the parallel package, let boot handle the parallelization
 require(parallel)
 options(boot.ncpus = 8)
 boot_both <- cluster.boot(m1, cbind(petersen$firm, petersen$year), parallel = TRUE)
 coeftest(m1, boot_both)
