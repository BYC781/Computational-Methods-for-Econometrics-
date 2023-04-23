library(tidyverse)
set.seed(2000)
data <- c(rnorm(2000, 20, 5), rnorm(3000, 10, 3))
df <- data.frame(value = data)

# Plot the histogram
ggplot(df, aes(x = value)) +
    geom_histogram(binwidth = 1, color = "black", fill = "gray") +
    labs(title = "Histogram of Shuffled Mixed Vector",
         x = "Value",
         y = "Count")


# parameters generator
global.par <- function(par){
    mu1 <<- par[1]; sigma1 <<- par[2]; pi1 <<- par[3]
    mu2 <<- par[4]; sigma2 <<- par[5]; pi2 <<- par[6]
}

Pr_zi.eaual.1 <- function(vec){
    pr_vec1 <- dnorm(vec, mean = mu1, sd = sigma1)
    pr_vec2 <- dnorm(vec, mean = mu2, sd = sigma2)
    
    pr_zi.eaual.1 <- pi1*pr_vec1 / (pi1*pr_vec1 + pi2*pr_vec2)
    return(pr_zi.eaual.1)
}

Pr_zi.eaual.2 <- function(vec){
    pr_vec1 <- dnorm(vec, mean = mu1, sd = sigma1)
    pr_vec2 <- dnorm(vec, mean = mu2, sd = sigma2)
    
    pr_zi.eaual.2 <- pi2*pr_vec2 / (pi1*pr_vec1 + pi2*pr_vec2)
    return(pr_zi.eaual.2)
}

obj.fun <- function(vec){
    pr_zi.eaual.1 <- Pr_zi.eaual.1(vec)
    pr_zi.eaual.2 <- Pr_zi.eaual.2(vec)
    q <- pr_zi.eaual.1* (log(pi1) + dnorm(vec, mu1, sigma1, log = TRUE)) + 
        pr_zi.eaual.2* (log(pi2) + dnorm(vec, mu2, sigma2, log = TRUE))
    neg.Q <- -sum(q)
    return(neg.Q)
}

init <- runif(6, 0, 1)
tolerance <- 1e-6; max_iter <- 1e6
count <- 0
record <- rep(NA, max_iter)
global.par(init)
logl <- obj.fun(vec=data)
record[1] <- obj.fun(vec=data)


while (TRUE) {
    count <- count + 1
    # E-step: calculate the posterior probabilities
    pr1 <- Pr_zi.eaual.1(data) ; pr2 <- Pr_zi.eaual.2(data)
    pi1 <- mean(pr1)
    pi2 <- mean(pr2)
    
    # M-step: update the parameters
    mu1 <- sum(data * pr1) / sum(pr1)
    mu2 <- sum(data * pr2) / sum(pr2)
    sigma1 <- sqrt(sum((data - mu1)^2 * pr1) / sum(pr1))
    sigma2 <- sqrt(sum((data - mu1)^2 * pr2) / sum(pr2))
    
    newlogl <- obj.fun(data)
    record[count+1] <- newlogl
    if (abs(newlogl - logl) < tolerance) {
        break
    } else {
        logl <- newlogl
    }

    if (count == max_iter){
        break
    }
}
record <- record[!is.na(record)]
