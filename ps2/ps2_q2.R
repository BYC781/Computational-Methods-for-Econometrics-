library(tidyverse)
# DGP
set.seed(123)
x <- c(rnorm(3000, 20, 5), rnorm(2000, 10, 3))
df <- data.frame(x = x)
ggplot(df, aes(x = x)) +
    geom_histogram(binwidth = 1, color = "gray", fill = "gray") +
    labs(title = "Histogram of Mixed Vector",
         x = "Value",
         y = "Count")

# init
mu <- c(8, 15)
sigma <- c(5, 5)
pi <- c(0.5, 0.5)


# E-step
estep <- function(x, mu, sigma, pi) {
    n <- length(x)
    k <- length(mu)
    post <- matrix(0, n, k)
    for (i in 1:n) {
        for (j in 1:k) {
            post[i, j] <- dnorm(x[i], mu[j], sigma[j]) * pi[j]
        }
        post[i, ] <- post[i, ] / sum(post[i, ])
    }
    return(post)
}


# M-step
mstep <- function(x, post) {
    n <- nrow(post)
    k <- ncol(post)
    mu <- numeric(k)
    sigma <- numeric(k)
    pi <- numeric(k)
    for (j in 1:k) {
        mu[j] <- sum(post[, j] * x) / sum(post[, j])
        sigma[j] <- sqrt(sum(post[, j] * (x - mu[j])^2) / sum(post[, j]))
        pi[j] <- sum(post[, j]) / n
    }
    return(list(mu = mu, sigma = sigma, pi = pi))
}



# EM
em <- function(x, mu, sigma, pi, tol = 1e-6, maxiter = 100) {
    loglik <- numeric(maxiter)
    for (iter in 1:maxiter) {
        # E-step
        post <- estep(x, mu, sigma, pi)
        # M-step
        params <- mstep(x, post)
        # update parameters
        mu <- params$mu
        sigma <- params$sigma
        pi <- params$pi
        # calculate log-likelihood
        loglik[iter] <- sum(post[,1]* log(pi[1]) + dnorm(x, mu[1], sigma[1], log=TRUE)+
                                post[,2]* log(pi[2]) + dnorm(x, mu[2], sigma[2], log=TRUE))
        # check convergence
        if (iter > 1 && abs(loglik[iter] - loglik[iter - 1]) < tol) {
            break
        }
    }
    return(list(mu = mu, sigma = sigma, pi = pi, loglik = loglik[1:iter]))
}

# EM algorithm
result <- em(x, mu, sigma, pi)

# output
cat("mu:", result$mu, "\n")
cat("sigma:", result$sigma, "\n")
cat("pi:", result$pi, "\n")


df$y <- result$pi[1]*dnorm(df$x, mean = result$mu[1], sd = result$sigma[1])
df$z <- result$pi[2]*dnorm(df$x, mean = result$mu[2], sd = result$sigma[2])

ggplot(df, aes(x = x)) + 
    geom_line(aes(y = y), color = "blue") +
    geom_line(aes(y = z), color = "red") + 
    geom_histogram(aes(y = ..density..), bins = 50, alpha = 0.5) + 
    theme_minimal()

