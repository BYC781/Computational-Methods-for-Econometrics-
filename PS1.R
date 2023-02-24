set.seed(7218)
## R=1, N=400
b1 <- 0.5; b2 <- -0.5; n =400
x1 <- rnorm(n, 0, 1)
x2 <- rchisq(n, 1)
u1 <- rgumbel(n)
u2 <- rgumbel(n)

p <- exp(b1*x1 - b2*x2) / (1+exp(b1*x1 - b2*x2))

# construct y
y <- rep(0, 400)
y <- as.numeric((x1 + u1) > (x2 + u2)))

## log_likelihood function
loglik <- function(beta1, beta2){
    index <- beta1*x1 - beta2 * x2
    sum(y*(index - log(1+exp(index))) - (1-y)* log(1+exp(index)))
}

## grid search
beta1_grid <- seq(from = -5, to = 5, by = 0.5)
beta2_grid <- seq(from = -5, to = 5, by = 0.5)
max_lik <- -1000
max_beta <- c(0,0)
for (i in beta1_grid){
    for(j in beta2_grid){
        temp <- loglik(i,j)
        if (temp >= max_lik){ 
            max_lik <- temp
            max_beta <- c(i,j)
        }
    }
}
max_beta



## R=100, N=400
R=100;N=400
X1 <- matrix(rnorm(R*N, 0, 1), nrow=N)
X2 <- matrix(rchisq(R*N, 1), nrow=N)
U1 <- matrix(rgumbel(R*N), nrow=N)
U2 <- matrix(rgumbel(R*N), nrow=N)
Y <- matrix(0, nrow=N, ncol=R)
Y <- as.matrix((X1 + U1) > (X2 + U2))
Y <- ifelse(Y, 1, 0)
