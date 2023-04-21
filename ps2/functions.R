get_XX_and_Z.mat_and_Y <- function(...){
    X1 <- data$marketdistance
    X2 <- data$marketsize
    XX <<- cbind(rep(1, n.mkt), X1, X2)
    
    Z1 <- cbind(data[ ,17], data[,23])
    Z2 <- cbind(data[ ,18], data[,24])
    Z3 <- cbind(data[ ,19], data[,25])
    Z4 <- cbind(data[ ,20], data[,26])
    Z5 <- cbind(data[ ,21], data[,27])
    Z6 <- cbind(data[ ,22], data[,28])
    Z.mat <<- cbind(Z1, Z2, Z3, Z4, Z5, Z6)
    
    Y <<- data$total.N
}

draw_u <- function(...){
    u_ik <<- matrix(rnorm(n.mkt * 6), ncol=6)
    u_i0 <<- rnorm(n.mkt)
}



like_oprobit<- function(init){
    beta.mat <- init[1:3]
    delta <- init[4]
    f <- 0
    N <- length(Y)
    for (i in 1:N){
        if (Y[i] == 0){
            p = pnorm(-t(XX[i,])%*% beta.mat)
        }
        else if (Y[i] == 6){
            p = 1 - pnorm(-t(XX[i,]) %*% beta.mat + delta*log(Y[i]))
        }
        else{
            p = pnorm(-t(XX[i,]) %*% beta.mat + delta*log(Y[i]+1)) - pnorm(-t(XX[i,]) %*% beta.mat + delta*log(Y[i]))
        }
    }
    f <- f -log(p)
}



single.sim.process <- function(A, B, d, rho){
    n_pred <- matrix(nrow=n.mkt, ncol=7)
    for (mkt.i in 1:n.mkt){ # 1:n.mkt
        is.in <- 0
        Z <- Z.mat[mkt.i, ] %>% matrix(ncol = 2, byrow = TRUE)
        
        for (n.try in 0:6){
            profit <- Z %*% A + sqrt(1-rho^2) * u_ik[mkt.i, ] + as.numeric(XX[mkt.i, ] %*% B) + rho * u_i0[mkt.i] - log(n.try)*d
            
            total.in <- sum(profit > 0)
            
            if (total.in < n.try){
                n_pred[mkt.i, 7] <- n.try - 1
                n_pred[mkt.i, 1:6] <- is.in
                break
            }
            is.in <- as.numeric(profit > 0)
            if (total.in == 6 & n.try == 6){
                n_pred[mkt.i, 7] <- n.try
                n_pred[mkt.i, 1:6] <- is.in
            }
        }
        
    }
    return(n_pred)
}

get.n_hat <- function(A,B,d,rho){
    container <- matrix(0, nrow=n.mkt, ncol=7)
    for (t in 1:T){
        draw_u()
        container <- container + single.sim.process(A,B,d,rho)
    }
    n_hat<- container / T
    return(n_hat)
}

pred.error <- function(n_hat){
    v <- true.n - n_hat
    return(v)
}

g_n.fn <- function(v){
    v <- as.matrix(v)
    new.Z.mat <- cbind(rep(1,nrow(Z.mat)),Z.mat)
    E_v0X  = v[,7] * XX             
    E_v1Z  = v[,1] * new.Z.mat[,c(1,2,3)]  
    E_v2Z  = v[,2] * new.Z.mat[,c(4,5)] 
    E_v3Z  = v[,3] * new.Z.mat[,c(6,7)]      
    E_v4Z  = v[,4] * new.Z.mat[,c(8,9)]
    E_v5Z  = v[,5] * new.Z.mat[,c(10,11)]
    E_v6Z  = v[,6] * new.Z.mat[,c(12,13)]
    
    v_i_hat <- cbind(E_v0X, E_v1Z, E_v2Z, E_v3Z, E_v4Z, E_v5Z, E_v6Z)
    g_n <- colMeans(v_i_hat)
    return(g_n)
}

get_ABDrho <- function(para){
    A <- para[1:2]
    B <- para[3:5]
    d <- para[6]
    rho <- para[7]
}

obj.fn <- function(para){
    get_ABDrho(para)
    n_hat <- get.n_hat(A,B,d,rho)
    v <- pred.error(n_hat)
    g_n <- g_n.fn(v)
    mom <- t(g_n) %*% g_n
    return(mom)
}



