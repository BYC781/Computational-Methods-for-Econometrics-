get_XX_and_Z.mat_and_Y <- function(...){
    XX <<- dat2[, market.regressors]
    
    Z1 <<- dat2[, firm.regressors.i]
    Z2 <<- dat2[, firm.regressors.i +1]
    Z3 <<- dat2[, firm.regressors.i +2]
    Z4 <<- dat2[, firm.regressors.i +3]
    Z5 <<- dat2[, firm.regressors.i +4]
    Z6 <<- dat2[, firm.regressors.i +5]
    Z.mat <<- cbind(Z1, Z2, Z3, Z4, Z5, Z6)
    
    Y <<- dat2$N
    
    true.n <<- dat2 %>% select(airlineAA:N)
}

draw_u <- function(...){
    u_ik <<- matrix(rnorm(n.mkt * 6), ncol=6)
    u_i0 <<- rnorm(n.mkt)
}



like_oprobit<- function(init){
    XX <- XX %>% as.matrix()
    beta.mat <- init[1:3] %>% as.matrix(, ncol=1)
    delta <- init[4]
    f <- 0
    for (i in 1:n.mkt){
        if (Y[i] == 0){
            p = pnorm(-XX[i,]%*% beta.mat)
        }
        else if (Y[i] == 6){
            p = 1 - pnorm(-XX[i,] %*% beta.mat + delta*log(Y[i]))
        }
        else{
            p = pnorm(-XX[i,] %*% beta.mat + delta*log(Y[i]+1)) - 
                pnorm(-XX[i,] %*% beta.mat + delta*log(Y[i]))
        }
    }
    f <- f -log(p)
    return(f)
}


sim.process <- function(mkt.i, A, B, d, rho){
    n_pred <- matrix(nrow=1, ncol=7)
    is.in <- 0
    Z <- Z.mat[mkt.i, ] %>% as.numeric() %>%  matrix(ncol = 2, byrow = TRUE)
    rho <- 1 / (1 + exp(-rho)) # keep rho between 0 and 1
    
    for (n.try in 0:6){
        profit <- Z %*% A + sqrt(1-rho^2) * u_ik[mkt.i, ] + as.numeric((XX[mkt.i, ] %>% as.numeric()) %*% B) + rho * u_i0[mkt.i] - log(n.try)*d
        
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
    return(n_pred)
}

single.sim.process <- function(A, B, d, rho){
    n_pred <- matrix(nrow=n.mkt, ncol=7)
    for (mkt.i in 1:n.mkt){ # 1:n.mkt
        is.in <- 0
        Z <- Z.mat[mkt.i, ] %>% as.numeric() %>%  matrix(ncol = 2, byrow = TRUE)
        rho <- 1 / (1 + exp(-rho)) # keep rho between 0 and 1
        
        for (n.try in 0:6){
            profit <- Z %*% A + sqrt(1-rho^2) * u_ik[mkt.i, ] + as.numeric((XX[mkt.i, ] %>% as.numeric()) %*% B) + rho * u_i0[mkt.i] - log(n.try)*d
            
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
    E_v0X  = v[,7] * XX             
    E_v1Z  = v[,1] * cbind(XX, Z1)
    E_v2Z  = v[,2] * cbind(XX, Z2)
    E_v3Z  = v[,3] * cbind(XX, Z3)
    E_v4Z  = v[,4] * cbind(XX, Z4)
    E_v5Z  = v[,5] * cbind(XX, Z5)
    E_v6Z  = v[,6] * cbind(XX, Z6)
    
    v_i_hat <- cbind(E_v0X, E_v1Z, E_v2Z, E_v3Z, E_v4Z, E_v5Z, E_v6Z)
    g_n <- colMeans(v_i_hat)
    return(g_n)
}

obj.fn <- function(para){
    A <- para[1:2] %>% as.matrix()
    B <- para[3:5] %>% as.matrix()
    d <- para[6]
    rho <- para[7]
    n_hat <- get.n_hat(A,B,d,rho)
    v <- pred.error(n_hat)
    g_n <- g_n.fn(v)
    mom <- t(g_n) %*% g_n
    return(mom)
}


