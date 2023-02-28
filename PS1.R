set.seed(1234)
library(evd)
## R=1, N=400
b1 <- 0.5; b2 <- -0.5; n =400
x1 <- rnorm(n, 0, 1)
x2 <- rchisq(n, 1)
u1 <- rgumbel(n)
u2 <- rgumbel(n)

# construct y
y <- as.numeric((x1 + u1) > (x2 + u2))

## log_likelihood function
loglik <- function(beta1, beta2){
    index <- beta1*x1 - beta2 * x2
    sum(y*(index - log(1+exp(index))) - (1-y)* log(1+exp(index)))
}

##### grid search #####
beta1_grid <- seq(from = -5, to = 5, by = 0.1)
beta2_grid <- seq(from = -5, to = 5, by = 0.1)
max_lik <- -10000000
argmax_beta <- c(0,0)
for (i in beta1_grid){
    for(j in beta2_grid){
        temp <- loglik(i,j)
        if (temp >= max_lik){ 
            max_lik <- temp
            argmax_beta <- c(i,j)
        }
    }
}
argmax_beta


##### BHHH 跑不出來啦幹#####
## R=100, N=400
R=100;N=400
X1 <- matrix(rnorm(R*N), nrow=N)
X2 <- matrix(rchisq(R*N, 1), nrow=N)
U1 <- matrix(rgumbel(R*N), nrow=N)
U2 <- matrix(rgumbel(R*N), nrow=N)
Y <- as.matrix((X1 + U1) > (X2 + U2)) %>% ifelse(1,0)

beta_seq <- matrix(0, nrow=R, ncol=2)
for (i in 1:R){
    gradient <- function(beta1, beta2) {
        p <- exp(beta1*X1[,i] - beta2*X2[,i]) /(1+exp(beta1*X1[,i] - beta2*X2[,i]))
        g1 <- sum(Y[,i]*X1[,i] - X1[,i]*p)
        g2 <- sum(Y[,i]*X2[,i] - X2[,i]*p)
        return(matrix(c(g1, g2), ncol=1))
    }
    hessian <- function(beta1, beta2) {
        p <- exp(beta1*X1[,i] - beta2*X2[,i]) /(1+exp(beta1*X1[,i] - beta2*X2[,i]))
        H11 <- sum(x1^2*p*(1-p))
        H12 <- sum(x1*x2*p*(1-p))
        H22 <- sum(x2^2*p*(1-p))
        return(matrix(c(H11, H12, H12, H22), ncol=2))
    }
    beta <- c(0, 0)
    tol <- 1e-6
    maxiter <- 100
    for (j in 1:maxiter) {
        g <- gradient(beta[1], beta[2])
        H <- hessian(beta[1], beta[2])
        BHHH <- g %*% t(g)
        if (max(abs(g)) < tol) {
            break
        }
        beta <- beta + solve(BHHH) %*% g
    }
    beta_seq[i, 1] <- beta[1]
    beta_seq[i, 2] <- beta[2]
}


##################### q2-1 #####################
library(dplyr)
df <- readxl::read_xlsx('cps09mar.xlsx')
df$married <- ifelse(df$marital %in% c(1, 2, 3), 1, 0)
blk_wm_midwest <- df %>% 
    filter(race==2, region == 2, female == 1)
blk_wm_midwest_logit <- glm(married ~ age + I(age^2) + education, family = binomial, data = blk_wm_midwest)


summary(blk_wm_midwest_logit)$coefficients[,1:2] ## coef and std. error


##################### q2-2 #####################
# 設定 Bootstrap 樣本的大小和重複次數
n <- nrow(blk_wm_midwest)
B <- 1000

# 創建一個空矩陣來保存 Bootstrap 標準誤
se_boot_vec <- matrix(NA, nrow=B, ncol=4)

# 進行 Bootstrap
for (i in 1:B) {
    # 從原始資料集中隨機取樣 n 個觀察值
    sample_idx <- sample(1:n, size = n, replace = TRUE)
    sample_data <- blk_wm_midwest[sample_idx, ]
    
    # 使用取樣的資料跑logit
    fit <- glm(married ~ age + I(age^2) + education, data = sample_data, family = binomial())
    
    # 提取sd
    for(j in 1:4){
        se_boot_vec[i,j] <- summary(fit)$coefficients[j, 2]
    }
}

# 計算 Bootstrap 標準誤
se_boot <- apply(se_boot_vec, 2, sd)

# 印出 Bootstrap 標準誤
se_boot

##################### q2-3 not sure#####################

# 計算 - (b1) / (2*b2) 的值
beta_hat <- coef(blk_wm_midwest_logit)
result <- - beta_hat[2] / (2 * beta_hat[3])

# 計算偏導數
d1 <- -1 / (2 * beta_hat[3])
d2 <- beta_hat[2] / (2 * beta_hat[3]^2)

# 計算標準誤
vcov_matrix <- vcov(blk_wm_midwest_logit)
delta_se <- sqrt(d1^2 * vcov_matrix[2,2] + d2^2 * vcov_matrix[3,3] + 2 * d1 * d2 * vcov_matrix[2,3])

# 輸出結果
unname(delta_se)


##################### q2-4 ##############################
# 設定 Bootstrap 重複次數
B <- 1000

# 儲存 Bootstrap 估計量的向量
t_boot <- numeric(B)

# 執行 Bootstrap
for (i in 1:B) {
    # 從原始資料集中抽樣，構建 Bootstrap 樣本
    sample_data_boot <- blk_wm_midwest[sample(nrow(blk_wm_midwest), replace = TRUE), ]
    
    # 在 Bootstrap 樣本上擬合迴歸模型
    fit_boot <- glm(married ~ age + I(age^2) + education, data = sample_data_boot, family = binomial())
    
    # 取得迴歸係數
    coef_boot <- coef(fit_boot)
    
    # 計算 Bootstrap 估計量
    t_boot[i] <- -coef_boot[2] / (2 * coef_boot[3])
}

# 計算 Bootstrap 估計量的標準差
se_boots <- sd(t_boot)

# 印出 Bootstrap 標準差
cat("Bootstrap standard error is", se_boots, "\n")

