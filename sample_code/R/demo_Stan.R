## Demonstration of Rstan with logit regression ##
## reference: https://mc-stan.org/docs/2_22/stan-users-guide/index.html ##
## To run Rstan, Rtools (download from https://cran.r-project.org/bin/windows/Rtools/rtools42/rtools.html) are also needed ##

install.packages("rstan")
install.packages("bayesplot")
library("bayesplot")
library("rstan")


N=1000;
beta_0=-1;
beta_1=0.5;
X   <- rnorm(N,mean=0,sd=2);
Y   <- rep(0, N)
for (i in 1:N ) {
    Xb=beta_0+beta_1*X[i]
    p=exp(Xb)/(1+exp(Xb))
    if (runif(1)<=p){
    Y[i]=1
    } 
}


smoke.stan = "
data {
  int N;          // number of trials
  vector[N] x;    // predictor matrix
  int y[N];       // outcome vector
}

parameters {
  real alpha;
  real beta; 
}

model {
  y ~ bernoulli_logit(beta*x+alpha); // likelihood
}
"

fit = stan(model_code=smoke.stan, 
      data=list(y=Y,x=X,N=N), iter=3000)

print(fit, probs=c(0.1, 0.9))

# Extracting the posterior draws
posterior = as.array(fit)

dimnames(posterior)

# Plotting MCMC draws using bayesplot package
# referece: https://mc-stan.org/bayesplot/articles/plotting-mcmc-draws.html

mcmc_trace(posterior ,  pars = c("alpha","beta"), n_warmup = 300)

color_scheme_set("red")

mcmc_intervals(posterior, pars = c("alpha", "beta"))

mcmc_areas(posterior,
  pars =  c("alpha", "beta"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)




