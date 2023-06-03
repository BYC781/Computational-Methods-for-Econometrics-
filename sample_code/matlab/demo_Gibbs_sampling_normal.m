% Econ 7218
% Demonstration of Gibbs sampling for mean and variance from a normal
% distribution
clear;
ssd = 12345;
rng(ssd);

%% model and estimation setting 

N=500;     % sample size 
T=1000;    % length of MCMC
burn=100;  % length of burn-in

%% DGP parameters

mu=2.0;
sigma=1.0;

%% generate exogenous regressors


Y = randn(N,1)*sqrt(sigma)+mu;


%% start Gibbs sampling

mu_T=zeros(T,1);
sigma_T=zeros(T,1);

mu_0=0;
rho_0=0.1;
delta_0=2.2;
lambda_0=1;

mu_T(1)=10;
sigma_T(1)=20;


for t=2:T
    
    % THE SAMPLING OF MU FROM PROSTERIOR DISTRIBUTION %
    V=sigma_T(t-1);
    mu_T(t)=mvnrnd(((lambda_0^(-1)+N*V^(-1))\(lambda_0^(-1)*mu_0+N*V^(-1)*mean(Y))),inv(lambda_0^(-1)+N*V^(-1)));
    
    % THE SAMPLING OF SIGMA^2 FROM PROSTERIOR DISTRIBUTION %
    
    ep=Y-mu_T(t);
    delta_1=delta_0+N;
    sigma_T(t)=(ep'*ep+rho_0)/chi2rnd(delta_1);
        
end


subplot(2,1,1)
plot(mu_T(1:t));figure(gcf);
title('\mu')
subplot(2,1,2)
plot(sigma_T(1:t));figure(gcf);
title('\sigma^2')

figure;

subplot(2,1,1)
histogram(mu_T(1:t));figure(gcf);
title('\mu')
subplot(2,1,2)
histogram(sigma_T(1:t));figure(gcf);
title('\sigma^2')





