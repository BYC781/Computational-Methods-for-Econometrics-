% Econ 7218
% Demonstration of Gibbs sampling for linear regression model

ssd = 12345;
rng(ssd);

%% model and estimation setting 

N=500;     % sample size 
T=1000;    % length of MCMC
burn=100;  % length of burn-in

%% DGP parameters

beta=[0.5; 0.5];
sigma=0.5;

%% generate exogenous regressors

x1 = randn(N,1);
x2 = chi2rnd(1,[N,1]);
ep = randn(N,1)*sqrt(sigma);
X = [x1 x2];
Y = X*beta+ep;


%% start Gibbs sampling

beta_T=zeros(2,T);
sigma_T=zeros(T,1);

beta_0=zeros([1,2]);
B_0=eye(2)*3;
rho_0=2.2;
eta_0=0.1;

beta_T(:,1)=20;
sigma_T(1)=10;


for t=2:T
    
    % THE SAMPLING OF BETA FROM PROSTERIOR DISTRIBUTION %
    V=sigma_T(t-1)*eye(N);
    ZVX=X'/V*X;
    ZVY=X'/V*Y;
    beta_T(:,t)=mvnrnd(((inv(B_0)+ZVX)\((B_0)\beta_0'+ZVY)),inv(inv(B_0)+ZVX));
    
    % THE SAMPLING OF SIGMA_E^2 FROM PROSTERIOR DISTRIBUTION %
    
    ep=Y-X*beta_T(:,t);
    rho_1=rho_0+length(ep);
    sigma_T(t)=(ep'*ep+eta_0)/chi2rnd(rho_1);
    
    
end


subplot(2,3,1)
plot(beta_T(1,2:t));figure(gcf);
title('\beta_1')
subplot(2,3,2)
plot(beta_T(2,2:t));figure(gcf);
title('\beta_2')
subplot(2,3,3)
plot(sigma_T(2:t));figure(gcf);
title('\sigma^2')
subplot(2,3,4)
histogram(beta_T(1,2:t));figure(gcf);
title('\beta_1')
subplot(2,3,5)
histogram(beta_T(2,2:t));figure(gcf);
title('\beta_2')
subplot(2,3,6)
histogram(sigma_T(2:t));figure(gcf);
title('\sigma^2')




disp(mean(beta_T(1,burn:t)));
disp(mean(beta_T(2,burn:t)));
disp(mean(sigma_T(burn:t)));

disp(std(beta_T(1,burn:t)));
disp(std(beta_T(2,burn:t)));
disp(std(sigma_T(burn:t)));


