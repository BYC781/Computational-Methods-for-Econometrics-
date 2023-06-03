% Econ 7218
% Demonstration of Gibbs sampling for linear regression model

ssd = 12345;
rng(ssd);

%% model and estimation setting

N=300;    % sample size
T=1000;    % length of MCMC
burn=100;  % length of burn-in

%% DGP parameters

beta=[0.5; 0.5];


%% generate exogenous regressors

x1 = randn(N,1)*2;
x2 = chi2rnd(2,[N,1]);
X = [x1 x2];
P = normcdf(X*beta);
Y = rand(N,1)<=P;


%% Gibbs sampling with latent variable augmentation

tic;

Z=zeros(N,1);
beta_T=zeros(2,T);

beta_0=zeros([1,2]);
B_0=eye(2)*3;

beta_T(:,1)=0;

for t=2:T
    
    % THE SAMPLING OF LATENT VARIABLE Z FROM CONDITIONAL PROSTERIOR DISTRIBUTION %
    
    for i=1:N
        mu=X*beta_T(:,t-1);
        
        if Y(i)==1
            temp=trandn((0-mu(i)),(inf-mu(i)));
            Z(i)=mu(i)+temp;
        else
            temp=trandn((-inf-mu(i)),(0-mu(i)));
            Z(i)=mu(i)+temp;
        end
    end
    
    % THE SAMPLING OF BETA FROM CONDITIONAL PROSTERIOR DISTRIBUTION %
    
    XX=X'*X;
    XZ=X'*Z;
    beta_T(:,t)=mvnrnd(((inv(B_0)+XX)\((B_0)\beta_0'+XZ)),inv(inv(B_0)+XX));
    
end

toc;
subplot(2,1,1)
plot(beta_T(1,burn:t));figure(gcf);
title('\beta_1')
subplot(2,1,2)
plot(beta_T(2,burn:t));figure(gcf);
title('\beta_2')


disp(mean(beta_T(1,burn:t)));
disp(mean(beta_T(2,burn:t)));

disp(std(beta_T(1,burn:t)));
disp(std(beta_T(2,burn:t)));

%{
%% M-H sampling 
tic;

beta_T=zeros(2,T);

acc_rate=zeros(T,1);
acc=0;

beta_0=zeros([1,2]);
B_0=eye(2)*3;

beta_T(:,1)=0;

for t=2:T
         
    % THE SAMPLING OF BETA FROM CONDITIONAL PROSTERIOR DISTRIBUTION %
    
    beta_1=mvnrnd(beta_T(:,t-1)',eye(2)*0.01); % the variance of proposal is important
        
    p1=0;
    p2=0;
    for i=1:N
        p1=p1+Y(i)*log(normcdf(X(i,:)*beta_1'))+(1-Y(i))*log(1-normcdf(X(i,:)*beta_1'));
        p2=p2+Y(i)*log(normcdf(X(i,:)*beta_T(:,t-1)))+(1-Y(i))*log(1-normcdf(X(i,:)*beta_T(:,t-1)));
    end
    
    if log(rand(1))<=(p1-p2)
        beta_T(:,t)=beta_1;
        acc=acc+1;
    else
        beta_T(:,t)=beta_T(:,t-1);
    end
    acc_rate(t)=acc/t;
end

toc;

figure;
subplot(2,2,1)
plot(beta_T(1,burn:t));figure(gcf);
title('\beta_1')
subplot(2,2,2)
plot(beta_T(2,burn:t));figure(gcf);
title('\beta_2')
subplot(2,2,3)
histogram(beta_T(1,burn:t));figure(gcf);
title('\beta_1')
subplot(2,2,4)
histogram(beta_T(2,burn:t));figure(gcf);
title('\beta_2')

disp(mean(beta_T(1,burn:t)));
disp(mean(beta_T(2,burn:t)));

disp(std(beta_T(1,burn:t)));
disp(std(beta_T(2,burn:t)));
%}

