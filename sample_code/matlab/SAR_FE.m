%*******************************************************
% ESTIMATE SAR MODEL USING BAYESIAN ESTIMATION
% To run this program, one needs the Matlab package 
% jplv7 from James Lesage, which is downloaded from 
% https://www.spatial-econometrics.com/
%*******************************************************

clear;
load DGP_I;

% if matlabpool('size') > 0
%     matlabpool close;
% end;
% matlabpool open 4
% format short;
%

%ssd = 20200322;
%rng(ssd);
addpath('Users/bychen/Documents/ECON-7218/sample_code/matlab')


N=30;
G=50;
T=5000;  % number of iterations during Markov process
R=1;     % number of replication for Monte Carlo experiment

lambda_R=zeros(1,R);
beta_R=zeros(2,R);
sige_R=zeros(1,R);



    
   
    tic;
        
    %% assign parameter in prior distributions %%
    
    beta_0=zeros([1,2]);
    B_0=eye(2)*3;
    ALPHA_0=1;
    sig0=[1 0.5];
    rho_0=2.2;
    eta_0=0.1;
    
    c_1=1; c_2=0.1; c_3=0.1; c_4=1;
    acc_1=0; acc_2=0; acc_3=0; acc_4=zeros(G,1);
    
    acc_rate1=zeros(T,1); % we only need to focus on block 1 %
    acc_rate2=zeros(T,1); % can drop%
    acc_rate3=zeros(T,1); % can drop%
    acc_rate4=zeros(G,T); % can drop%
    
    lambda_T =zeros([1,T]);  % save for lambda
    beta_T   =zeros([2,T]);  % save for beta
    alpha_T  =zeros([G,T]);  % save for alpha
    sige_T   =zeros([1,T]);  % save for Sigma_e^2
    
    %% starting value of draw %%
    sige_T(1)=1;
        
    for t=2:T
        
        accept=0;
        %while accept==0;                                        % propose lambda^*
            if t<2
                lambda_1=mvnrnd(lambda_T(t-1),eye(1)*0.1^2);
            else
                lambda_1=mvnrnd(lambda_T(t-1),cov(lambda_T(1:t-1))*2.38^2)*0.95+mvnrnd(lambda_T(t-1),eye(1)*0.1^2)*0.05;
            end
            %if lambda_1>-1/(min(max(max(sum(w(:,:,:,r),1))),max(max(sum(w(:,:,:,r),2))))) && lambda_1<1/(min(max(max(sum(w(:,:,:,r),1))),max(max(sum(w(:,:,:,r),2)))));
            %    accept=1;
            %end;
        %end;
        
        pp_l=0;
        V=sige_T(t-1)*eye(N);
        
        for g=1:G
            S_1=eye(N)-lambda_1*W{r}(:,:,g);
            S_2=eye(N)-lambda_T(t-1)*W{r}(:,:,g);
            ep_1=S_1*Y{r}(:,g)-[X{r}(:,g),W{r}(:,:,g)*X{r}(:,g)]*beta_T(:,t-1)-ones([N,1])*alpha_T(g,t-1);
            ep_2=S_2*Y{r}(:,g)-[X{r}(:,g),W{r}(:,:,g)*X{r}(:,g)]*beta_T(:,t-1)-ones([N,1])*alpha_T(g,t-1);
            like_1=log(det(S_1))-(1/2)*(ep_1)'/V*(ep_1);
            like_2=log(det(S_2))-(1/2)*(ep_2)'/V*(ep_2);
            pp_l=pp_l+like_1-like_2;
        end
        pp_l=min(exp(pp_l),1);
        if rand(1)<=pp_l         % accept-reject decision 
            lambda_T(t)=lambda_1;
            acc_2=acc_2+1;
        else
            lambda_T(t)=lambda_T(t-1);
        end
        acc_rate2(t,1)=acc_2/t;
        
        
        % THE SAMPLING OF BETA FROM PROSTERIOR DISTRIBUTION %
        ZVY=0;
        ZVX=0;
        for g=1:G
            SS=eye(N)-lambda_T(t)*W{r}(:,:,g); % use new (t), since lambda already updated
            YY=SS*Y{r}(:,g)-ones([N,1])*alpha_T(g,t-1); % use old (t-1)
            ZZ=[X{r}(:,g) W{r}(:,:,g)*X{r}(:,g)];
            ZVX=ZVX+ZZ'/V*ZZ;
            ZVY=ZVY+ZZ'/V*YY;
        end
        beta_T(:,t)=norm_rnd(inv(inv(B_0)+ZVX))+((inv(B_0)+ZVX)\((B_0)\beta_0'+ZVY));
        
        
        % THE SAMPLING OF SIGMA_E^2 FROM PROSTERIOR DISTRIBUTION %
        ep_v=zeros(N*G,1);
        for g=1:G
            SS=eye(N)-lambda_T(t)*W{r}(:,:,g);
            ep=SS*Y{r}(:,g)-[X{r}(:,g) W{r}(:,:,g)*X{r}(:,g)]*beta_T(:,t)-ones([N,1])*alpha_T(g,t-1);
            ep_v((g-1)*N+1:g*N)=ep;
        end
        rho_1=rho_0+length(ep_v);
        sige_T(t)=(ep_v'*ep_v+eta_0)/chis_rnd(1,rho_1);
        
        % THE SAMPLING OF ALPHA_G FROM PROSTERIOR DISTRIBUTION %
        dd=(ALPHA_0^(-1)+(sige_T(t))^(-1)*ones([1,N])/(eye(N))*ones([N,1]))^(-1);
        for g=1:G
            SS=eye(N)-lambda_T(t)*W{r}(:,:,g);
            YY=SS*Y{r}(:,g);
            XX=[X{r}(:,g) W{r}(:,:,g)*X{r}(:,g)];
            alpha_T(g,t)=sige_T(t)^(-1)*dd*ones([1,N])/(eye(N))*(YY-XX*beta_T(:,t))+randn(1)*sqrt(dd);
        end
                
        if (t/100)-round(t/100)==0
            disp(t);
            disp(lambda_T(t));
            disp(beta_T(:,t)');
            disp(sige_T(t));
            
            subplot(3,1,1)
            plot(lambda_T(1:t));figure(gcf);
            title('\lambda')
            subplot(3,1,2)
            plot(beta_T(1,1:t));figure(gcf);
            title('\beta_1')
            subplot(3,1,3)
            plot(beta_T(2,1:t));figure(gcf);
            title('\beta_2')
            drawnow;
        end
                
    end
    time=toc;
    
    disp('lambda');disp(mean(lambda_T(1000:t)));
    disp('beta');disp(mean(beta_T(:,1000:t),2));
    disp('sigma_e^2');disp(mean(sige_T(1000:t)));
    disp('time'); disp(time);
    
    lambda_R(r)=mean(lambda_T(1000:t));
    beta_R(:,r)=mean(beta_T(:,1000:t),2);
    sige_R(r)=mean(sige_T(1000:t));
    
    info.q = 0.025;
    info.r = 0.005;
    info.s = 0.95;
    info.p1 = 0.3;
    info.p2 = 0.5;
    vnames = strvcat('lambda');
    coda(lambda_T(1000:5:t)',vnames,info);
    plot(lambda_T(1000:5:t));
    
    %pause;
     

disp('lambda');disp(mean(lambda_R,2));
disp(std(lambda_R,0,2));
disp('beta');disp(mean(beta_R,2));
disp(std(beta_R,0,2));
disp('sigma_e^2');disp(mean(sige_R,2));
disp(std(sige_R,0,2));



