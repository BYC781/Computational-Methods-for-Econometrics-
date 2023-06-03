%***********************************%
% monte carlo simulation for SC-SAR %
% Data Generating Process           %  
%***********************************%

clear;
tic;
ssd = 12345;
rng(ssd);


%******* DGP Paremeters **********%
gamma=[-1.5; 0.5; 1];
beta=[0.5; 0.5];
Sigma=[1 0; 0 1];
sigma_a=0.05;
lambda=0.05;
delta=0.3;
mu=[0 0];

% choose group size and number of groups %
N=30;       
G=50;
R=5;


C=cell(R,1);
W=cell(R,1);
X=cell(R,1);
Y=cell(R,1);

for r=1:R
    
    z=zeros(N,G);
    e=zeros(N,G);
    c_1=rand([N,G]);
    c1=rand([N,N,G]);
    for g=1:G
        Mat=mvnrnd(mu,Sigma,N);
        for i=1:N
            for j=1:N
                if c_1(i,g)>=0.7 && c_1(j,g)>=0.7
                    c1(i,j,g)=1;
                elseif c_1(i,g)<=0.3 && c_1(j,g)<=0.3
                    c1(i,j,g)=1;
                else
                    c1(i,j,g)=0;
                end                
                if i==j
                    c1(i,j,g)=0;                    
                end               
            end
            z(i,g)=Mat(i,1);
            e(i,g)=Mat(i,2);
        end       
    end
    C{r}=c1;    
    p=zeros(N,N,G);
    w=zeros(N,N,G);
    for g=1:G
        for i=1:N
            for j=1:N                   
                psi=gamma(1)+gamma(2)*c1(i,j,g)-gamma(3)*abs(z(i,g)-z(j,g));
                p(i,j,g)=exp(psi)/(1+exp(psi));      % Logit setting
                if rand(1)<=p(i,j,g)
                    w(i,j,g)=1;
                else
                    w(i,j,g)=0;
                end
                if i==j
                    w(i,j,g)=0;
                end
            end
        end
    end
    W{r}=w;
    x=randn(N,G)*2;
    alpha=mean(x,1)*delta+randn(1,G)*sqrt(sigma_a);
    y=zeros(N,G);
    for g=1:G
        S=eye(N)-lambda*w(:,:,g);
        y(:,g)=S\([x(:,g) w(:,:,g)*x(:,g)]*beta+ones([N,1])*alpha(g)+e(:,g));
    end
    X{r}=x;
    Y{r}=y;
end

save DGP_I C W X Y;

toc;

