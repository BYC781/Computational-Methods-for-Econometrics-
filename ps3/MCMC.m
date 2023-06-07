clear;
rng(7414);
addpath( '/Users/bychen/Documents/ECON-7218/ps3/problem_set_3_sample/group/')
addpath( '/Users/bychen/Documents/ECON-7218/ps3/problem_set_3_sample/network/')
addpath( '/Users/bychen/Documents/ECON-7218/ps3/')

%% Read data

% Initialize a cell array to store the tables
data = cell(1, 76);
W = cell(1, 76);
X = cell(1, 76);
Y = cell(1, 76);
N = zeros(1, 76); 

% Loop through each file
for i = 1:76
    % Read the file into a table
    filename = sprintf('/group%d.dat', i);
    data{i} = table2array(readtable(filename));
    filename2 = sprintf('/network%d.dat', i);
    W{i} = table2array(readtable(filename2));
    X{i} = data{i}(:, 1:17);
    Y{i} = data{i}(:, 18);
    N(i) = size(W{i}, 1);
end
N = N';

%% set prior param
b_0 = zeros(34, 1);
B_0 = eye(34)*3;
a_0 = 4; 
A_0 = 0.5;
k_0 = 0.5; 
v_0 = 0.5;
lambda_prime = 0;

%% set MCMC savers
T = 20000;

lambda_T = zeros(1, T);
alpha_T = zeros(76,T);
beta_T = zeros(34,T);
sige_T = zeros(1,T);
%% starting value
lambda_T(1,1) = 0.0397;
beta_T(:,1) = zeros(34,1);
alpha_T(:,1) = zeros(76,1);
sige_T(1,1) = 1;
%% main process

for t = 2:T
        accept=0;
        while accept==0                                       % propose lambda^*
            if t<2
                lambda_1=mvnrnd(lambda_T(t-1),eye(1)*0.1^2);
            else
                lambda_1=mvnrnd(lambda_T(t-1),cov(lambda_T(1:t-1))*2.38^2)*0.95+mvnrnd(lambda_T(t-1),eye(1)*0.1^2)*0.05;
            end
            if lambda_1>-1/10 && lambda_1<1/10
                accept=1;
            end
        end
        pp_l=0;
        
        for g=1:76
            S_1=eye(N(g))-lambda_1*W{g};
            S_2=eye(N(g))-lambda_T(t-1)*W{g};
            ep_1=S_1*Y{g}-[X{g},W{g}*X{g}]*beta_T(:,t-1)-ones(N(g), 1)*alpha_T(g,t-1);
            ep_2=S_2*Y{g}-[X{g},W{g}*X{g}]*beta_T(:,t-1)-ones(N(g), 1)*alpha_T(g,t-1);
            V=sige_T(t-1)*eye(N(g));
            like_1=log(det(S_1))-(1/2)*(ep_1)'/V*(ep_1);
            like_2=log(det(S_2))-(1/2)*(ep_2)'/V*(ep_2);
            pp_l=pp_l+like_1-like_2;
        end

        pp_l=min(exp(pp_l),1);
        if rand(1)<=pp_l         % accept-reject decision 
            lambda_T(t)=lambda_1;
        else
            lambda_T(t)=lambda_T(t-1);
        end


        for g = 1:76
            
end
