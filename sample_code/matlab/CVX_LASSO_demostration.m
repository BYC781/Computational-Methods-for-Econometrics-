clear;
addpath C:\Users\ssunr\Dropbox\Dropbox\teaching_NTU\Econ7218\matlab_codes\cvx
addpath C:\Users\ssunr\Dropbox\teaching_NTU\Econ7218\matlab_codes\cvx\structures
addpath C:\Users\ssunr\Dropbox\teaching_NTU\Econ7218\matlab_codes\cvx\lib
addpath C:\Users\ssunr\Dropbox\teaching_NTU\Econ7218\matlab_codes\cvx\functions
addpath C:\Users\ssunr\Dropbox\teaching_NTU\Econ7218\matlab_codes\cvx\commands
addpath C:\Users\ssunr\Dropbox\teaching_NTU\Econ7218\matlab_codes\cvx\builtins

% set up cvx if needed % 
cvx_setup;

format short; 
rng('default');
rng(1);


%% Data generating process 

N = 400;
K = 250;
X = [ones(N,1) randn(N,K-1)];
B = 1+randn(K,1);
B(rand(K,1)<0.7) = 0; % set 70% of coefficients to zero

Y = X*B + randn(N,1);
lambda = 0.1;


%% run fminsearch (or fminunc)

fun = @(B_c)(1/N)*(Y - X*B_c)'*(Y - X*B_c) + lambda*norm(B_c,1);

init=zeros(K,1);

% options = optimset('fminsearch');
% options = optimset(options, 'LargeScale', 'on', 'Display','final',...
%        'FunValCheck','off','MaxIter',5e06,'MaxFunEvals',5e06,'Tolfun', 1e-6, 'GradObj','off', 'Hessian', 'off');
% tic;
% [est,~,~] = fminsearch(fun,init,options);
% toc;

options = optimset('fminunc');
options = optimset(options, 'Display','final',...
        'FunValCheck','off','MaxIter',5e06,'MaxFunEvals',5e06,'Tolfun', 1e-6, 'GradObj','off', 'Hessian', 'off');
tic;
[est,fval] = fminunc(fun,init,options);
toc;


figure;
subplot(1,2,1)
hist(B,20);
subplot(1,2,2)
hist(est,20);

%% run CVX

%cvx_solver mosek
cvxp = cvx_precision( 'high' );
tic;
cvx_begin
    variable B_c(K)
    minimize((1/N)*sum_square(Y - X*B_c) + lambda*norm(B_c,1))
cvx_end
toc

figure;
subplot(1,2,1)
hist(B,20);
subplot(1,2,2)
hist(B_c,20);



