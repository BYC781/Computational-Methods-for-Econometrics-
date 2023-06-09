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
    data{i} = readtable(filename);
    filename2 = sprintf('/network%d.dat', i);
    W{i} = readtable(filename2);
    X{i} = data{i}(:, [1:6 8:14]);
    Y{i} = data{i}(:, 18);
    N(i) = size(W{i}, 1);
end
N = N';

%% set init param
init_lambda = 0.029;
init_a_g = repmat(4.0, 1, 76);
init_b = [-0.193, -0.111, -0.096, -0.012, 0.032, -0.117, -0.145, ...
          0.145, -0.081, 0.125, -0.042, -0.004, 0.029, ...
          -0.005, -0.0029, -0.029, 0.040, -0.044, -0.076, -0.04, ...
          0.005, 0.038, 0.028, -0.007, 0.018, 0.046];
init_sigma2 = 0.459;
init_param = [init_lambda, init_a_g, init_b, init_sigma2]';

%% MLE

options = optimset('fmincon');
lb = [-0.1; -Inf(103, 1)];
ub = [0.1; Inf(103, 1)];
options = optimset(options, 'Display', 'final','MaxIter',5e06,'MaxFunEvals',5e06,'TolFun', 1e-6, 'GradObj', 'off','Hessian', 'off');

%%
[est,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(param) -log_likelihood(param, Y, X, W, N), ...
    init_param,[], [], [], [], lb, ub, [], options);
%% coefficient and standard error
lambda_coef = est(1);
a_g = est(2:77);
b1 = est(78:90);
b2 = est(91:103);
sigma2 = est(104);

se = sqrt(diag(inv(hessian)));

lambda_se = se(1);
a_g_se = se(2:77);
b1_se = se(78:90);
b2_se = se(91:103);
sigma2_se = se(104);

%% output
fprintf('coefficient of lambda=%.4f, standard error of lambda = %.4f\n', lambda_coef,  lambda_se );
fmt=['coefficient of alpha_g =' repmat(' %.4f',1,numel(a_g))];
fprintf(fmt,a_g);
fprintf('\n');
fmt=['standard error of alpha_g =' repmat(' %.4f',1,numel(a_g_se))];
fprintf(fmt,a_g_se);
fprintf('\n');
fmt=['coefficient of beta1 =' repmat(' %.4f',1,numel(b1))];
fprintf(fmt,b1);
fprintf('\n');
fmt=['standard error of beta1 =' repmat(' %.4f',1,numel(b1_se))];
fprintf(fmt,b1_se);
fprintf('\n');
fmt=['coefficient of beta2 =' repmat(' %.4f',1,numel(b2))];
fprintf(fmt,b2);
fprintf('\n');
fmt=['standard error of beta2 =' repmat(' %.4f',1,numel(b2_se))];
fprintf(fmt,b2_se);
fprintf('\n');
fprintf('coefficient of sigma2=%.4f, standard error of sigma2 = %.4f\n', sigma2,  sigma2_se );

%%
save('mle_result.mat','a_g','a_g_se', "b1", "b1_se", "b2", "b2_se", "lambda", "lambda_se", "sigma2", "sigma2_se")
