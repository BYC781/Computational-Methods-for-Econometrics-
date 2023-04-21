%% This code compares the computation speed of different matrix decomposition methods 

n = 5000;
A = rand(n, n);
A = A*A.';
B = randn(n,1);

%% L-U decomposition 
disp('L-U');
tic;
dA = decomposition(A,'lu');
x1= dA\B;
toc;

%% QR decomposition
disp('Q-R');
tic;
dA = decomposition(A,'qr');
x2= dA\B;
toc;

%% Cholesky decomposition
disp('Cholesky');
tic;
dA = decomposition(A,'chol');
x3= dA\B;
toc;


%% MATLAB backslash
disp('MATLAB backslash');
tic;
x4= A\B;
toc;
