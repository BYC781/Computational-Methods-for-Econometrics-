addpath('/Users/bychen/Documents/ECON-7218/sample code/compecon2011/CEtools');
addpath('/Users/bychen/Documents/ECON-7218/sample code/compecon2011/CEdemos');

%% Computes multivariate trapezoid rule quadrature nodes and weights

n=4;
a=-1;
b=1;

[x,w] = qnwtrap(n,a,b);

%% Computes multivariate Simpson quadrature nodes and weights

n=5;
a=-1;
b=1;

[x,w] = qnwsimp(n,a,b);


%% Computes multivariate Guass-Legendre quadrature nodes and weights
n=2;
a=-1;
b=1;
[x,w] = qnwlege(n,a,b);

n=3;
a=-1;
b=1;
[x,w] = qnwlege(n,a,b);


