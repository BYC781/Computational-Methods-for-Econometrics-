R = 100;
N = 400;
X1 = randn(N*R, 1);
X1 = reshape(X1, N, R);
X2 = chi2rnd(1, N*R, 1);
X2 = reshape(X2, N, R);
U1 = gumbelrnd(0, 1, N*R, 1);
U1 = reshape(U1, N, R);
U2 = gumbelrnd(0, 1, N*R, 1);
U2 = reshape(U2, N, R);
Y = zeros(N, R);
for i = 1:R
    Y(:, i) = (X1(:, i) + U1(:, i)) > (X2(:, i) + U2(:, i));
end
beta_seq = zeros(R, 2);
for i = 1:R
    gradient = @(beta) [sum(Y(:,i).*X1(:,i) - X1(:,i).*(exp(beta(1)*X1(:,i) - beta(2)*X2(:,i))./(1+exp(beta(1)*X1(:,i) - beta(2)*X2(:,i))))); 
                        sum(Y(:,i).*X2(:,i) - X2(:,i).*(exp(beta(1)*X1(:,i) - beta(2)*X2(:,i))./(1+exp(beta(1)*X1(:,i) - beta(2)*X2(:,i)))))];
    beta = [0, 0]';
    tol = 1e-6;
    maxiter = 100;
    for j = 1:maxiter
        g = gradient(beta);
        BHHH = g * g';
        if max(abs(g)) < tol
            break;
        end
        beta = beta + (BHHH \ g)';
    end
    beta_seq(i, 1) = beta(1);
    beta_seq(i, 2) = beta(2);
end
