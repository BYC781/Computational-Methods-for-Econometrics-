%% obj fn
function ll = log_likelihood(param, Y, X, W, N)
    lambda = param(1);
    a_g = param(2:77);
    b1 = param(78:90);
    b2 = param(91:103);
    sigma2 = param(104);

    sigma2term = zeros(76,1);
    e_gT_e_g = zeros(76,1);
    logdetterm = zeros(76,1);
    log2pi = zeros(76,1);
    for i = 1:76
        log2pi(i) = -log(2*pi)*(N(i) / 2);
        S_g = eye(size(W{i}, 1)) - lambda * table2array(W{i});
        Y_g = table2array(Y{i});
        X_g = table2array(X{i});
        W_g = table2array(W{i});
        l_g = ones(N(i), 1);
        e_g = S_g * Y_g - X_g * b1 - W_g * X_g * b2 - a_g(i) * l_g;
        e_gT_e_g(i) = e_g' * e_g;
        sigma2term(i) = - log(sigma2) * (N(i) / 2);
        logdetterm(i) = log(det(S_g));
    end    

   
    ll = sum(log2pi) + sum(sigma2term) + sum(logdetterm) +  -(2 * sigma2)^-1 * sum(e_gT_e_g);
end