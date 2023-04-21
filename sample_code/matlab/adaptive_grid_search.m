% This code demonstrates the adaptive grid search method

m1=6; m2=6;
nR = 0;         % number of grid refinements
int_length = 5; % initial interval length

% construct grid

G1 = linspace(0,5,m1);
G2 = linspace(0,5,m2);
G1 = log(3*G1+1);  % Non-equally spaced grid
G2 = log(3*G2+1);  % Non-equally spaced grid

% initial value 
x1_init = -1e+10;
x2_init = -1e+10;
f_max  = -1e+10;
diff   =  1e+10;


while diff>=1e-3
    
    for i = 1:m1
        for j = 1:m2
            f = G1(i)-0.2*G1(i)^2 + G2(j) - 0.3*G2(j)^2;
            if f >= f_max
                f_max = f;
                x1_max = G1(i);
                x2_max = G2(j);
            end
        end
    end
    
    diff = norm([x1_max x2_max]-[x1_init x2_init]);
    
    % update grid
    int_length = int_length/2;
    G1 = linspace(x1_max-int_length/2, x1_max+int_length/2, m1);
    G2 = linspace(x2_max-int_length/2, x2_max+int_length/2, m2);
    
    x1_init=x1_max;
    x2_init=x2_max;
    
    nR=nR+1;
    disp(nR);
    
end
