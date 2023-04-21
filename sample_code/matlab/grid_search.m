% This code demonstrates the grid search method

m1=6; m2=6;

% construct grid

G1 = linspace(0,5,m1);
G2 = linspace(0,5,m2);

G1 = log(3*G1+1);  % Non-equally spaced grid
G2 = log(3*G2+1);  % Non-equally spaced grid

% initial value, to be very small 
x1_max = -1e+10;
x2_max = -1e+10;
f_max  = -1e+10;

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

disp([x1_max,x2_max]);

% surface plot


[X1,X2] = meshgrid(0:0.1:5);
Z = X1-0.2*X1.^2+X2-0.3*X2.^2;
figure
mesh(X1,X2,Z);



