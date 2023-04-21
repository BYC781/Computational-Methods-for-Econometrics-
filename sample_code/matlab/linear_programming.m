% This code solve the LP problem 
%   min  -2x-y
%   s.t. 4x +y <= 400
%         x +y <= 300 
%        2x+5y <= 1000 
%            x >= 0 
%            y >= 0

A=[4 1 
  1 1
  2 5
  -1 0
  0 -1];

B=[400 300 1000 0 0];

f=[-2 -1];

x=linprog(f,A,B);