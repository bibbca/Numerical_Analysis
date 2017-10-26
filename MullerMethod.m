function [ p ] = MullerMethod( C )
%This function runs one iteration of Muller's Method.
%   input C which is a vector of coefficients, output p which is a one
%   iteration close guess.
%   Made By Caleb Bibb
syms x
% Arbitrary Initial Guesses, originally the first 3 coefficients, but this
% caused a division by 0 for some functions.
p0 = 2;
p1 = 38;
p2 = -12;
f = poly2sym(C);

h1 = p1 - p0;
h2 = p2 - p1;
delta1 = (double(subs(f,x,p1))-double(subs(f,x,p0)))/h1;
delta2 = (double(subs(f,x,p1))-double(subs(f,x,p0)))/h1;
d = (delta2 - delta1)/(h2 + h1);

b = delta2 + h2*d;
D = (b^2-4*double(subs(f,x,p2))*d)^(1/2);

if abs(b - D) < abs(b + D)
    E = b + D;
else
    E = b - D;
end

h = -2*double(subs(f,x,p2))/E;
p = p2 + h;

end

