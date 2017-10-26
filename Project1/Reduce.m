function [ C ] = Reduce( C, p )
%Reduce takes your p and C and puts out a vector of coefficients one degree
%fewer.
% Made by Caleb Bibb
%   C is a vector of coefficients, p is where your 0 is.
n = length(C);
o(1) = C(1);
    for j=2:n
        o(j) = C(j)+o(j-1)*p;
    end
C = o(1:n-1);
end

