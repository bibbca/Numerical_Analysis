function [ p ] = NewtonHornerMethod( C, TOL, p, N )
%This is Newton's Method that uses Synthetic Division to find f(x) and
%f'(x)
%    Horner's Method evaluates P(x) = C(1)x^n +C(2)x^(n-1)+...
%    Evaluates P(x) at x=p; 
%    Made by Caleb Bibb
i = 1;
while i < N

n = length(C);
b(1) = C(1);
    for j=2:n
        b(j) = C(j)+b(j-1)*p;
    end
p1 = b(n); % P(x) at x=p;
D = b(1:n-1);
e(1) = D(1);
    for j=2:n-1
        e(j) = D(j)+e(j-1)*p;
    end
p2 = e(n-1); % P'(x) at x=p;

% Newton's Method
i = 1;
p0 = p;
p = p0 - (p1/p2);
    if abs(p - p0) < TOL
        break
    else
        i = i+1;
    end
end
end

