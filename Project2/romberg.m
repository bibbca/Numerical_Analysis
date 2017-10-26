% Create a function called romberg, that accepts the function, the interval,
%and the tolerance and uses Romberg integration to return the approx value
%and the step size required. 
% Written By Caleb Bibb
function [h_last, Approx] = romberg (LeftEndpoint, RightEndpoint,  Tolerance, Function)
% Test Vars
a = LeftEndpoint; % Left Endpoint
b = RightEndpoint; % Right Endpoint
f = Function; % Anonymous Function
TOL = Tolerance;

% Placeholders
Check1 = 1;
Check2 = 1;
k = 1;
N = 1;
i = 0;
% Actual Algorithm
while i == 0
        h(k) = (b-a)/2^(k-1);
    
    if k>1 %Loop deals with k = 1 (cause g=1:2^(-1) is a problem)
        for g = 1:2^(k-2)
            Sums(g) = f(a+(2*g-1)*h(k));
        end
        SUM = sum(Sums);
        R(k,1) = 1/2*(R(k-1,1)+h(k-1)*SUM);
    else
        R(k,1) = (b-a)/2*(f(a)+f(b));
    end

    for j = 2:k
        R(k,j) = R(k,j-1)+(R(k,j-1)-R(k-1,j-1))/(4^(j-1)-1);
    end
    if k>1
        Check1 = abs(R(k,k)- R(k-1,k-1));
    end
    if k>2
        Check2 = abs(R(k-1,k-1)- R(k-2,k-2));
    end
    if N == 1000
        disp('Romberg Iteration Limit Reached. Romberg Aborted.')
    end
    if Check1 < TOL && Check2 < TOL
        i = 1;
    end
    k = k+1;
    N = N+1;
end
h_last = h(end);
Approx = R(end);
end
