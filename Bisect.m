%Caleb Bibb 
%Math 341 
%H2

%Bisection Method
%f is our function, a and b are our interval bounds, TOL is our tolerance, I is our
%maximum number of iterations. This program finds the zero within the
%bounds

syms x
f = 2*x+3*cos(x)-exp(x);
a = 0;
b = 2;
TOL = 10^-4;
I = 1000;


%start off with the first iteration
n = 1;
while n < I
    g = (b-a)/2;
    m = a+g;
    FM = single(subs(f, x, m));
    FA = single(subs(f, x, a));
    if FM == 0
        Ans = ['The answer is ', num2str(m),' After ',num2str(n),' Iterations.'];
        disp(Ans);
        break
    elseif g < TOL
        Ans = ['The answer is ', num2str(m),' After ',num2str(n),' Iterations.'];
        disp(Ans);
        break
    elseif FA*FM>0
        a = m;
    else
        b = m;
    end
    n = n+1;
end
if n==I
Oops = ['Failed after ', num2str(I),' iterations, We failed! :('];
disp(Oops);
end