%Caleb Bibb 
%Math 341 
%H3

%Secant Method
%f is our function, a and b are our initial guesses, TOL is our tolerance, I is our
%maximum number of iterations. This program finds the zero.

syms x
f = exp(x)+2^(-x)+2*cos(x)-6;
a = 1;
b = 2;
TOL = 10^-4;
I = 1000;
FA = single(subs(f, x, a));
FB = single(subs(f, x, b));
%Where the magic starts
n = 2;
while n < I
    
    E = b-FB*(b-a)/(FB-FA);
    g = abs(a-E);
    if FA == 0
        Ans = ['The answer is ', num2str(E),' After ',num2str(n),' Iterations.'];
        disp(Ans);
        break
    elseif FB == 0
        Ans = ['The answer is ', num2str(E),' After ',num2str(n),' Iterations.'];
        disp(Ans);
        break
    elseif E == 0
        Ans = ['The answer is ', num2str(E),' After ',num2str(n),' Iterations.'];
        disp(Ans);
        break
    elseif g < TOL
        Ans = ['The answer is ', num2str(E),' After ',num2str(n),' Iterations.'];
        disp(Ans);
        break
    end
    n = n+1;
    a = b;
    FA = FB;
    b = E;
    FB = single(subs(f, x, b));
end
if n==I
Oops = ['Failed after ', num2str(I),' iterations, We failed! :('];
disp(Oops);
end