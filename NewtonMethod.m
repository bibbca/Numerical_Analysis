%Caleb Bibb 
%Math 341 
%H3

%Newton's Method
%f is our function, a is our initial guess, TOL is our tolerance, I is our
%maximum number of iterations. This program finds the zero.

syms x
f = exp(x)+2^(-x)+2*cos(x)-6;
a = 1;
TOL = 10^-4;
I = 1000;

%Where the magic starts
n = 1;
FD = diff(f);
while n < I
    FA = single(subs(f, x, a));
    FAD = single(subs(FD, x, a));
    E = a-(FA/FAD);
    g = abs(a-E);
    if FA == 0
        Ans = ['The answer is ', num2str(E),' After ',num2str(n),' Iterations.'];
        disp(Ans);
        break
    elseif g < TOL
        Ans = ['The answer is ', num2str(E),' After ',num2str(n),' Iterations.'];
        disp(Ans);
        break
    end
    n = n+1;
    a = E;
end
if n==I
Oops = ['Failed after ', num2str(I),' iterations, We failed! :('];
disp(Oops);
end