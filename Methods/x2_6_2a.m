clear all
syms x
f = x^4+5*x^3-9*x^2-85*x-136; %Our function
a = 1; %Initial Guess
TOL = 10^-5; %Tolerance 
I = 1000;

%Where the magic starts
n = 1;
FD = diff(f);
FD2 = diff(FD);
while n < I
    FA = single(subs(f, x, a));
    FAD = single(subs(FD, x, a));
    FAD2 = single(subs(FD2, x, a));
    E = a-((FA*FAD)/(FAD^2-FA*FAD2));
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
Oops = ['Modified Newtons Failed after ', num2str(I),' iterations, We failed! :('];
disp(Oops);
end
f1 = f/(x-E);
