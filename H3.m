% %Caleb Bibb
% %Math 341
% %H3
% 
%Section 2.3 #6a
% %Use Newton's Method to find solutions accurate to within 10^-4 for
% %exp(x)+2^(-x)+2*cos(x)-6=0 for [1,2]
clear all
syms x
f = exp(x)+2^(-x)+2*cos(x)-6;
a = 1;
TOL = 10^-5;
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

%% Section 2.3 #8a
%Repeat #6a using the Secant Method (I'm using initial guesses .50 and .75)
clear all
syms x
f = exp(x)+2^(-x)+2*cos(x)-6;
a = 1;
b = 2;
TOL = 10^-5;
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

%% Section 2.3 #15
%The following describes Newton’s method graphically: Suppose that f'(x) exists on [a, b] 
%and that f'(x)!= 0 on [a,b]. Further, suppose there exists one p ? [a, b] such that f(p) = 0,
%and let p0 ? [a, b] be arbitrary. Let p1 be the point at which the tangent line to f at
%p0, f(p0)) crosses the x-axis. For each n ? 1, let pn be the x-intercept of the line tangent 
%to f at (pn?1, f ( pn?1)). Derive the formula describing this method.
%
% Newton's method makes tangent lines.  To do this method you take p0 and
% find where f(p0) crosses the x-axis.  This means that your slope for p_n
% is f'(p_(n+1)) for any p_n.  If we reindex this we get the slope equals
% f'(p_n) for any p_(n-1). But what will the x-coordinate be when
% f'(p_n)=0? Well, that is of course p_(n-1).  Based on this we can derive
% y=mx+b => y(p_n)=f'(p_n)*x+p_(n-1).
%% Section 2.3 #19
%The iteration equation for the Secant method can be written in the simpler
%form. Explain why, in general, this iteration equation is likely to be
%less accurate than the one given in Algorithm 2.4

%As p_(n-1) and p_(n-2) get closer and closer together (nearly equal) we
%can see that the numerator and the denominator both have subtraction of
%these nearly equal numbers.  This leads to reduced accuracy.
%% Section 2.4 #1a
% %Use Newton's method to find solutions accurate to within 10^(-5) for x^2-2*xe^(-x)+e^(-2*x) = 0; for [0, 1];
clear all
syms x
f = x^2-2*x*exp(-x)+exp(-2*x);
a = .75;
TOL = 10^-5;
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
% 
 %% Section 2.4 #3a
% % Repeat #1a using modified Newton's Method (Eq. 2.13). Is there an improvement in speed and accuracy over exercise 1?
clear all
syms x
f = x^2-2*x*exp(-x)+exp(-2*x);
a = 1;
TOL = 10^-5;
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
Oops = ['Failed after ', num2str(I),' iterations, We failed! :('];
disp(Oops);
end
% %In fact we do see a improvement in speed.

%% Section 2.4 #11
%Show that the Bisection Algorithm 2.1 gives a sequence with an error bound
%that converges linearly to 0.

%We know that the error bound for the bisection method is (b-a)/2^n by Theorem 2.1.So,
%((b-a)/(2^(n+1)))/((b-a)/(2^n))=1/2 after some canceling. We can see that this error bound converges only linearly. 
%% Section 2.6 #2a
%Remember to use Newton's method combined with Horner's method and
%deflation
%Find approximations to within 10^?5 to all the zeros of each of the following polynomials 
%by first finding the real zeros using Newton’s method and then reducing to polynomials of lower degree to determine 
%any complex zeros.
%f(x)=x^4+5x^3-9x^2-85x-136




