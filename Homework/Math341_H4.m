%Caleb Bibb
%MATH 341 H4

%% Section 2.6 #4a
%Fund approximations to within 10^(-5) to all the zeros of each of the
%following polynomials using Muller's method.
clear all
syms x
f=x^4+5*x^3-9*x^2-85*x-136; % Our functions
p0 = -4; % Initial guess 1
p1 = -1; % Initial guess 2
p2 = 4; % Initial guess 3
TOL = 10^(-5); % Our Tolerance
N = 1000; % Maximum Number of iterations

i = 3;
h1 = p1-p0;
h2 = p2-p1;
Delta1 = (double(subs(f, x, p1))-double(subs(f, x, p0)))/h1;
Delta2 = (double(subs(f, x, p2))-double(subs(f, x, p1)))/h2;
a = (Delta2-Delta1)/(h2-h1);

for k = 1:4
    while i<N
    b = Delta2-h2*a;
    D =(b^2-4*a*(double(subs(f, x, p2))))^(1/2);
    if abs(b-D)<abs(b+D)
        E = b+D;
    else
        E = b-D;
    end 
    h = (double(subs(f, x, p2))*(-2))/E;
    p = p2+h;
    if abs(h)<TOL
        disp(['One of the roots is ' num2str(p) ' in #4a'])
        break 
    end
    p0 = p1;
    p1 = p2;
    p2 = p;
    h1 = p1-p0;
    h2 = p2-p1;
    Delta1 = (double(subs(f, x, p1))-double(subs(f, x, p0)))/h1;
    Delta2 = (double(subs(f, x, p2))-double(subs(f, x, p1)))/h2;
    a = (Delta2-Delta1)/(h2-h1);
    i = i+1;
    end
    div = p;
    f = simplify(f/(x-div));
    k = k+1;
end 
%% Section 2.6 #5a
% Use Newton's Method to find, within 10^-3 the zeros and critical points
% of the following function. f(x) = x^3-9*x^2+12
clear all
syms x
f = x^3-9*x^2+12; % Our function
a = 1; % Our initial guess
TOL = 10^-3; % Our Tolerance
I = 1000; % Maximum number of iterations
n = 1; % Current Iteration Number
h = f;
for k = 1:3
    dh = diff(h);
    while n < I
    b = a - double(subs(h, x, a))/double(subs(dh, x, a));
    if abs(a-b)<TOL
        Ans = ['For 5a, There is a zero at ', num2str(b)];
        disp(Ans)
        break
    end
    a = b;
    i = i+1;
    end
    h = simplify(h/(x-b)); 
    k = k+1;
end
f = x^3-9*x^2+12; % Our function
j = diff(f);
for k = 1:2
    dh = diff(j);
    while n < I
        b = a - single(subs(j, x, a))/single(subs(dh, x, a));
        if abs(a-b)<TOL
            Ans = ['For 5a, There is a critical point at ', num2str(b)];
            disp(Ans)
            break
        end
        a = b;
        i = i+1;
    end
    j = simplify(j/(x-b));
    k = k+1;
end

%% Section 3.1 #1a
% f(x) = cos(x), x0 = 0, x1 = 0.6, x2 = 0.9
% approx f(0.45) using interpolation polynomials of of degree 1 and 2 and
% find the absolute error

% f(x) = cos(x):    at  x=0,    x=1,    x=2
% yi = f(xi):       is  f(x0)   f(x1)   f(x2)
%        =>         cos(0)     cos(.6)  cos(.9)
clear all
syms x
f = cos(x);
a = double(subs(f, x, 0))
b = double(subs(f, x, .6))
c = double(subs(f, x, .9))
% This makes our three points: (0,1), (.6, .8253), (.9, .6216)
% We need to find L1,0(x) and L1,1(x) to find P1(x)
L10 = (x-.6)/(0-.6)
L11 = (x-0)/(.6-0)
P1 = a*L10 + b*L11
% Absolute Error = |f(x) -P1(x)|
AbsoluteError = abs(double(subs(f, x, .45))-double(subs(P1, x, .45)))

%% Section 3.1 #3a
% Use Theorem 3.3 to find an error bound for #1a
%
%
%% Section 3.1 #5a
% Use Lagrange interpolating polynomials of degrees one, two, and three to approximate
% f(8.4) if f(8.1) = 16.94410, f(8.3) = 17.56492, f(8.6) = 18.50515, f(8.7) = 18.82091
%
%% Section 3.1 #17
% Suppose you need to construct eight-decimal-place tables for the common, or base-10, logarithm function from x=1 to x=10 in such a way that linear interpolation is accurate to within 10^(-6).  Determine a bound for the step size for this table. What choice of step size would you make to ensure that x=10 is included in the table?
%
%
%% Programing Problem - Dived Difference Coefficients
% function [D, P] = divDiffTable(xv,y)
% xv and y are vectors where f(xv) = y.

%Output: D = matrix containing divided-difference coefficients in
%its lower triangle.
% syms x;
% n = length(y);
% if length(xv)~=n 
%     error('xv and y are not the same length'); 
% end
% 
% D = zeros(n,n);
% P = sym(zeros(n,n));
% D(:,1) = y(:); % First column is zeroth order difference, f[x_i] = y_i
% for j=2:n
%     for i=j:n
%         D(i,j) = (D(i,j-1)-D(i-1,j-1))/(xv(i)-xv(i-j+1));
%     end
% end
% for i=1:n
%     if i == 1
%         P(1)=D(1,1);
%     else
%     P(i) = D(i,i)*prod(x-xv(1:i))+P(i-1);
%     end
% end
%% Section 3.3 #1a
% Use Eq. (3.10) or Algorithm 3.2 to construct interpolating polynomials of degree one, two, and three for f(8.4) if f(8.1) = 16.94410, f(8.3) = 17.56492, f(8.6) = 18.50515, f(8.7) = 18.82091
% Time to take out our new program for a spin.
clear all
[D,P] = divDiffTable([8.1 8.3 8.6 8.7], [16.94410 17.56492 18.50515 18.82091])
% We can see that P1 = 16.9440999999999998+3.10409999(x-8.1)
clear all
syms x
P0 = 16.9440999999999988
P1 = P0+3.10409999999999320*(x-8.1);
P2 = P1+0.0600000000000333600*(x-8.1)*(x-8.3);
P3 = P2 +-0.00208333333334478552*(x - 8.1000000)*(x - 8.300)*(x - 8.6000)

P184 = double(subs(P1, x, 8.4));
P284 = double(subs(P2, x, 8.4));
P384 = double(subs(P3, x, 8.4));
disp(['P1(8.4) = ', P184,' P2(8.4) = ', P284,' P3(8.4) = ', P384]);
%% Section 3.3 #7 (a and b)
% Use Algorithm 3.2 to construct the interpolating polynomial of degree three for the unequally spaced points given in the following table: f(-.1)=5.30000, f(0)=2.00000, f(.2)=3.19000, f(.3)=1.00000
% b) Add f(.35)=.97260 to the table, and construct the interpolating polynomial of degree four.
[D,P] = divDiffTable([-.1 0 .2 .3],[5.30000 2.00000 3.19000 1.00000])
% We can see that P3(x) = 5.2999999+(-32.9999999)*(x
% +0.1000)+129.83333333*(x + 0.1000)*x-556.666666*(x+.1000)*x*(x-.20000)
[D,P] = divDiffTable([-.1 0 .2 .3 .35],[5.30000 2.00000 3.19000 1.00000 .97260])
% We can see that P4(x) = P3(x)+
% 2730.24338624338679(x+.100)(x)(x-.200)(x-.300)
%% Section 3.3 #11
% Show that the cubic polynomials both interpolate the data

x=-10:0.5:10;
P= 3-2.*(x+1)+0.*(x+1).*(x)+(x+1).*(x).*(x-1);
Q= -1+4.*(x+2)-3.*(x+2).*(x+1)+(x+2).*(x+1).*(x);
figure(3);plot(x,P,x,Q,'--');grid on;
title('Showing that both polynomials interpolate the same data');
%% Section 3.3 #17
% For a function f, the forward-divided differences are given by the table.
%
% What we're given: x2 = .7, f[x2] = 6, f[x1,x2] = 10, f[x0, x1, x2] =
% 50/7, x1 = .4, x0 = 0.
%
%
%
%
%
%
