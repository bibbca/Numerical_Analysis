% Caleb Bibb
% Math 341 H7

%% Section 4.3 #1a
% Approximate Int{0.5,1} of x^4 dx using the Trapezoidal rule.
% Trapezoidal Rule: Int{a,b} f(x)dx = h/2[f(x0)+f(x1)] - h^3/12*f"(xi)
clear all
syms x
f = x^4;
x0 = .5;
x1 = 1;
h = x1 - x0;

trap = h/2 * (double(subs(f, x, x0)) + double(subs(f, x, x1)));
disp(trap);

%% Section 4.3 #3a
% Find a bound for the error in #1a using the error formula, and compare
% this to the actual error.
clear all
syms x
f = x^4;
x0 = .5;
x1 = 1;
h = x1 - x0;

% x^4 is maximized at 1 so xi = 1
xi = 1;
df = diff(f);
dff = diff(df);
err = (h^(3)/12)*double(subs(dff,x,xi));
act = double(int(f,x,x0,x1));
trap = h/2 * (double(subs(f, x, x0)) + double(subs(f, x, x1)));
ActualError = double(act - trap);
Ans = ['Our Method Error is ', num2str(err), ' and our actual error is ' num2str(ActualError)];
disp(Ans);
%% Section 4.3 #5a
% Approximate Int{0.5,1} of x^4 dx using Simpson's rule.
% Simpson's Rule: Int{a,b} f(x)dx = h/3[f(x0)+4f(x1)+f(x2)] - h^5/90*f""(xi)
clear all
syms x
f = x^4;
a = .5;
b = 1;
n = 2;

h = (b-a)/n;
xi0 = double(subs(f,x,a)) + double(subs(f,x,b));
xi1 = 0;
xi2 = 0;
for i = 1:n-1
   X = a+i*h;
   if mod(X,2) == 0
       %Even
       xi2 = xi2 +double(subs(f,x,X));
   else
       %Odd
       xi1 = xi1+double(subs(f,x,X));
   end
   xi = h*(xi0+2*xi2+4*xi1)/3;
end
Ans = ['Using Composite Simpsons Rule we get ',num2str(xi)];
disp(Ans);
%% Section 4.3 #7a
% Find a bound for the error in #5a using the error formula, and compare
% this to the actual error.
clear all
syms x
f = x^4;
x0 = .5;
x1 = .75;
x2 = 1;
h = (x2 - x0)/2;

% x^4 is maximized at 1 so xi = 1
xi = 1;
df = diff(f);
dff = diff(df);
df3 = diff(dff);
df4 = diff(df3);
err = (h^(5)/90)*double(subs(df4,x,xi));
act = integral(@(x) x.^4,x0,x1);
simps = h/3 * (double(subs(f, x, x0))+ 4*double(subs(f, x, x1)) + double(subs(f, x, x2)));
ActualError = double(act - simps);
Ans = ['Our Method Error is ', num2str(err), ' and our actual error is ' num2str(ActualError)];
disp(Ans);
%% Section 4.3 #15
% Find the degree of precision of the quadrature formula Int{-1,1} f(x)dx =
% f(-sqrt(3)/3) + f(sqrt(3)/3)

%    Int{-1,1} f(x)dx = f(-sqrt(3)/3) + f(sqrt(3)/3)
% => f(x) = x
% Int(f(x)dx){-sqrt(3)/3,sqrt(3)/3} = 0;
% How about f(x) = x^2
% Int(f(x)dx){-sqrt(3)/3,sqrt(3)/3} = 0
% So we can see that our degree of precision is 1.

%% Section 4.4 #1a
clear all
syms x
f = x*log(x)
a = 1;
b = 2;
n = 4;

h = (b-a)/n;
g = (double(subs(f,x,a))+double(subs(f,x,b)))*h/2;
for i = 1:n-1
    g = h*double(subs(f,x,(a+i*h)))+g;
end
disp(['The Composite Trapezoidal rule gives ', num2str(g)]);
%% Section 4.4 #3a
g = h/3 * (double(subs(f,x,a))+double(subs(f,x,b)));
for i = 1:n-1
   if mod(i,2)==0
       %even
       g = h/3*2*double(subs(f,x,(a+i*h)))+g;
   else
       %odd
       g = h/3*4*double(subs(f,x,(a+i*h)))+g;
   end
end
disp(['The Composite Simpsons rule gives ', num2str(g)]);
%% Section 4.4 #11a
% Determine the values of n and h required to approx Int{0,2}
% exp(2x)*sin(3x) dx to within 10^-4.  Use Composite Trapezoidal Rule.

% Recall that the error term for Composite Trapezoidal Rule is given by:
% E<=|(b-a)/(12)*h^2*f"(mu)|
% We want to find h and n so that E<= 10^-4. So let's compute some
% derivates:
clear all
syms x
f = exp(2*x)*sin(3*x);
df = diff(f);
df2 = diff(df,x);
b = 2;
a = 0;
% We can see that f" is maximized at x=2 within x=[0,2].
% This yields a value:
yield = double(subs(df2,x,2));
% Now, if E<=|(b-a)/(12)*h^2*f"(mu)| = 2/12*h^2*f"(mu)<=(1/6)*705.3601*h^2
%       <=117.5600*h^2
% So if 117.56h^2 = 10^-4 =>  8.5063*10^-7=h^2  =>
h = sqrt( (10^-4) / (yield*(1/6)));

% We know that (b-a)/h = n => 2/h = n
n = ceil((b - a)/h);
Ans = ['h = ', num2str(h), ', n = ', num2str(n)];
disp(Ans);
%% Section 4.4 #11b
% Determine the values of n and h required to approx Int{0,2}
% exp(2x)*sin(3x) dx to within 10^-4.  Use Composite Simpson's Rule.

% Recall that the error term for Composite Trapezoidal Rule is given by:
% E<=|(b-a)/(180)*h^4*f""(mu)|
% We want to find h and n so that E<= 10^-4. So let's compute some
% derivates:
clear all
syms x
f = exp(2*x)*sin(3*x);
df = diff(f);
df2 = diff(df);
df3 = diff(df2);
df4 = diff(df3);
b = 2;
a = 0;

range = linspace(0,2);
for i=1:length(range)
   d4f(i) = double(subs(df4, x, range(i))); 
end
figure(1);plot(range,d4f)
% We can see that f"" is maximized at x=pi/2 within x=[0,2].
% This yields a value:
x1 = pi/2;
yield = double(subs(df4,x,x1));
Test = ['f""(', num2str(x1), ') = ', num2str(yield)];
disp(Test);
% Now, if E<=|(b-a)/(180)*h^4*f"(mu)| =
% 2/180*h^4*f'(mu)<=(1/90)*2844.7844*h^4
%       <=31.6087*h^4
% So if 31.6087h^4 = 10^-4 =>  3.1637*10^-6=h^4  =>
h =( (10^-4) / (yield*(1/90)))^(1/4);

% We know that (b-a)/h = n => 2/h = n
n = ceil((b - a)/h);
Ans = ['We should use h = ', num2str(h), ', and n = ', num2str(n),'.'];
disp(Ans);


%% Section 4.4 #22
% A car laps a race track in 84 seconds. The speed of the car at each
% 6-second interval is determined by using a radar gun and is given from
% the beginning of the lap, in feet/second, by the entries in the following
% table: (0,124), (6,134), (12,148), (18,156), (24,147), (30,133),
% (36,121), (42,109), (48,99), (54,85), (60,78), (66,89), (72,104),
% (78,116), (84,123)
% How long is the track?
% We can see that the estimated distance is the Int{0,84} of s(t) dt where
% s = speed and t = time.  We have 14 intervals, so h = (84-0)/14 = 6.
% Simpson's Rule (with n=14):
h = (84-0)/14;
yay = (h/3)*(124+(4*134)+(2*148)+(4*156)+(2*147)+(4*133)+(2*121)+(4*109)+(2*99)+(4*85)+(2*78)+(4*89)+(2*104)+(4*116)+123);
Ans = ['We can see that the track is ', num2str(yay), ' feet long using Simpsons Rule'];
disp(Ans);

%% Section 4.5 #1b
% Use equation 4.34 to calculate R{k,1} for k >=2
% Use Romberg integration to compute R_{3,3} for the following intergral:
% Int{0,1} x^2*exp(-x) dx
clear all
syms x k
f = x^2*exp(-x); % Function
a = 0; % Left Endpoint
b = 1; % Right Endpoint

h(1) = (b-a)/1;
R(1,1) = h(1)/2*(double(subs(f,x,a))+double(subs(f,x,b)));
h(2) = (b-a)/2;
R(2,1) = (1/2)*(R(1,1)+h(1)*double(subs(f,x,(a+h(2)))));
R(2,2) = R(2,1)+(1/(4^(1)-1))*(R(2,1)-R(1,1));
h(3) = (b-a)/4;
R(3,1) = (1/2)*(R(2,1)+h(2)*double(subs(f,x,(a+h(3))))+double(subs(f,x,(a+3*h(3)))));
R(3,2) = R(3,1)+(1/(4^(1)-1))*(R(3,1)-R(2,1));
R(3,3) = R(3,2)+(1/(4^(2)-1))*(R(3,2)-R(2,2));
disp(['We can see that R(3,3) = ',num2str(R(3,3))]);
%% Section 4.5 #3b
%Use equation 4.34 to calculate R{k,1} for k >=2
clear all
syms x
f = x^2*exp(-x); % Function
a = 0; % Left Endpoint
b = 1; % Right Endpoint

h(1) = (b-a)/1;
R(1,1) = h(1)/2*(double(subs(f,x,a))+double(subs(f,x,b)));
h(2) = (b-a)/2;
R(2,1) = (1/2)*(R(1,1)+h(1)*double(subs(f,x,(a+h(2)))));
R(2,2) = R(2,1)+(1/(4^(1)-1))*(R(2,1)-R(1,1));
h(3) = (b-a)/4;
R(3,1) = (1/2)*(R(2,1)+h(2)*double(subs(f,x,(a+h(3))))+double(subs(f,x,(a+3*h(3)))));
R(3,2) = R(3,1)+(1/(4^(1)-1))*(R(3,1)-R(2,1));
R(3,3) = R(3,2)+(1/(4^(2)-1))*(R(3,2)-R(2,2));
h(4) = (b-a)/8;
V = [.00366958 .0291456 .0714468 .123581];
R(4,1) = (1/2)*(R(3,1)+h(4)*sum(V));
R(4,2) = R(4,1)+(1/(4^(1)-1))*(R(4,1)-R(3,1));
R(4,3) = R(4,2)+(1/(4^(2)-1))*(R(4,2)-R(3,2));
R(4,4) = R(4,3)+(1/(4^(3)-1))*(R(4,3)-R(3,3));
disp(['We can see that R(4,4) = ',num2str(R(4,4))]);
%% Section 4.5 #5b
%Use equation 4.34 to calculate R{k,1} for k >=2
% Use Romberg integration to approx the integrals in Exercise 1 to within
% 10^-6. Compute the Romberg table until eith |R(n-1,n-1)-R(n,n)|<10^-6 or
% n = 10.  Compare your results to the exact value of the integral.
clear all

f = @(x) x.^2.*exp(-x); % Function
a = 0; % Left Endpoint
b = 1; % Right Endpoint
act = integral(f,a,b);

h(1) = (b-a)/1;
R(1,1) = h(1)/2*(f(a)+f(b));
h(2) = (b-a)/2;
R(2,1) = (1/2)*(R(1,1)+h(1)*f(a+h(2)));
R(2,2) = R(2,1)+(1/(4^(1)-1))*(R(2,1)-R(1,1));
Check(1) = abs(R(2,2) - R(1,1));
h(3) = (b-a)/4;
R(3,1) = (1/2)*(R(2,1)+h(2)*f(a+h(3))+f(a+3*h(3)));
R(3,2) = R(3,1)+(1/(4^(1)-1))*(R(3,1)-R(2,1));
R(3,3) = R(3,2)+(1/(4^(2)-1))*(R(3,2)-R(2,2));
Check(2) = abs(R(3,3)-R(2,2));
h(4) = (b-a)/8;
R(4,1) = (1/2)*(R(3,1)+h(3)*f(a+h(4))+f(a+3*h(4))+f(a+5*h(4))+f(a+7*h(4)));
R(4,2) = R(4,1)+(1/(4^(1)-1))*(R(4,1)-R(3,1));
R(4,3) = R(4,2)+(1/(4^(2)-1))*(R(4,2)-R(3,2));
R(4,4) = R(4,3)+(1/(4^(3)-1))*(R(4,3)-R(3,3));
Check(3) = abs(R(4,4)-R(3,3));
h(5) = (b-a)/16;
R(5,1) = (1/2)*(R(4,1)+h(4)*f(a+h(5)));
R(5,2) = R(5,1)+(1/(4^(1)-1))*(R(5,1)-R(4,1));
R(5,3) = R(5,2)+(1/(4^(2)-1))*(R(5,2)-R(4,2));
R(5,4) = R(5,3)+(1/(4^(3)-1))*(R(5,3)-R(4,3));
R(5,5) = R(5,4)+(1/(4^(4)-1))*(R(5,4)-R(4,4));
Check(4) = abs(R(5,5)-R(4,4));
disp(Check)
error = abs(R(5,5)-act);
disp(['We can see that the Rhomberg method has completed it within 5 iterations.  Our actual error is ', num2str(error),'.']);