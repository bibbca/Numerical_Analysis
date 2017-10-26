% Caleb Bibb
% Math 341
% H8

%% Section 4.6 # 1d
% Compute the Simpson's rule approx S(a,b), S(a, (a+b)/2), and S((a+b)/2, b) for Int{0, pi/4} x^2*sin(x) dx
clear all
% S(a,b) = h/3[f(a) + 4f(a+h) +f(b) = S1
% S(a, (a+b)/2) = h/6*[f(a)+4f(a+h/2)+f(a+h)] = S2
% S((a+b)/2, b) = h/6*[f(a+h)+4f(a+3h/2)+f(b)] = S3
f = @(x) x.^2*sin(x);
a = 0;
b = pi/4;
h = (b-a)/2;

FA = f(a);
FAH = f(a + h);
FB = f(b);

S1 = h/3 * (FA+4*FAH+FB)
S2 = h/6 * (FA+4*f(a+h/2)+FAH)
S3 = h/6 * (FAH+4*f(a+3*h/2)+FAH)
%% Section 4.6 #2d
% Use adaptive quadrature to find approxs within 10^-3 for #1d.  Don't use an algorithm to do it.
clear all
f = @(x) x.^2.*sin(x);
a = 0;
b = pi/4;
h = (b-a)/2;

FA = f(a)
FAH = f(a + h)
FB = f(b)

S11 = h/3 * (FA+4*FAH+FB);
S12 = h/6 * (FA+4*f(a+h/2)+FAH);
S13 = h/6 * (FAH+4*f(a+3*h/2)+FAH);

y = h/6*(FA+4*f(a+h/2)+2*FAH+4*f(a+3*h/2)+FB);
act = integral(f, a, b);
err = act - y;
disp(['The estimation is ',num2str(y),' the error is ',num2str(err),'< 10^-3']);
%% Section 4.6 #3a
% Use Adaptive curvature to approx Int{1, 3} exp(2x)*sin(3x)
clear all
f = @(x) exp(2.*x).*sin(3.*x);
a = 0;
b = pi/4;
h = (b-a)/2;

FA = f(a)
FAH = f(a + h)
FB = f(b)

S11 = h/3 * (FA+4*FAH+FB);
S12 = h/6 * (FA+4*f(a+h/2)+FAH);
S13 = h/6 * (FAH+4*f(a+3*h/2)+FAH);

y = h/6*(FA+4*f(a+h/2)+2*FAH+4*f(a+3*h/2)+FB);


act = integral(f, a, b);
err = act - y;
disp(['The estimation is ',num2str(y),' the error is ',num2str(err),'< 10^-3']);
%% Section 5.1 #1a
% Use Theorem 5.4 to show that the following IVP y' = y*cos(t), 0<=t<=1, y(0)=1
%       Please see attached page.

%% Section 5.1 #2a
% Show that the following IVP has a unique solution and find the solution.  Can Thm 5.4 be used in each case?
%       Please see attached page.
%% Section 5.1 #3a
% For f(t,y) = t^(2)y +1, does f satisfy a Lipschitz condition on D={(t,y)|0<=t<=1, -inf<y<Inf}?
% Can Thm 5.6 be used to show that the IVP y' = f(t,y), 0<=t<=1, y(0)=1 is well-posed?
%       Please see attached page.
%% Section 5.2 #1c
% Use Euler's method to approx solutions for y' = 1+y/t, 2<=t<=3, y(1)=2, h = .25
EulersMethod(@(t,y) 1+y/t,1,2,4,2)
%% Section 5.2 #9
% Given the IVP y' = (2/t)y+t^(2)e^t, 1<=t<=2, y(1)=0, with exact solution y(t) = t^2(e^t - e^1):
% a) Use Euler's method with h=.1 to approx the solution and compare it with the actual values of y
EulersMethod(@(t,y) 2/t.*y+t.^2.*exp(t),1,2,10,0)
% b) Use the answers generated in part a) and linear interpolation to approx the following values of y, and compare them to the actual values: y(1.04), y(1.55), y(1.97)
y100 = 0;
y110 = .27183;
y104 =  y100 + 4/10*y110;
err104 = abs(y104 - .119987497);

y150 = 3.1874;
y160 = 4.6208;
y155 = y150 + 5/10*y160;
err155 = abs(y155 - 4.7886350208);

y190 = 11.748;
y200 = 15.398;
y197 = y190 + 7/10*y200;
err197 = abs(y197 - 17.2792984355576628);

T = [1.04;1.55;1.97];
Approx = [y104;y155;y197];
Error = [err104;err155;err197];
Table = table(T,Approx,Error);
disp(Table);
% c) Compute the value of h necessary for |y(t_i)-w_i| <= .1, Using Eqn 5.10
        % |y(ti)-w(i)|<=(h*M)/(2L)*[e^(L(t(i)-a))-1]
% We can see that f(t,y) is continuous and satisfies a Lipschitz condition
% such that D = {(t,y)|a<=t<=b and -inf<y<inf}.  We can see that
% |y"(t)|<=3+sqrt(3)=M    =>|y(ti)-w(i)|<=(h*3+sqrt(3))/(2*L)*[e^(L(t(i)-a))-1].  If we find
% the Lipchitz Condition we see that the partial derivative with respect to
% y is 2/t which is maximized at 2 over our interval so that is our
% Lipchitz constant.
% => .1>=(h*(3+sqrt(3))/(2*2)*[e^(2(t(i)-1))-1]
% => .1>=(h*(3+sqrt(3))/4*[e^(2t(i)-2)-1] =>
L = 2;
C = ((.1*2*L)/(exp(.97)-1))/(sqrt(3)+3);
disp(C);
%% Section 5.2 #11
% Given the IVP y'=-y+t+1, 0<=t<=5, y(0)=1, with exact solution y(t) = exp(-t) +t:
% a) Approximate y(5) using Euler's method with h =.2, h = .1, and h=.05.
EulersMethod(@(t,y) -y+t+1,0,5,25,1)
EulersMethod(@(t,y) -y+t+1,0,5,50,1)
EulersMethod(@(t,y) -y+t+1,0,5,100,1)
% b) Determind the optimal value of h to use in computing y(5), assuming delta = 10^-6 and that Eqn 5.14 is valid
%   Since the exact solution is y(t) = exp(-t)+t, we can show that y"(t) =
%   exp(-t).  Therefore, |y"(t)|<=1=M.
%   h = sqrt(2*delta/M)
G = sqrt(2*10^-6/1);
disp(G);