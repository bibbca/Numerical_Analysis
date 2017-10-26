% Caleb Bibb
% Numerical H5

%% Section 5.6 #1a
% Use Adams-Bashforth 4 method to approx the solutions to the
% following IVP and compare with actual values
%   y' = te^(3t)-2y, 0<=t<=1, y(0) = 0, with h = .2; 
%   actual solution y(t) = (1/5)*te^(3t) - (1/25)*e^(3t) + (1/25)*e^(-2t)
f = @(t,w) t.*exp(3.*t)-2.*w;
a = 0;
b = 1;
h = .2;
alpha = 0;
N = (b-a)/h;
actual = @(t) (1/5).*t.*exp(3.*t) - (1/25).*exp(3.*t) + (1/25).*exp(-2.*t);
%Runge Kutta Method
 [t, ws, act]=RK4(f, a, .2*3, .2, alpha, actual);
t = [a t'];
w = [alpha ws];
act = [alpha act];
%AB4 Explicit Method
for i = 4:N
     t(i+1) = t(i)+h;
     w(i+1) = w(i)+h*(55*f(t(i),w(i))-59*f(t(i-1),w(i-1))+37*f(t(i-2),w(i-2))-9*f(t(i-3),w(i-3)))/24;
end
for i = 1:length(t)
   act(i) = actual(t(i)); 
end
t = t';
AB4 = w';
Actual = act';
TableAB4 = table(t,AB4,Actual);
disp(TableAB4);
%% Section 5.6 #4a
% Use Adams-Moulton 3 method to approx the solutions to the
% following IVP and compare with actual values. Use exact starting values
% and explicitly solve for w(i+1).
%   y' = te^(3t)-2y, 0<=t<=1, y(0) = 0, with h = .2; 
%   actual solution y(t) = (1/5)*te^(3t) - (1/25)*e^(3t) + (1/25)*e^(-2t)
f = @(t,w) t.*exp(3.*t)-2.*w;
a = 0;
b = 1;
h = .2;
alpha = 0;
N = (b-a)/h;
actual = @(t) (1/5).*t.*exp(3.*t) - (1/25).*exp(3.*t) + (1/25).*exp(-2.*t);
%Runge Kutta Method
 [Table, t, ws, act]=RK4(f, a, b, .2, alpha, actual);
t = [a t'];
w = [alpha ws];
act = [alpha act];
%AM3 Implicit Method
for i = 3:N
     t(i+1) = t(i)+h;
     w(i+1) = w(i)+h*(9*f(t(i+1),w(i+1))+19*f(t(i),w(i))-5*f(t(i-1),w(i-1))+f(t(i-2),w(i-2)))/24;
end
for i = 1:length(t)
   act(i) = actual(t(i)); 
end
t = t';
AB3 = w';
Actual = act';
TableAB3 = table(t,AB3,Actual);
disp(TableAB3);
%% Section 5.6 #5a
% Use Algorithm 5.4 to approx:
%   y' = te^(3t)-2y, 0<=t<=1, y(0) = 0, with h = .2; 
% Initial Values
f = @(t,w) t.*exp(3.*t)-2.*w;
a = 0;
b = 1;
h = .2;
N = (b-a)/h;
t(1) = a;
w(1) = 0;
% Magic Ensues
for i = 2:4
   K(1) = h*f(t(i-1),w(i-1)); 
   K(2) = h*f(t(i-1)+h/2,w(i-1)+K(1)/2);
   K(3) = h*f(t(i-1)+h/2,w(i-1)+K(2)/2);
   K(4) = h*f(t(i-1)+h,w(i-1)+K(3));
   w(i) = w(i-1)+(K(1)+2*K(2)+2*K(3)+K(4))/6;
   t(i) = a+(i)*h;
   A(i) = t(i);
   B(i) = w(i);
end
for i = 5:N+1
    t(i) = a+i*h;
    w(i) = w(4)+h*(55*f(t(4),w(4))-59*f(t(3),w(3))+37*f(t(2),w(2))-9*f(t(1),w(1)))/24;
    w(i) = w(4)+h*(9*f(t(i),w(i))+19*f(t(4),w(4))-5*f(t(3),w(3))+f(t(2),w(2)))/24;
    A(i) = t(i);
    B(i) = w(i);
    for j= 1:3
        t(j) = t(j+1);
        w(j) = w(j+1);
    end
    t(3) = t(i);
    w(3) = w(i);
end
A = A'; % Row vector -> Column Vector (Transpose)
Adams4 = B'; % Row vector -> Column Vector (Transpose)
AB4 = table(A,Adams4); % Make Table
disp(AB4); % Display Table
