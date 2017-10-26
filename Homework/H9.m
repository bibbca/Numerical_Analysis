% Caleb Bibb
% Math 341 H9

%% Section 5.3 #1c
% Use Taylor's method of order two to approx the solution for y' = 1+y/t, 1<=t<=2, y(1)=2, with h = .25
% T^(2) (t(i),w(i)) = f(t(i),w(i)) + (h/2)*f'(t(i),w(i))+...+h^(n-1)/n! * f^(n-1)(t(i),w(i))
% So, f'(t(i),y(t)) = d/dt(1+y/t) = -y'/t^2 = -(1+y/t)/t^2
% Then, T^(2)(t(i),w(i)) = f(t(i),w(i)) + (h/2)f'(t(i),w(i)) = 
f = @(t,w) 1+w./t+(.25/2)*(-(1+w./t)/t.^2); % T^(2)
N = 4; % 4 t-values see below
w(1) = 2; % Initial Condition
t = 1:.25:2; % t is in [1,2] in .25 increments (h = .25)
for i = 1:N
 w(i+1) = w(i)+.25*(f(t(i),w(i)));
end

t = t'; % Row vector -> Column Vector (Transpose)
Taylor2 = w'; % Row vector -> Column Vector (Transpose)
T2 = table(t,Taylor2); % Make Table
disp(T2); % Display Table

%% Section 5.3 #3c
% Repeat #1c using Taylor's method of order 4.
% T^(4) = f(t(i),w(i)) + (h/2)f'(t(i),w(i)) + (h/6)f"(t(i),w(i)) + (h/24)f'"(t(i),w(i))
clear all
f = @(t,w) 1+w./t+(.25/2)*(-(1+w./t)/t.^2)+(.25/6)*(2*(1+w./t)/t.^3)+(.25/24)*(-6*(1+w./t)/t.^4); % T^(4)
N = 4; % 4 t-values see below
w(1) = 2; % Initial Condition
t = 1:.25:2; % t is in [1,2] in .25 increments (h = .25)
for i = 1:N
 w(i+1) = w(i)+.25*(f(t(i),w(i)));
end

t = t'; % Row vector -> Column Vector (Transpose)
Taylor4 = w'; % Row vector -> Column Vector (Transpose)
T4 = table(t,Taylor4); % Make Table
disp(T4); % Display Table
%% Section 5.4 #1c
% Use the modified Euler method to approx the solution and compare:
% y' = 1+y/t, 1<=t<=2, y(1)=2, with h = .25; actual solution: y(t) = t*ln(t) + 2
clear all
f = @(t,w) 1.+w./t;
g = @(t) t.*log(t)+2.*t;
N = 4; % 4 t-values see below
w(1) = 2; % Initial Condition
t = 1:.25:2; % t is in [1,2] in .25 increments (h = .25)
h = 0.25;

for i = 1:N
 w(i+1) = w(i)+(h/2)*(f(t(i),w(i))+f(t(i+1),w(i)+h*f(t(i),w(i))));
end

act = g(t); % Actual Values
t = t'; % Row vector -> Column Vector (Transpose)
ModifiedEuler = w'; % Row vector -> Column Vector (Transpose)
Actual = act'; % Row vector -> Column Vector (Transpose)
ME = table(t,ModifiedEuler,Actual); % Make Table
disp(ME); % Display Table
%% Section 5.4 #5c
% Repeat using the Midpoint Method
clear all
f = @(t,w) 1.+w./t;
g = @(t) t.*log(t)+2.*t;
N = 4; % 4 t-values see below
w(1) = 2; % Initial Condition
t = 1:.25:2; % t is in [1,2] in .25 increments (h = .25)
h = 0.25;

for i = 1:N
 w(i+1) = w(i) + h*f(t(i)+h/2,w(i)+h/2*f(t(i),w(i)));
end

act = g(t); % Actual Values
t = t'; % Row vector -> Column Vector (Transpose)
Midpoint = w'; % Row vector -> Column Vector (Transpose)
Actual = act'; % Row vector -> Column Vector (Transpose)
MP = table(t,Midpoint,Actual); % Make Table
disp(MP); % Display Table
%% Section 5.4 #13c
% Repeat Exercise 1 using Runge-Kutta method of order 4:
%RK4 (Function, LeftEndpoint, RightEndpoint, Iterations, InitialVal, ActualSoln)
RK4(@(t,w)1.+w./t,1,2,4,2,@(t)t.*log(t)+2.*t);

%% Section 5.9 #1a
% Use the Runge-Kutta for Systems Algorithm to approx the solutions of the
% following differential equation and compare the result to the actual
% solution.
% u1' = 3u1+2u2-(2t^2+1)exp(2t), u1(0) = 1;
% u2' = 4u1 + u2 + (t2 + 2t?4)exp(2t) , u2(0) = 1; 0 ? t ? 1; h = 0.2;
%actual solutions u1(t) = 1/3*exp(5t) ? 1/3*exp(?t) + exp(2t) and u2(t) =
%1/3*exp(5t) + 2/3 exp(?t) + t^2exp(2t)
RKSystem ({@(t,u1,u2) 3.*u1+2.*u2-(2.*t.^2+1).*exp(2.*t);@(t,u1,u2) 4.*u1+u2+(t.^2+2.*t-4).*exp(2.*t)}, {@(t) (1/3.*exp(5.*t)-1/3.*exp(-t)+exp(2.*t));@(t) (1/3.*exp(5.*t)+2/3.*exp(-t)+t.^2.*exp(2.*t))}, 0, 1, .2, [0;0], 2);

%% Section 5.9 #3a
% Use the Runge-Kutta for Systems Algorithm to approx the solutions of the
% following differential equation and compare the result to the actual
% solution.
% y" -2y' +y = te^t-t, 0<=t<=1, y(0)=0, y'(0)=0, with h = .1;
% Actual Solution: y(t) = 1/6*t^3*exp(t)+2*exp(t)-t-2
%   u1(t) = y(t), u2(t) = y'(t).  Then u1' = y'(t) = u2. u2' = y" =
%   2y'-y+te^t-t
%       This is a 1st order system of 2 DEs with initial y(0)=0, y'(0)=0
clear all
RKSystem({@(t,u1,u2)u2;@(t,u1,u2)2.*u2-u1+t.*exp(t)-t},{@(t)1/6*t.^3.*exp(t)-t.*exp(t)+2*exp(t)-t-2;@(t)t},0,1,.1,[0;0],2)

%% Section 5.9 #9
% This problem is a predator/prey model.  With x1'(t) = k1x1(t)-k2x2(t) for
% the prey and x2'(t) = k3x1(t)x2(t)-k4x2(t).  We want to solve for
% 0<=t<=4.  Assuming Initial prey = 1000 and preditors = 500 and k1 = 3, k2
% = .002, k3 = .0006, and k4 = .5.  Substituting gets us:
%   x1'(t) = 3*x1(t)-.002x2(t), x2'(t) = .0006x1(t)*x2(t)-.5x2(t)
clear all
[t,C,Act1,Act2]=RKSystem({@(t,x1,x2)3.*x1-.002.*x1.*x2;@(t,x1,x2).0006.*x1.*x2-.5.*x2},{@(t)t;@(t)t;},0,4,.1,[1000;500],2);
figure
plot(t,C(:,1),t,C(:,2));
legend('Prey', 'Predator', 'Location', 'NorthEast')
xlabel('Time')
ylabel('Population')
title('Predator/Prey Model')

% We can see that at t=4, the number of prey is approximately 25 and the
% number of predators is approximately 1258.  There are two "stable
% solutions", one at x1'=x2'=0 (When there are no predators or prey), and
% k1-k2*k2(t)=0 and k3x1(t)-k4=0 which can mean that x1(t)=k4/k3 and
% x2(t)=k1/k2. This occurs ar x1(t)=833.333333 and x2(t)=1500.00. (this
% can't ever happen since we can't have a third of an animal).