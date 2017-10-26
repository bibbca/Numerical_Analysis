% RK4 Algorithm
% This is Runge-Kutta (Order Four) Algorithm
function [t, ws, act] = RK4 (Function, LeftEndpoint, RightEndpoint, StepSize, InitialVal, ActualSoln)
a = LeftEndpoint; % Left Endpoint
b = RightEndpoint; % Right Endpoint
h = StepSize; % Iteration Number
alph = InitialVal; % Initial Condition
f = Function; % Our ODE.
g = ActualSoln; % Actual Solution to the ODE (if we have it);

N = (b-a)/h;
t = a;
w = alph;

T = zeros(N:1);
W = zeros(N:1);

for i = 1:N
   K1 = h*f(t,w);
   K2 = h*f(t+h/2, w+K1/2);
   K3 = h*f(t+h/2, w+K2/2);
   K4 = h*f(t+h, w+K3);
   w = w + (K1+2*K2+2*K3+K4)/6;
   t = a+i*h;
   T(i) = t;
   W(i) = w;
   act(i) = g(t);
end

t = T'; % Row vector -> Column Vector (Transpose)
RK4 = W'; % Row vector -> Column Vector (Transpose)
Actual = act'; % Row vector -> Column Vector (Transpose)
Table = table(t,RK4, Actual);
disp(Table);
ws = W;
end