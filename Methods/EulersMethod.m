% Euler's Method
function [Ans] = EulersMethod (Function, LeftEndpoint, RightEndpoint, Iterations, InitialVal)
% W is our approximate solution at time T.

a = LeftEndpoint; % Left Endpoint
b = RightEndpoint; % Right Endpoint
N = Iterations; % Max # of Iterations
alp = InitialVal; % Initial Condition
f = Function; %function of y and t

h = (b-a)/N;
t = a;
w = alp;

for i = 1:N
    w = w+h*f(t,w);
    t = a + i*h;
    T(i) = t;
    W(i) = w;
end
if N == 4
    T = [T(1);T(2);T(3);T(4)];
    W = [W(1);W(2);W(3);W(4)];
    I = [1;2;3;4];
    Ans = table(I,T,W);
elseif N == 10
    fy = @(t) t.^2.*(exp(t)-exp(1));
    for i = 1:N
       y =  fy(T(i));
       Y(i) = y;
    end
    T = [T(1);T(2);T(3);T(4);T(5);T(6);T(7);T(8);T(9);T(10)];
    W = [W(1);W(2);W(3);W(4);W(5);W(6);W(7);W(8);W(9);W(10)];
    I = [1;2;3;4;5;6;7;8;9;10];
    Y = [Y(1);Y(2);Y(3);Y(4);Y(5);Y(6);Y(7);Y(8);Y(9);Y(10)];
    Ans = table(I,T,W,Y);
elseif N==25
    T = T(25);
    W = W(25);
    Ans = table(T,W);
elseif N==50
    T = T(50);
    W = W(50);
    Ans = table(T,W);
elseif N==100
    T = T(100);
    W = W(100);
    Ans = table(T,W);
else
    disp('Oops, we didnt expect this N value');
end
end
