function [APP] = AdaptiveQuadrature (Function, LeftEndpoint, RightEndpoint, Tolerance, LevelLimit);
%This function calculates the Adaptive Quadrature
a = LeftEndpoint; % Left Endpoint
b = RightEndpoint; % Right Endpoint
TOL = Tolerance; % Tolerance
N = LevelLimit; % Limit N to number of levels
f = Function; % Our Function

APP = 0;
i = 1;
TOL(i) = 10*TOL;
a(i) = a;
h(i) = (b-a)/2;
FA(i) = f(a);
FC(i) = f(a+h(i));
FB(i) = f(b);
S(i) = h(i)*(FA(i)+4*FC(i)+FB(i))/3;
L(i) = 1;

    while i > 0
    FD = f(a(i)+h(i)/2);
    FE = f(a(i)+3*h(i)/2);
    S1 = h(i)*(FA(i)+4*FD+FC(i))/6;
    S2 = h(i)*(FC(i)+4*FE+FB(i))/6;
    v(1) = a(i);
    v(2) = FA(i);
    v(3) = FC(i);
    v(4) = FB(i);
    v(5) = h(i);
    v(6) = TOL(i);
    v(7) = S(i);
    v(8) = L(i);
    i = i-1;
    if abs(S1+S2-v(7))<v(6)
        APP = APP +(S1+S2); 
    else
        if v(8) >= N
            disp('Level Exceeded :(');
            return;
        else
            i = i+1;
            a(i) = v(1)+v(5);
            FA(i) = v(3);
        end
    end
    end
   disp(APP);
end