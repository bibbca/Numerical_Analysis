function [I, h_final, x] = adaptive_quad (LeftEndpoint, RightEndpoint,  Tolerance, Function, plotChoice, RomSteps)
% This function calculates the integral I = Integral from [a,b] of f(x) dx.
% Adapted by Caleb Bibb from Algorithm 4.3
% Inputs
a = LeftEndpoint; % Left Endpoint
b = RightEndpoint; % Right Endpoint
TOL = Tolerance; % Tolerance
N = 1000; % Number of Levels (Arbitrary Large Positive Integer)
f = Function; % Anonymous function

% Start

APP = 0;
i = 1;

TOL(i) = 10*TOL;
a(i) = a; % Left Endpoint of our interval
h(i) = (b-a)/2; % Stepsize
xvals1 = [a a+h(i) b];

FA(i) = f(a); % function at a.
FC(i) = f(a+h(i)); % function at a+h.
FB(i) = f(b); % function at b.

S(i) = h(i)*(FA(i)+4*FC(i)+FB(i))/3; 
%Approx from Simpson's Method for entire interval
L(i) = 1;

while i>0
    FD = f(a(i) + h(i)/2);
    FE = f(a(i) + 3/2*h(i));
    % Approxs from Simpson's Method for halves of subintervals.
    S1 = h(i)*(FA(i)+4*FD+FC(i))/6;
    S2 = h(i)*(FC(i)+4*FE+FB(i))/6;
    
    %Save data at this level:
    v(1) = a(i);
    v(2) = FA(i);
    v(3) = FC(i);
    v(4) = FB(i);
    v(5) = h(i);
    v(6) = TOL(i);
    v(7) = S(i);
    v(8) = L(i);
    
    i = i-1; % Delete the level
    
    if abs(S1 + S2 -v(7))< v(6)
        APP = APP + (S1+S2);
    else
        if v(8) >=N
            disp('Level Exceeded, Procedure Fails');
            break;
        else
            % Data for Right Half Subinterval
            i = i+1;
            a(i) = v(1)+v(5);
            
            FA(i) = v(3);
            FC(i) = FE;
            FB(i) = v(4);
            
            h(i) = v(5)/2;
            TOL(i) = v(6)/2;
            S(i) = S2;
            L(i) = v(8)+1;
            
            % Data for Left Half Subinterval
            i = i+1;
            a(i) = v(1);
            
            FA(i) = v(2);
            FC(i) = FD;
            FB(i) = v(3);
            h(i) = h(i-1);
            xvals = [xvals1, a+h(i)/2, a+3*h(i)/2];
            TOL(i) = TOL(i-1);
            S(i) = S1;
            L(i) = L(i-1);
        end
    end
end
%Outputs:
x = sort(xvals);% Let's output our Stepsizes to the next function
h_final = h(end); %h_final is the step size required.
I = APP; % APP approxs the Actual Value to within Tolerance 
adaptive_plot (LeftEndpoint, RightEndpoint, x, Function, plotChoice);
end