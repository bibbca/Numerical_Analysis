% Adaptive Quadrature Plot
function adaptive_plot (LeftEndpoint, RightEndpoint, Xvals, Function, plotChoice)
% This function plots our Anonymous Functions and plots points at each step
%  Of course this is only if desired, written by Caleb Bibb.
 
f = Function;      % Original Function
a = LeftEndpoint;  % Original Left Endpoint
b = RightEndpoint; % Original Right Endpoint
x = Xvals;         % Stepsizes from adaptive_quad

if plotChoice == 'Y'
    y = f(x);           % Maps Adaptive Quadrature Points to an appropriate codomain
    xinterval = [a b];  % Sets interval of the domain
    figure
    fplot(f,xinterval); % Plots Function (now the easiest way to plot a function)
    hold on
    scatter(x,y);       % Adds in Adaptive Quadrature Points

    % Sets Axis when applicable and appropriate
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    hold off
end
end