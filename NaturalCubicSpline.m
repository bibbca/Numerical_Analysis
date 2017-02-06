% Natural Cubic Spline

xvals = [1 2 3 4 5]; % xvals as vector
yvals = [3 6 9 12 15]; % yvals as vector

if length(xvals) ~= length(yvals)
    error('xvals and yvals are not the same length');
end

n = length(xvals)-1; %Sets n = # of points - 1
for i = 0:n
