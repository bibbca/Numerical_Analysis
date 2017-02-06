function [Y] =DividedDiff(funct, xvals)
syms x

for i = 1:numel(xvals) % This iteration finds f(x)
    val(i) = double(subs(funct, x, xvals));
end
D1 = val(1);

for i = 2:numel(xvals) % Finds 1st Divided Difference
    divdiff1(i-1) = (val(i)-val(i-1))/(xvals(i)-xvals(i-1));
end
D2 = (x-xvals(1))*divdiff1(1)+D1;

for i = 3:numel(xvals) %Finds 2nd Divided Difference
    divdiff2(i-2) = (divdiff1(i)-divdiff1(i-2))/(xvals(i)-xvals(i-2));
end
D3 = prod(x-xvals(1:3))*divdiff2(1)+D2;

for i = 4:numel(xvals) % 3rd Divided Difference
    divdiff3(i-3) = ((divdiff2(i)-divdiff2(i-3))/(xvals(i)-xvals(i-3)));
end
D4 = prod(x-xvals(1:4))*divdiff3(1)+D3;

for i = 5:numel(xvals) % 4th Divided Difference
    divdiff4(i-4) = ((divdiff3(i)-divdiff3(i-4))/(xvals(i)-xvals(i-4)));
end
D5 = prod(x-xvals(1:5))*divdiff4(1)+D4;

for i = 6:numel(xvals) % 5th Divided Difference
    divdiff5(i-5) = ((divdiff4(i)-divdiff4(i-5))/(xvals(i)-xvals(i-5)));
end
D6 = prod(x-xvals(1:5))*divdiff5(1)+D5;

Y = [D1 D2 D3 D4 D5 D6]; %Throws it all together
end