function [  ] =DivDiff(fun,  xval)

syms x

for i = 1:numel(x)
    fval(i) = double(subs(fun, x, xval(i)));%getting function values
end
P1 = fval(1);
for i = 2:numel(x)  %First Difference 
    diff1(i-1) = (fval(i)-fval(i-1))/(xval(i)-xval(i-1));
end
P2 = P1+diff1(1)*(x-xval(1));
for i = 3:numel(x) %Second Difference
    diff2(i-2) = (diff1(i-1)-diff1(i-2))/(xval(i)-xval(i-2));
end
P3 = P2+diff2(1)*prod(x-xval(1:3));
for i = 4:numel(x)  %third difference
    diff3(i-3) = (diff2(i-2)-diff2(i-3))/(xval(i)-xval(i-3));
end

P4 = P3+diff3*prod(x-xval(1:4));
P5 = 0;

if numel(xval) == 5 %chechking for a need for a fifth difference
    for i = 5:numel(x)
        diff4(i-4) = (diff3(i-3)-diff3(i-4))/(xval(i)-xval(i-4));
        P5 = P+diff4(1)*prod(x-xval(1:5));
    end
end
P = [P1 P2 P3 P4 P5];
end