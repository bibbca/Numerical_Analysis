%Test of summation method

f = @(x) x.^2+23;
a = 1;
h = .5;
n = 5;

for k = 1:2^(n-2)
    z(k) = f(a+(k-.5)*h);
end
disp(z);
g = sum(z);
disp(g);