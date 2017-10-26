% Caleb Bibb
% Programing Project

% Collecting Data / Seperating smart users from dumb
% These increase user experience by allowing them to input values
C = input('Got some coefficients for me? They must be in square brackets. '); % We need at least 3
TOL = input('What is our tolerance? ');
N = input('How many iterations should we do in our Newtons Method? ');
% Would have used apostraphe but... problems.
n = length(C);
if n<3
    disp('Our vector needs to be at least 3 elements');
end
Ans = zeros(n-1,1);
for v = 1:n-3
p = MullerMethod(C);
p = NewtonHornerMethod(C,TOL,p,N);
C = Reduce(C,p);
Ans(v) = p;
end
a = C(1);
b = C(2);
c = C(3);
Ans(v+1) = 2*c/(-b-sqrt(b^2-4*a*c));
Ans(v+2) = 2*c/(-b+sqrt(b^2-4*a*c));
disp(Ans);