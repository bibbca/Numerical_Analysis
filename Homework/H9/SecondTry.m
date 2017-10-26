
a = 0; % Left Endpoint
b = 1; % Right Endpoint
N = 5; % Iteration Number
alpha = [1;1]; % Vector of Initial Conditions
m = 2; % Number of Equations

Z = {@(t,u1,u2) 3.*u1+2.*u2-(2.*t.^2+1).*exp(2.*t);@(t,u1,u2) 4.*u1+u2+(t.^2+2.*t-4).*exp(2.*t)};
G = {@(t) (1/3.*exp(5.*t)-1/3.*exp(-t)+exp(2.*t));@(t) (1/3.*exp(5.*t)+2/3.*exp(-t)+t.^2.*exp(2.*t))};

h = (b-a)/N;
t = a;

for j = 1:m
    w(j) = alpha(j);
    A = w;
end

for i = 1:N
   for j = 1:m
       k(1,j) = h*Z{j}(t,w(1),w(2));
   end
   for j = 1:m
      k(2,j) = h*Z{j}(t+h/2,w(1)+1/2*k(1,1),w(2)+1/2*k(1,2));
   end
   for j = 1:m
      k(3,j) = h*Z{j}(t+h/2,w(1)+1/2*k(2,1),w(2)+1/2*k(2,2));
   end
   for j = 1:m
      k(4,j) = h*Z{j}(t+h,w(1)+k(3,1),w(2)+k(3,2));
   end
   for j= 1:m
      w(j) = w(j) + (k(1,j)+2*k(2,j)+2*k(3,j)+k(4,j))/6; 
   end
   C(i+1,:) = w;
   t = a+i*h;
end
t = a:h:b;
t = t';
Act1 = G{1}(t);
Act2 = G{2}(t);
C(1,1) = alpha(1);
C(1,2) = alpha(2);

Ans = table(t,C(:,1),Act1,C(:,2),Act2);
Ans.Properties.VariableNames = {'t' 'w1' 'Actual1' 'w2' 'Actual2'};
disp(Ans);
