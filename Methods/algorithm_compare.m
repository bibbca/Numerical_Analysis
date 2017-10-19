%This code uses the bisection method, fixed point method, and Newton's
%method to solve the same problem.
function algorithm_compare(f,fprime,g,a_end,b_end,N,varargin) %no outputs.
%inputs are the functions f, f', g, left endpoint, right endpoint, number
%of iterations, and an optional initial guess for the fixed point algorithm
%and Newton's method.  f, f', and g must be function handles or anonymous
%functions

%Bisection method
bis=zeros(N+1,1); %vector to store iterates from the bisection method
a=a_end; %left endpoint to start bisection method
b=b_end; %right endpoint for bisection
fa=f(a);
for i=1:N+1
    bis(i)=a+(b-a)/2;
    fp=f(bis(i));
    if fa*fp>0
        a=bis(i);
    else
        b=bis(i);
    end
end

%Fixed point method
fixed=zeros(N+1,1); %vector to store iterates from fixed point method
if nargin==7 %check to see if optional initial guess was entered
    fixed(1)=varargin{1}; %optional initial guess. Note varargin is a cell array
else
    rand_guess=a_end+(b_end-a_end)*rand; %random number between a and b
    fixed(1)=rand_guess; 
end
for i=1:N
    fixed(i+1)=g(fixed(i));
end

%Newton's Method
newt=zeros(N+1,1);
if nargin==7
    newt(1)=varargin{1};
else
    newt(1)=rand_guess;
end
for i=1:N
    newt(i+1)=newt(i)-f(newt(i))/fprime(newt(i));
end

%Secant method
sec=zeros(N+1,1);
if nargin==7
    sec(1)=varargin{1};
else
    sec(1)=rand_guess;
end
sec(2)=newt(2); %steal our second value from the results of Newton's method
q0=f(sec(1));
q1=f(sec(2));
for i=2:N
    sec(i+1)=sec(i)-q1*(sec(i)-sec(i-1))/(q1-q0);
    if isnan(sec(i+1))
        break
    end
    q0=q1;
    q1=f(sec(i+1));
end

%the following code displays the outputs from the iteration in the command
%window.  See the matlab documentation for the function fprint f, along
%with the information about formatSpec
fprintf('\n\n%10s %20s %20s %20s %20s\n\n','N','Bisection','Fixed Point','Newton','Secant')
for i=1:N+1
    fprintf('%10d %20.15f %20.15f %20.15f %20.15f\n',i-1,bis(i),fixed(i),newt(i),sec(i))
end

end
