%Caleb Bibb
%Math 341
%H4

%% Programming Problem 1
%Implement Divided Difference Coefficients
% 
% x = [0 1 2 3];
% y = [0 1 2 3];
% 
% n = length(y); 
% if length(x)~= n, error('x and y are not paired up'); end
% 
% %magic
% T = zeros(n,n);
% for a = 2:n
%     for b = a:n
%         T(b,a) = (T(b,a-1)-T(b-1,a-1))/(x(b)-x(b-a+1))
%     end
% end
% disp(T(b,a))
function D = divDiffTable(x,y)
% divdiffTable Construct a table of divided-difference coefficients
%
% Synopsis: D = divDiffTable(x,y)
%
% Input: x,y = vectors containing the tabulated y = f(x) data
%
% Output: D = matrix containing divided-difference coefficients in
% its lower triangle. Diagonal entries are coefficients of the Newton polynomial
n = length(y);
if length(x)~=n, error('x and y are not compatible'); end

D = zeros(n,n);
D(:,1) = y(:); % First column is zeroth order difference, f[x_i] = y_i
for j=2:n
    for i=j:n
        D(i,j) = (D(i,j-1)-D(i-1,j-1))/(x(i)-x(i-j+1));
    end
end