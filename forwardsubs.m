function [x]=forwardsubs(L,b)

% Forward substitution method for lower triangular matrices L

[n,m] = size(L);
x = [];
if n ~= m 
    warning('L is not a square matrix');
    return;
end
if min(abs(diag(L)))==0
    error('Singular matrix');
    return;
end
x(1)=b(1)/L(1,1);

for i=2:n
    x(i)=(b(i)-L(i,1:i-1)*x(1:i-1)')/L(i,i);
end

return