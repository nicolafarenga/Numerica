function [x]=backwardsubs(U,b)

% Backward substitution method for upper triangular matrices U

[n,m] = size(U);
x = [];
if n ~= m 
    warning('U is not a square matrix');
    return;
end
if min(abs(diag(U)))==0
    error('Singular matrix');
    return;
end
x(n)=b(n)/L(n,n);

for i=n-1:1
    x(i)=(b(i)-U(i,i+1:n)*x(1:i-1))/U(i,i);
end

end