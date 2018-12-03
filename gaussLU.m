function [x,L,U] = gaussLU(A,b)

% Gauss Elimination Method as a LU factorization of A

[n,m] = size(A);
x = [];

if n ~= m 
    warning('A is not a square matrix');
    return;
end

if abs(A(1:n,1:n))<(sum(A(1:n,:))-A(1:n,1:n))
    warning('A is a singular matrix');
    return;
end

if abs(A(1:n,1:n))==0
    warning('null pivot');
    return;
end

M = eye(n,n);
Ak = A;
bk = b;

for k=1:n-1
    M(k+1:n,k) = Ak(k+1:n,k)/Ak(k,k);
    Ak(k+1:n,k+1:n) = Ak(k+1:n,k+1:n) - M(k+1:n,k)*Ak(k,k+1:n);
end

U = triu(Ak);
L = M;

y = forwardsubs(L,b);
x = backwardsubs(U,y);
return