function [x,L,U] = meg(A,b)

% Gauss Elimination Method

[n,m] = size(A);
x = [];
if n ~= m 
    warning('A is not a square matrix');
    return;
end
Ak = A;
bk = b;
for k=1:n-1
    M(k,k+1:n) = Ak(k+1:n,k)/Ak(k,k);
    Ak(k+1:n,k+1:n)=Ak(k+1:n,k+1:n) - M(k+1:n,k)*Ak(k,k+1:n);
    bk(k+1:n)=bk(k+1:n)-M(k+1:n,k)*bk(k);
end
