function [H] = cholesky(A)

% Cholesky factorization of A

[n,m] = size(A);
x = [];

if n ~= m 
    warning('A is not a square matrix');
    return;
end

if A ~= A'
    warning('A is not a simmetrical matrix');
    return;
end

if abs(A(1:n,1:n))<(sum(A(1:n,:))-A(1:n,1:n))
    warning('A is a singular matrix');
    return;
end

H = A;
H(1,1)=sqrt(A(1,1));

for j=2:n
    for i=1:j-1
    H(i,j) = (A(i,j)- (H(1:i-1,i))'*H(1:i-1,j))/H(i,i);
    end
    H(j,j) = sqrt(A(j,j)-(H(1:j-1,j))'*H(1:j-1,j));
end
H=triu(H);
return