function [Q,R] = qrgram(A)

% QR factorization by using the modified gram-schmidt orthogonalization algo
% Q has orthogonal columns, so if r(Q)=n then Q's columns form an orthogonal basis for col(A).
% That's why QR factorization could also be seen as way for calculating an 
% orthogonal basis of a set of vectors.

[m,n]=size(A);
Q=zeros(m,n); Q(:,1) = A(:,1); 
R=zeros(n); R(1,1)=1;

for k=1:n
    R(k,k) = norm(A(:,k));
    Q(:,k) = A(:,k)/R(k,k);
    R(k,k+1:n) = Q(:,k)'*A(:,k+1:n);
    A(:,k+1:n) = A(:,k+1:n)-Q(:,k)*R(k,k+1:n);
end

return