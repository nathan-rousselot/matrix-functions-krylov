function f = fAb(fun,A,b,k)
    q1 = b/norm(b);
    [Q,H] = arn_it(A,q1,k);
    H = H(1:k,1:k);
    Q = Q(:,1:k);
    e = zeros(k,1);
    e(1) = 1;
    f = norm(b)*Q*fun(H)*e;
end

function [Q,H] = arn_it(A,q1,k)
%ARNOLDI    Arnoldi iteration
%   [Q,H] = ARNOLDI(A,q1,M) carries out M iterations of the
%   Arnoldi iteration with N-by-N matrix A and starting vector q1
%   (which need not have unit 2-norm).  For M < N it produces
%   an N-by-(M+1) matrix Q with orthonormal columns and an
%   (M+1)-by-M upper Hessenberg matrix H such that
%   A*Q(:,1:M) = Q(:,1:M)*H(1:M,1:M) + H(M+1,M)*Q(:,M+1)*E_M',
%   where E_M is the M'th column of the M-by-M identity matrix.
    
    n = length(A);
    if nargin < 3, k = n; end
    q1 = q1/norm(q1);
    Q = zeros(n,k); Q(:,1) = q1;
    H = zeros(min(k+1,k),n);
    
    for i=1:k
        z = A*Q(:,i);
        for j=1:i
            H(j,i) = Q(:,j)'*z;
            z = z - H(j,i)*Q(:,j);
        end
        if i < n
           H(i+1,i) = norm(z);
           if H(i+1,i) == 0, return, end
           Q(:,i+1) = z/H(i+1,i);
       end
    end
end