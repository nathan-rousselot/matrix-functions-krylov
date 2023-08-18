function f = fAb(fun,A,b,k)
% fAb - Computes the action of a matrix function on a vector using the Arnoldi method.
%
% Input:
%   fun  : A handle to the matrix function to be evaluated.
%          For example, to compute the matrix sine, fun could be @(x) sinm(x).
%   A    : A square matrix (NxN) to which the function should be applied.
%   b    : A column vector to which the matrix function is applied.
%   k    : Number of Arnoldi iterations (typically much smaller than N).
%
% Output:
%   f    : Result of fun(A)*b, approximately computed using k Arnoldi iterations.
%
% Description:
% 1. Normalizes vector b and sets it as the initial vector for Arnoldi iteration.
% 2. Constructs a Krylov subspace using the Arnoldi iteration.
% 3. Computes the action of the matrix function on the vector using the constructed Krylov subspace.
%
% Note:
% - The function operates by expanding the matrix-vector product fun(A)*b into a Krylov subspace, 
%   allowing for approximate but faster computations. The accuracy improves with increased k.
% - It is significantly faster than using Matlab's naive way: funm(A)*b (if
%   k is chose appropriately)
%
% Example:
% f = fAb(@(x) expm(x), [1, 2; 3, 4], [1; 0], 2);
%
% This computes the action of the matrix exponential on the vector [1; 0] 
% using 2 Arnoldi iterations and stores the result in f. 
%
% Author: Nathan Rousselot
    v1 = b/norm(b,2);
    [V,H] = arnoldi_iteration(A,v1,k);
    e = zeros(k,1);
    e(1) = 1;
    f = norm(b,2)*V*fun(H)*e;
end

function [V,H] = arnoldi_iteration(A,v1,k)
    n = length(v1);
    [V,H]= deal(zeros(n,n));
    V(:,1) = v1;
    for i = 1:k
        w = A*V(:,i);
        for j = 1:i
            H(j,i) = V(:,j)'*w;
            w = w - H(j,i)*V(:,j);
        end
        H(i+1,i) = norm(w,2);
        if H(i+1,i) == 0
            H = H(1:k,1:k);
            V = V(1:n,1:k);
            return;
        end
        V(:,i+1) = w/H(i+1,i);
    end
    H = H(1:k,1:k);
    V = V(1:n,1:k);
end