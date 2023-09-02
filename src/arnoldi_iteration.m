function [V,H] = arnoldi_iteration(A,v1,k)
% ARNOLDI_ITERATION - This function performs Arnoldi iteration to construct
% an orthonormal basis for the Krylov subspace spanned by {v1, A*v1, ..., A^(k-1)*v1}.
% Given a square matrix A and an initial vector v1, it iteratively builds
% matrices V and H such that A*V = V*H + h_{k+1,k}*v_{k+1}*e_k', where H is upper
% Hessenberg and V has orthonormal columns.
%
% Mathematical Background:
% Arnoldi iteration is an important algorithm in numerical linear algebra for
% constructing an orthonormal basis of the Krylov subspace. It can be employed
% in eigenvalue problems, solving linear systems, and model order reduction.
%
% Parameters:
% - A: Input square matrix (nxn). A is the matrix for which we are interested
%      in generating the Krylov subspace.
% - v1: Initial vector (nx1). v1 should be non-zero and of unit l2-norm
% - k: Number of Arnoldi iterations to be performed. Determines the dimension
%      of the Krylov subspace (spanned by {v1, A*v1, ..., A^(k-1)*v1}).
%
% Outputs:
% - V: Orthonormal basis matrix (nxk). Columns of V form an orthonormal basis
%      for the Krylov subspace.
% - H: Upper Hessenberg matrix (kxk). H represents the matrix A in the
%       Krylov subspace and satisfies the Arnoldi relation A*V = V*H + h_{k+1,k}*v_{k+1}*e_k'.
%
% Usage:
% >> [V, H] = arnoldi_iteration(A, v1, k)
%
% Algorithm Complexity:
% The function has a computational complexity of O(n*k^2), dominated by the
% Gram-Schmidt process for orthogonalizing the vector w against the columns of V.
%
% Assumptions and Limitations:
% - The function assumes A to be square and v1 to be non-zero; error
%   checking is not implemented.
% - If the algorithm encounters h_{k+1,k} == 0, it means the Krylov subspace
%   has been fully spanned before reaching k iterations and the function returns
%   early.
%
% Author: Nathan Rousselot
%
% References:
% [1] "Matrix Computations" by Golub and Van Loan

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