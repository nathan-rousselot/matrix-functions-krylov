function [V,H,sigma] = shifted_arnoldi_iteration(A,v1,k,sigma)
% SHIFTED_ARNOLDI_ITERATION - This function performs a shifted Arnoldi iteration
% to construct an orthonormal basis for the Krylov subspace spanned by
% {v1, (A-sigma*I)^{-1}*v1, ..., (A-sigma*I)^{-(k-1)}*v1}.
% The function iteratively builds matrices V and H such that (A - sigma*I)*V = V*H + h_{k+1,k}*v_{k+1}*e_k',
% where H is upper Hessenberg and V has orthonormal columns.
%
% Mathematical Background:
% Shifted Arnoldi iteration is a variant of standard Arnoldi iteration used
% for accelerated convergence in eigenvalue problems, particularly those close to
% a shift point 'sigma'. The function leverages LU decomposition for solving 
% the shifted system efficiently.
%
% Parameters:
% - A: Input square matrix (nxn). A is the matrix for which the Krylov
%      subspace is being approximated.
% - v1: Initial vector (nx1). v1 should be non-zero for a meaningful Krylov
%       subspace.
% - k: Number of Arnoldi iterations to be performed. Determines the dimension
%      of the Krylov subspace.
% - sigma: Scalar shift applied to the matrix A.
%
% Outputs:
% - V: Orthonormal basis matrix (nxk). Columns of V form an orthonormal basis
%      for the Krylov subspace.
% - H: Upper Hessenberg matrix (kxk). H represents the matrix A in the shifted
%       Krylov subspace.
% - sigma: The shift parameter used in the iteration.
%
% Usage:
% >> [V, H, sigma] = shifted_arnoldi_iteration(A, v1, k, sigma)
%
% Algorithm Complexity:
% The function has a computational complexity of O(n*k^2), similar to
% the standard Arnoldi, with an additional overhead for LU decomposition,
% which is O(n^3). The LU decomposition, however, is only performed once.
%
% Assumptions and Limitations:
% - The function assumes A to be square and v1 to be non-zero; error checking
%   is not implemented.
% - If h_{k+1,k} == 0, the Krylov subspace has been fully spanned before reaching
%   k iterations, and the function returns early.
%
% Author: Nathan Rousselot
%
% References:
% [1] "Matrix Computations" by Golub and Van Loan

    n = length(v1);
    [V,H]= deal(zeros(n,n));
    V(:,1) = v1;
    [L, U] = lu(A - sigma * eye(n));
    for i = 1:k
        w = U\(L\V(:,i));
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