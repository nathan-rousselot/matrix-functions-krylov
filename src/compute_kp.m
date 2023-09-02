function kp = compute_kp(A,k,N,lb,ub)
% COMPUTE_KP - This function estimates the number of positive eigenvalues of a matrix A.
% It employs Monte Carlo sampling along with shifted Arnoldi iterations to approximate 
% the spectral characteristics of the matrix.
%
% Parameters:
% - A: Input square matrix (nxn). A is the matrix whose positive eigenvalues are of interest.
% - k: Number of Arnoldi iterations, indicating the dimensionality of the Krylov subspace.
% - N: Number of Monte Carlo samples.
% - lb: Lower bound for the random shift (sigma).
% - ub: Upper bound for the random shift (sigma).
%
% Outputs:
% - kp: Estimated number of positive eigenvalues of A.
%
% Usage:
% >> kp = compute_kp(A, k, N, lb, ub)
%
% Assumptions and Limitations:
% - Assumes A to be square; error-checking is not implemented.
% - The method is stochastic and the quality of the estimate depends on N and k.
%
% Dependencies:
% - shifted_arnoldi_iteration: Custom function to perform the shifted Arnoldi iteration.
%
% Author: Nathan Rousselot

    q = zeros(N,1);
    n = length(A);
    for i = 1:N
        u = randn(n,1);
        sigma = lb+(ub-lb)*rand(1,1);
        u = u/norm(u,2);
        [~,H,~] = shifted_arnoldi_iteration(A,u,k,sigma);
        q(i) = trace(sign(H));
    end
    kp = mode(q);
    kp = real((kp+k)/2);
end