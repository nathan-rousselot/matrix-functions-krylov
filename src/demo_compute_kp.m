function demo_compute_kp(A,N,k)
% DEMO_COMPUTE_KP - This function demonstrates the Monte Carlo estimation of the 
% number of positive eigenvalues (k+) of a matrix A. The function compares the 
% estimated k+ against the true value and plots both over multiple iterations.
%
% Parameters:
% - A: Input square matrix (nxn) whose positive eigenvalues are to be estimated.
% - N: Total number of external sampling runs for compute_kp function (default=100).
% - k: Number of Arnoldi iterations (default=20).
%
% Outputs:
% A single figure showing:
%   1. The estimated k+ at each of the N runs (black line).
%   2. The exact number of positive eigenvalues (dashed line).
%
% Usage:
% >> demo_compute_kp(A, N, k)
% >> demo_compute_kp(A, N)       % Uses default k=20
% >> demo_compute_kp(A)          % Uses default N=100, k=20
%
% Assumptions and Limitations:
% - Assumes A to be square; error-checking is not implemented.
% - Assumes availability of compute_kp function and MATLAB's built-in eig and plot functions.
% - Due to stochastic nature, variance in estimates can occur.
%
% Dependencies:
% - compute_kp: Custom function to perform the Monte Carlo estimation of k+.
%
% Author: Nathan Rousselot

    if nargin < 3
        k = 20;
    end
    if nargin < 2
        N=100;
    end
    lambda = eig(full(A));
    kp_exact = length(lambda(lambda>0));
    lb = min(lambda);
    ub = max(lambda);
    kp = zeros(N,1);
    for j = 1:N
        kp(j) = compute_kp(A,k,j,lb,ub);
        disp(j)
    end
    figure;
    plot(kp,'k')
    hold on
    plot([1 N], [kp_exact kp_exact],'k--')
    xlabel('Number of Sampling')
    ylabel('Number of Positive Eigenvalues')
    legend('Estimation','True Value')

end