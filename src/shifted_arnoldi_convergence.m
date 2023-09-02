function shifted_arnoldi_convergence(A,K)
% SHIFTED_ARNOLDI_CONVERGENCE - This function examines the convergence behavior of 
% shifted Arnoldi iteration for approximating the eigenvalues of a matrix A.
% It plots the relative error between the Ritz values and the actual eigenvalues 
% as a function of the dimension of the Krylov space.
%
% Parameters:
% - A: Input square matrix (nxn). A is the matrix whose eigenvalues are of interest.
% - K: Maximum dimension of the Krylov space.
%
% Outputs:
% Two figures are generated:
%   Figure 1: A semi-logarithmic plot of the relative error against the dimension
%             of the Krylov space.
%   Figure 2: A plot comparing the actual eigenvalues (circles) and Ritz values (crosses)
%             in the complex plane.
%
% Usage:
% >> shifted_arnoldi_convergence(A, K)
%
% Assumptions and Limitations:
% - Assumes A to be square; error checking not implemented.
%
% Dependencies:
% - shifted_arnoldi_iteration: Custom function to perform the shifted Arnoldi iteration.
%
% Author: Nathan Rousselot

    err = zeros(K,1);
    n = length(A);
    u = randn(n,1);
    u = u/norm(u,2);
    true_lambda = eigs(A);
    sigma = (true_lambda(1)+true_lambda(2))/2;
    for k = 1:K
        [~,H,sigma] = shifted_arnoldi_iteration(A,u,k,sigma);
        mu = eigs(H);
        lambda = sigma+1./mu;
        err(k) = norm(lambda-true_lambda(1:length(lambda)),2)/norm(true_lambda(1:length(lambda)),2);
        disp(k)
    end
    close all
    figure(1)
    semilogy(err,'k','LineWidth',0.9)
    xlabel('Dimension of the Krylov space')
    ylabel('Relative Error')
    grid on

    figure(2)
    plot(complex(true_lambda),'ko')
    hold on
    plot(complex(lambda),'kx')
    xlabel('Real')
    ylabel('Imaginary')
    legend('Eigenvalues','Ritz Values')
end