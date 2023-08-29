function A_est = shifted_arnoldi_convergence(A,K)
    err = zeros(K,1);
    n = length(A);
    u = randn(n,1);
    u = u/norm(u,2);
    true_lambda = eigs(A);
    [~,H] = arnoldi_iteration(A,u,K);
    sigma = max(eigs(H));
    for k = length(true_lambda):K
        [V,H,sigma] = shifted_arnoldi_iteration(A,u,k,sigma);
        lambda_shifted = eigs(H);
        lambda = sigma+1./lambda_shifted;
        err(k) = norm(lambda-true_lambda,2)/norm(true_lambda,2);
        disp(k)
    end
    close all
    figure(1)
    semilogy(err,'k','LineWidth',0.9)
    xlabel('Dimension of the Krylov space')
    ylabel('Relative Error')
    grid on
    % matlab2tikz('sign_shift_bad.tex')

    figure(2)
    plot(complex(true_lambda),'ko')
    hold on
    plot(complex(lambda),'kx')
    xlabel('Real')
    ylabel('Imaginary')
    legend('Eigenvalues','Ritz Values')
    % matlab2tikz('ritz_shift_bad.tex')
end