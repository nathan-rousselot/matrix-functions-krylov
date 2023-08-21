function arnoldi_convergence(A,K)
    err = zeros(K,1);
    n = length(A);
    u = randn(n,1);
    u = u/norm(u,2);
    for k = 99:K
        [V,H] = arnoldi_iteration(A,u,k);
        A_est = V*H*V';
        err(k) = norm(A_est-A,2)/norm(full(A),2);
        disp(k)
    end
    figure(1)
    semilogy(err,'k','LineWidth',0.9)
    xlabel('Dimension of the Krylov space')
    ylabel('Relative Error')
    grid on
end