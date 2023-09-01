function demo_compute_kp(A,N,k)
    if nargin < 3
        k = 20;
    end
    if nargin < 2
        N=100;
    end
    % rng(42);
    lambda = eig(full(A));
    % sigma = (lambda(1)+lambda(2))/2;
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