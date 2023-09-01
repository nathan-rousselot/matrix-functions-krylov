function kp = compute_kp(A,k,N,lb,ub)
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