function kp = compute_kp(A,k,N)
    q = zeros(N,1);
    n = length(A);
    for i = 1:N
        u = randn(n,1);
        u = u/norm(u,2);
        [~,H] = arnoldi_iteration(A,u,k);
        q(i) = trace((sign(H)));
    end
    kp = mode(q);
    kp = (kp+k)/2;
end