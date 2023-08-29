function [V,H,sigma] = shifted_arnoldi_iteration(A,v1,k,sigma)
    n = length(v1);
    [V,H]= deal(zeros(n,n));
    V(:,1) = v1;
    for i = 1:k
        [L,U] = lu(A-sigma*eye(n));
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
        H_temp = H(1:k,1:k);
        sigma = max(eigs(H_temp));
    end
    H = H(1:k,1:k);
    V = V(1:n,1:k);
end