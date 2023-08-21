function [V,H] = arnoldi_iteration(A,v1,k)
    n = length(v1);
    [V,H]= deal(zeros(n,n));
    V(:,1) = v1;
    for i = 1:k
        w = A*V(:,i);
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