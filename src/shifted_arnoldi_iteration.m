function [V,H,sigma] = shifted_arnoldi_iteration(A,v1,k,sigma)
    n = length(v1);
    % tic
    % [~,non_shifted] = arnoldi_iteration(A,v1,length(n/5));
    % sigma = median(eigs(non_shifted));
    % disp('Time to compute sigma:')
    % toc
    [V,H]= deal(zeros(n,n));
    V(:,1) = v1;
    [L, U] = lu(A - sigma * eye(n));
    for i = 1:k
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
    end
    H = H(1:k,1:k);
    V = V(1:n,1:k);
end