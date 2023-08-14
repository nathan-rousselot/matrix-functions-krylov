function F = matrix_function(fun,A)
    [~,lambda,m,~,~]=hess_and_phi(A);
    M = max(m)+1;
    h = 1e-6;
    lambda = lambda';
    fun = sym(fun);
    k = length(lambda);
    FdF = zeros(k,M);
    for i = 1:M
        for j = 1:k
            FdF(j,i)=subs(fun,lambda(j));
        end
        fun = diff(fun);
    end
    FdF = [lambda; FdF'];
    [P, ~, ~ ] = hermitePoly(FdF);
    F = polyvalm(P,A);
end
