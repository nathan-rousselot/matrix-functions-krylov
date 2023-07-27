function F = matrix_function(fun,A)
    [~,lambda,m,~,~]=hess_and_phi(A);
    FdF = zeros(max(m)+1,length(lambda));
    fun = sym(fun);
    for i = 1:max(max(m)+1)
        for j = 1:length(lambda)
            FdF(i,j)=subs(fun,lambda(j));
        end
        fun = diff(fun);
    end
    P = hermite_interp(lambda',FdF);
    F = P(A);
end
