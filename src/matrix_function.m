function F = matrix_function(fun,A)
    [~,lambda,~,~,~]=hess_and_phi(A);
    h = 1e-6;
    %lambda = uniquetol(lambda,h)';
    lambda = lambda';
    fun = sym(fun);
    k = length(lambda);
    FdF = zeros(k,k);

    for i = 1:k
        for j = 1:k
            FdF(j,i)=subs(fun,lambda(j));
        end
        fun = diff(fun);
        %fun = @(x)(fun(x+h)-fun(x-h))/(2*h); 
    end
    P = hermite_interp(lambda,FdF);
    P = @(x) arrayfun(@(z) P(z), x);
    F = P(A);
end
