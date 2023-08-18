function F = matrix_function(fun,A)
% MATRIX_FUNCTION: Computes a matrix function of a given matrix using the Hermite polynomial method.
%
% Input:
% - fun : The function handle that should be applied to the matrix A.
%         For example, to compute the matrix exponential, fun could be @(x) exp(x).
% - A   : Square matrix whose function value needs to be computed.
%
% Output:
% - F   : Matrix that represents the function of the input matrix A. That is, F = fun(A).
%
% Dependencies:
% This function uses `hess_and_phi` and `hermitePoly` functions, 
% which are provided with this code.
%
% Description:
% 1. Decomposes the matrix A using the `hess_and_phi` function to get its
% eigenvalues, and the minimal polynomial multiplicative coefficients
% 2. Prepares a matrix with the results of the function applied to the eigenvalues and their derivatives.
% 3. Constructs a Hermite polynomial using the prepared matrix.
% 4. Evaluates the matrix function using the Hermite polynomial and the input matrix A.
%
% Note:
% Ensure the provided function fun can be symbolically differentiated 
% and that all required dependencies are available.
%
% Example:
% F = matrix_function(@(x) sin(x), [1, 2; 3, 4]);
%
% This will compute the matrix sine of the matrix [1, 2; 3, 4] and store the result in F.
%
% Author: Nathan Rousselot

    [~,lambda,m,~,~]=hess_and_phi(A);
    M = max(m)+1;
    lambda = lambda';
    fun = sym(fun);
    k = length(lambda);
    FdF = zeros(M,k);
    for i = 1:M
        for j = 1:k
            FdF(i,j)=subs(fun,lambda(j));
        end
        fun = diff(fun);
    end
    FdF = [lambda; FdF];
    [P, ~, ~ ] = hermitePoly(FdF);
    F = polyvalm(P,A);
end