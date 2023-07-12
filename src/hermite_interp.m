function P = hermite_interp(X, FdF)
%   Hermite Interpolation by divided differences.
%   
%   X should be a row vector
%   FdF contains function F and its derivatives, corresponding to the
%   points X
%   Columns of FdF correspond to points in X

    [ ~, m ] = size(X); % polynomial order <= m-1
    p = zeros(m, m);
    FdF = FdF';
    p(:, 1) = FdF(:,1);
    for j = 2 : m
        for i = 1 : (m - j + 1)
            if ((X(i + j - 1) - X(i))==0)
                p(i,j) = FdF(i,j-1)/(factorial(j-1));
            else
            p(i,j) = (p(i + 1, j - 1) - p(i, j - 1)) / (X(i + j - 1) - X(i));
            end
        end
    end

    P = @(x) eval_diff(x,X,p(1,:)); %evaluating the divided differences requires care, off-loaded to separate function

end