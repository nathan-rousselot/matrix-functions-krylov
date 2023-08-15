function f_vec = generate_diagonal_matrices(n)
% GENERATE_DIAGONAL_MATRICES Generate a list of diagonal matrices with random entries.
%
%   f_vec = GENERATE_DIAGONAL_MATRICES(n) creates a cell array of n 
%   diagonal matrices. The i-th matrix in the list will be of size i x i 
%   with random diagonal entries uniformly disitributed between 0 and 1.
%
%   Input:
%       n       - Number of diagonal matrices to generate. 
%                 An integer where 1 <= n.
%
%   Output:
%       f_vec   - A 1 x n cell array where the i-th cell contains a diagonal 
%                 matrix of size i x i with random diagonal entries.
%
%   Example:
%       matrices = GENERATE_DIAGONAL_MATRICES(3);
%       This will generate a cell array containing 3 matrices:
%       1x1, 2x2, and 3x3, each with random diagonal entries.
%   Author: Nathan Rousselot

    f_vec = cell(1, n);
    for i = 1:n
        f_vec{i} = diag(rand(1, i));
    end
end