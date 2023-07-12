function [lambda,m] = minimal_polynomial_J(J)
    % Input: J is a square matrix in Jordan normal form
    
    lambda = unique(diag(J));
    m = arrayfun(@(lambda) max_jordan_block_size(J, lambda), lambda);
end

function block_size = max_jordan_block_size(J, lambda)
    % Input: J is a square matrix in Jordan normal form
    %        lambda is an eigenvalue of J
    
    % Find the locations of lambda in the diagonal of J
    locations = find(diag(J) == lambda);

    % Calculate the size of the largest Jordan block
    block_size = 1;
    for i = 1:numel(locations)
        j = locations(i);
        current_block_size = 1;
        while j < size(J, 2) && J(j, j + 1) ~= 0
            current_block_size = current_block_size + 1;
            j = j + 1;
        end
        block_size = max(block_size, current_block_size);
    end
end