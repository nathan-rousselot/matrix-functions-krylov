function p = construct_minimal_polynomial(lambda, m)
    p = 1;
    for i = 1:length(lambda)
        factor = [1, -lambda(i)];
        for j = 2:m(i)
            factor = conv(factor, [1, -lambda(i)]);
        end
        p = conv(p, factor);
    end
end