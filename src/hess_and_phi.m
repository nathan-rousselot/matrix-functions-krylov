function [H,lambda,m,Q,J]=hess_and_phi(A)
    [Q,H]=hess(full(A));
    J=jordan(H);
    [lambda,m]=minimal_polynomial_J(J);
end