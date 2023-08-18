% Test 1
% Here you should write your test routines, to verify that your
% implementation of the matrix function evaluation is correct
% You should test it on the matrix 'test1.mat'

% A first test is the following:
clear
clc
load Matrix1.mat

f_vec = {@(x)funm(x,@sin),@(x)funm(x,@cos),@(x)funm(x,@exp) };
F= f_vec{3};%@(x)funm(x,@exp);
err = zeros(3,1);
for i = 1:3
    f_A=f_vec{i}(A);
    f_A_test = matrix_function(f_vec{i},A);
    err(i) = norm(f_A-f_A_test,2);
end

% you can continue with your own tests here...

% random seed for reproducibility
rng(42);

n = 15;

%%% DIAGONAL MATRICES
test_matrices = generate_diagonal_matrices(n);
err_diag = zeros(length(test_matrices),length(f_vec));
for j = 1:length(test_matrices)
    for i = 1:length(f_vec)
        f_A=f_vec{i}(test_matrices{j});
        f_A_test = matrix_function(f_vec{i},test_matrices{j});
        err_diag(j,i) = norm(f_A-f_A_test,2)/norm(f_A,2);
    end
end

%%% SYMMETRIC MATRICES
test_matrices = generate_symmetric_matrices(n);
err_sym = zeros(length(test_matrices),length(f_vec));
for j = 1:length(test_matrices)
    for i = 1:length(f_vec)
        f_A=f_vec{i}(test_matrices{j});
        f_A_test = matrix_function(f_vec{i},test_matrices{j});
        err_sym(j,i) = norm(f_A-f_A_test,2)/norm(f_A,2);
    end
    disp(j)
end

%%% DENSE MATRICES

test_matrices = cell(1, n);
for i = 1:n
    test_matrices{i} = rand(i,i);
end
err_dens = zeros(length(test_matrices),length(f_vec));
for j = 1:length(test_matrices)
    for i = 1:length(f_vec)
        f_A=f_vec{i}(test_matrices{j});
        f_A_test = matrix_function(f_vec{i},test_matrices{j});
        err_dens(j,i) = norm(f_A-f_A_test,2)/norm(f_A,2);
    end
    disp(j)
end

figure;
semilogy(err_diag(:,1),'k--')
hold on
semilogy(err_sym(:,1),'k-*')
hold on
semilogy(err_dens(:,1),'k')
legend('diagonal','symmetric','dense')
xlabel('Size of the matrices')
ylabel('Relative Error')
title('Cosine function')
grid on

figure;
semilogy(err_diag(:,2),'k--')
hold on
semilogy(err_sym(:,2),'k-*')
hold on
semilogy(err_dens(:,2),'k')
legend('diagonal','symmetric','dense')
xlabel('Size of the matrices')
ylabel('Relative Error')
title('Cosine function')
grid on

figure;
semilogy(err_diag(:,3),'k--')
hold on
semilogy(err_sym(:,3),'k-*')
hold on
semilogy(err_dens(:,3),'k')
legend('diagonal','symmetric','dense')
xlabel('Size of the matrices')
ylabel('Relative Error')
title('Exponenial function')
grid on