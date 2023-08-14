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
    disp(f_A-f_A_test);
    err(i) = norm(f_A-f_A_test);
end

% you can continue with your own tests here...

%%% Generate Dense (random) Matrix

N = 5:2:20;
err_arr = zeros(3,length(N)); % Array of measured errors
for i = 1:length(N)
    A = randn(N(i),N(i)); % Generate gaussian distributed matrix of size N*N ~ N(0,1)
    for j = 1:1
        f_A=f_vec{j}(A);
        f_A_test = matrix_function(f_vec{j},A);
        err_arr(j,i) = norm(f_A-f_A_test);
    end
end

figure;
semilogy(N,err_arr(1,:));
hold on
semilogy(N,err_arr(2,:));
hold on
semilogy(N,err_arr(3,:));
legend('sin(A)','cos(A)','exp(A)')
xlabel('Taille de la matrice N')
ylabel('Erreur')