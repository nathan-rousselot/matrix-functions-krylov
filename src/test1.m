%% Test 1
% Here you should write your test routines, to verify that your
% implementation of the matrix function evaluation is correct
% You should test it on the matrix 'test1.mat'

% A first test is the following:
clear
clc
load test1.mat

f_vec = {@(x)funm(x,@sin),@(x)funm(x,@cos),@(x)funm(x,@exp) };
F= f_vec{3};%@(x)funm(x,@exp);
err = zeros(3,1);
for i = 1:3
    f_A=f_vec{i}(A);
    f_A_test = %your implementation goes here
    err(i) = norm(f_A-f_A_test);
end

% you can continue with your own tests here...