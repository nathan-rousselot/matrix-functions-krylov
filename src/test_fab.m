close all;clear;clc

rng(42);

load('test2.mat')
A = Problem.A;
n = length(A);
b = rand(length(A),1);

tic
fab_mat = expm(A)*b;
time_naive = toc;

tic
fab_perso = fAb(@expm,A,b,22);
time_smart = toc;

err = zeros(n,1);
for k = 1:100
   fab_perso = fAb(@expm,A,b,k);
   err(k) = norm(fab_perso-fab_mat,2)/norm(fab_mat,2);
   disp(k)
end

figure;
semilogy(err)
xlabel('Dimension of the Krylov Space')
ylabel('Relative Error')
grid on