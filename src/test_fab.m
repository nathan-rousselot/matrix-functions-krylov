close all;clear;clc

rng(42);

load('test2.mat')
A = Problem.A;
n = length(A);
b = rand(length(A),1);

tic
fab_mat_expm = funm(full(A),@exp)*b;
time_naive_expm = toc;

tic
fab_mat_cosm = funm(full(A),@cos)*b;
time_naive_cosm = toc;

tic
fab_mat_sinm = funm(full(A),@sin)*b;
time_naive_sinm = toc;

tic
fab_perso = fAb(@(x)funm(x,@exp),A,b,22);
time_smart = toc;

err = zeros(n,3);
for k = 1:40
   fab_perso = fAb(@(x)funm(x,@exp),A,b,k);
   err(k,1) = norm(fab_perso-fab_mat_expm,2)/norm(fab_mat_expm,2);

   fab_perso = fAb(@(x)funm(x,@cos),A,b,k);
   err(k,2) = norm(fab_perso-fab_mat_cosm,2)/norm(fab_mat_cosm,2);

   fab_perso = fAb(@(x)funm(x,@sin),A,b,k);
   err(k,3) = norm(fab_perso-fab_mat_sinm,2)/norm(fab_mat_sinm,2);
   disp(k)
end

figure;
semilogy(err(:,3),'k--')
hold on
semilogy(err(:,2),'k-*')
hold on
semilogy(err(:,1),'k')
xlabel('Dimension of the Krylov Space')
ylabel('Relative Error')
legend('Sine','Cosine','Exponential')
grid on