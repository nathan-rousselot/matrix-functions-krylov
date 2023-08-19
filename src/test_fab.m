close all;clear;clc

rng(42);

% load('test2.mat')
A = mmread('C:\Users\nathan.rousselot\Documents\matrix-functions-krylov\Data\bcspwr10.mtx');
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

optimal_k = zeros(3,2);

err = zeros(n,3);
for k = 1:100
   tic
   fab_perso = fAb(@(x)funm(x,@exp),A,b,k);
   time = toc;
   err(k,1) = norm(fab_perso-fab_mat_expm,2)/norm(fab_mat_expm,2);
   if optimal_k(1) == 0 && err(k,1) < 2e-14
        optimal_k(1,1) = k;
        optimal_k(1,2) = time;
   end

   tic
   fab_perso = fAb(@(x)funm(x,@cos),A,b,k);
   time = toc;
   err(k,2) = norm(fab_perso-fab_mat_cosm,2)/norm(fab_mat_cosm,2);
   if optimal_k(2) == 0 && err(k,2) < 2e-14
        optimal_k(2,1) = k;
        optimal_k(2,2) = time;
   end

   tic
   fab_perso = fAb(@(x)funm(x,@sin),A,b,k);
   time = toc;
   err(k,3) = norm(fab_perso-fab_mat_sinm,2)/norm(fab_mat_sinm,2);
   if optimal_k(3) == 0 && err(k,3) < 2e-14
        optimal_k(3,1) = k;
        optimal_k(3,2) = time;
   end
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