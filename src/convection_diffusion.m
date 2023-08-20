clear;
close;
clc;
n=50;
x=linspace(0,1,n);
y=linspace(0,1,n);
dx = x(2)-x(1);
e = ones(n,1);
L = spdiags([e -2*e e], -1:1, n, n)/(dx*dx);
L = kron(speye(n),L)+kron(L,speye(n));
D = spdiags([e -e], [-1,1], n, n)/(2*dx);
alpha = [.5,1];
D = alpha(1)*kron(speye(n),D)+alpha(2)*kron(D,speye(n));
epp = 1e-1;
A = epp*L + D;
Pe = normest(alpha)*dx/(2*epp);

U0 = ((16*((1-y).*y.^2).^2).')*(16*((1-x).*x.^2).^2);
u0=U0(:);
%%
tic
u_exact=zeros(n*n,10);
u_exact(:,1)=u0;
dt = .1;
for i = 1:9
u_exact(:,i+1) = expm(dt*A)*u_exact(:,i);
end
t_exact = toc;


err = zeros(40,1);
for k =1:100
    u=zeros(n*n,10);
    u(:,1)=u0;
    for i = 1:9
       u(:,i+1) = fAb(@expm,dt*A,u(:,i),k);
    end
    err(k) = norm(u-u_exact,2)/norm(u_exact,2);
    disp(k)
end

tic
u=zeros(n*n,10);
u(:,1)=u0;
for i = 1:9
u(:,i+1) = fAb(@expm,dt*A,u(:,i),27);
end
t_optim = toc;
%%
figure(1)
clf
for i = 1:9
subplot(3,3,i)
U = reshape(u(:,i+1),[n,n]);
contourf(U,10)
axis equal
end

figure(2)
clf
for i = 1:9
subplot(3,3,i)
U = reshape(u_exact(:,i+1),[n,n]);
contourf(U,10)
axis equal
end

figure(3)
semilogy(err,'k','LineWidth',0.9)
xlabel('Dimension of the Krylov space')
ylabel('Relative error')
grid on
