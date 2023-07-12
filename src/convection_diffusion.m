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
epp = 1e-2;
A = epp*L + D;
Pe = normest(alpha)*dx/(2*epp);

U0 = ((16*((1-y).*y.^2).^2).')*(16*((1-x).*x.^2).^2);
u0=U0(:);
%%
u_exact=zeros(n*n,10);
u_exact(:,1)=u0;
dt = .1;
for i = 1:9
u_exact(:,i+1) = expm(dt*A)*u_exact(:,i);
end

u=zeros(n*n,10);
u(:,1)=u0;
for i = 1:9
[H,Q] = arnoldi(A,u(:,i),20);
eH = expm(dt*H);
u(:,i+1) = Q*(eH*((Q'*u(:,i))));
end

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