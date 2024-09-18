clear;
clc;
close all;

addpath('C:\Users\linxi\Documents\MATLAB\chebfun-master\chebfun-master');
addpath('D:\eigtool-master\eigtool-master');

%% Task1

L=15;
N=500;
N1=(N-1)/2;
T=zeros(N);
for i=1:(2*N1+1)
    T(i,i)=-1j*(i-(N1+1));
end
u=zeros(1,31);
for i=-L:L
u(i+L+1)=-2+3*(-1j*i)+(-1j*i)*(-1j*i)+2*(-1j*i)*(-1j*i)*(-1j*i);
end

figure()
plot(real(u),imag(u));
figure()
scatter(real(u),imag(u));
axis([-15 15 -15 15]);


%% Task 2

N = chebop(@(u) -2*u+3*diff(u)+diff(u,2)+2*diff(u,3), [-pi pi], 'periodic');
figure()
[V, D] = eigs(N);
format long, sqrt(-diag(D)) 
plot(V)

%% Task 3

d = [-pi pi]; 
x = chebfun('x', d);
A1 = chebop(@(u) -2*u+3*diff(u)+diff(u,2)+2*diff(u,3), [-pi pi], 'periodic'); 
u0 = sin(x)*sin(x)+0.5*cos(x);
t = [0.1 0.5 2];
u1 = expm(A1, t, u0);
figure()
colr = zeros(6, 3);colr(:,1) = 0.85.^(0:5)';
clf, set(gcf, 'defaultaxescolororder', colr)
plot(chebfun(u1), 'linewidth', 2)
figure()
[V1, D1] = eigs(A1);
format long, sqrt(-diag(D1))  
plot(V1)


%% Task 4

d = [-pi pi]; 
x = chebfun('x', d);
A2 = chebop(@(u) 2*u+3*diff(u)+diff(u,2)+2*diff(u,3), [-pi pi], 'periodic');
u0 = sin(x)*sin(x)+0.5*cos(x);
t = [0.1 0.5 2];
u2 = expm(A2, t, u0);
figure()
colr = zeros(6, 3); colr(:,1) = 0.85.^(0:5)';
clf, set(gcf, 'defaultaxescolororder', colr)
plot(chebfun(u2), 'linewidth', 2);
figure()
[V2, D2] = eigs(A2);
format long, sqrt(-diag(D2))  % integers, to 14 digits
plot(V2)

%% Extra work

A=[0 1 0;0 0 1;1 -1 -0.5];
T=T_fun(L,A);
eigtool(T)

%% Functions

function T=T_fun(n,A)
N=size(A,1);
T=zeros(N*n,N*n);
for i=1:n
    for j=1:n
        if i>j
            T((i-1)*N+1:i*N,(j-1)*N+1:j*N)=A^(i-j-1);
        elseif i==j
            T((i-1)*N+1:i*N,(j-1)*N+1:j*N)=zeros(N);
        end
    end
end
end

