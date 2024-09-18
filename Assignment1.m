%%  The first method to compute the operator norm ||T||
clear all;
clc;
close all;
% The connection between Toeplitz operators and robust control
A=[0.5 3 2;0 -0.5 -1;0 0 0.2];
b=[1;0;-2];
c=[1;-1;0];
d=0.2;
%check eigenvalues of matrix A
e = eig(A)
%all the eigenvalues locate in the open unit disc: |lambda|<1

%show five dimensional matrix of T
T=T_fun(5,A,b,c,d)

%increase N to estimate ||T||
figure(1)
x=2:1000;
T_approx1=zeros(1,999);
N=1000;
for i=2:N
    T=T_fun(i,A,b,c,d);
%     T_maxgain=norm(conj(T.'));
    T_maxgain=norm(T);
    T_approx1(1,i-1)=T_maxgain;
end
% y1=T_maxgain*ones(1,100);
plot(x,T_approx1);
hold on
plot(x,T_approx1(N-1)*ones(1,N-1));
title('The lower bound as N grows large')
xlabel('N');
ylabel('The lower bound for ||T||');
T_approx1(N-1)
% The lower bound of operator norm ||T|| is 12.6888

%%  A second method to compute the operator norm ||T||

% frequency domain
addpath('C:\Users\linxi\Documents\MATLAB\chebfun-master\chebfun-master');
I=eye(size(A));
cc=c';
F = @(x) abs(d+exp(-1j*x)*cc*inv(I-exp(-1j*x)*A)*b);
f = chebfun(F,[0,2*pi]);
figure(2)
plot(f)
r = roots(diff(f));
hold on, plot(r,f(r),'.r'), grid on
max(f)
title('The uppper bound in frequency domain')
xlabel('\omega');
ylabel('The uppper bound for ||T||');

% use MATLAB norm function to verify
Ts=0.01;%choosing the sampling time as 0.01s
sys = ss(A,b,c',d,Ts);
[ninf,fpeak]  = norm(sys,Inf);
ninf
% The upper bound of operator norm ||T|| is 12.6889

%% functions
function T=T_fun(n,A,b,c,d)
cc=zeros(1,n);
cc(1)=d;
rr=zeros(1,n);
rr(1)=d;
rr(2)=1;
 for i=3:n
   rr(i)=c'*( A^(i-2))*b;
 end
 T=toeplitz(rr,cc);
end

% function H=H_fre(omega)
% A=[0.5 3 2;0 -0.5 -1;0 0 0.2];
% b=[1;0;-2];
% c=[1;-1;0];
% d=0.2;
% I=eye(size(A));
% H=abs(d+exp(-1j*omega)*c'*inv(I-exp(I-1j*omega)*A)*b);
% end

